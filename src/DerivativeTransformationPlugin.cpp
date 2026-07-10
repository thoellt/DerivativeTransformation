#include "DerivativeTransformationPlugin.h"

#include <PointData/PointData.h>

#include <actions/DecimalAction.h>
#include <actions/IntegralAction.h>

#include <QDebug>
#include <QDialog>
#include <QDialogButtonBox>
#include <QEventLoop>
#include <QFutureWatcher>
#include <QGridLayout>
#include <QThread>
#include <QTimer>
#include <QtConcurrent/QtConcurrent>

#include <algorithm>
#include <atomic>
#include <utility>
#include <vector>

Q_PLUGIN_METADATA(IID "studio.manivault.DerivativeTransformationPlugin")

using namespace mv;

using Kernel = DerivativeTransformationPlugin::Kernel;

const QMap<Kernel, QString> DerivativeTransformationPlugin::kernels = QMap<Kernel, QString>({
    { Kernel::Forward, "Forward differences" },
    { Kernel::Central, "Central differences" },
    { Kernel::Central5, "Central differences (5-point)" },
    { Kernel::SavitzkyGolay, "Savitzky-Golay" },
    { Kernel::Gaussian, "Gaussian derivative" }
});

namespace {

/**
 * Modal parameter dialog for the Savitzky-Golay and Gaussian kernels.
 * Returns false if the user cancelled.
 */
bool promptKernelParameters(Kernel kernel, int& sgWindowSize, int& sgPolynomialOrder, float& sigma)
{
    QDialog dialog;
    dialog.setWindowTitle(QString("Derivative Transformation: %1").arg(DerivativeTransformationPlugin::getKernelName(kernel)));

    auto layout = new QGridLayout(&dialog);

    mv::gui::IntegralAction* windowSizeAction = nullptr;
    mv::gui::IntegralAction* polynomialOrderAction = nullptr;
    mv::gui::DecimalAction* sigmaAction = nullptr;

    int row = 0;

    const auto addAction = [&](mv::gui::WidgetAction* action) {
        layout->addWidget(action->createLabelWidget(&dialog), row, 0);
        layout->addWidget(action->createWidget(&dialog), row, 1);
        ++row;
    };

    if (kernel == Kernel::SavitzkyGolay) {
        windowSizeAction = new mv::gui::IntegralAction(&dialog, "Window size", 3, 101, sgWindowSize);
        windowSizeAction->setToolTip("Sliding window size in samples (odd values; even input is rounded up)");
        addAction(windowSizeAction);

        polynomialOrderAction = new mv::gui::IntegralAction(&dialog, "Polynomial order", 1, 10, sgPolynomialOrder);
        polynomialOrderAction->setToolTip("Order of the fitted polynomial (must be smaller than the window size)");
        addAction(polynomialOrderAction);
    }

    if (kernel == Kernel::Gaussian) {
        sigmaAction = new mv::gui::DecimalAction(&dialog, "Sigma", 0.1f, 25.0f, sigma, 2);
        sigmaAction->setToolTip("Standard deviation of the Gaussian in samples");
        addAction(sigmaAction);
    }

    auto buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
    QObject::connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    layout->addWidget(buttons, row, 0, 1, 2);

    if (dialog.exec() != QDialog::Accepted)
        return false;

    if (windowSizeAction) {
        sgWindowSize = windowSizeAction->getValue() | 1; // round up to odd
        sgPolynomialOrder = std::min(polynomialOrderAction->getValue(), sgWindowSize - 1);
    }

    if (sigmaAction)
        sigma = sigmaAction->getValue();

    return true;
}

} // namespace

DerivativeTransformationPlugin::DerivativeTransformationPlugin(const PluginFactory* factory) :
    TransformationPlugin(factory),
    _kernel(Kernel::Central),
    _sgWindowSize(7),
    _sgPolynomialOrder(2),
    _sigma(1.0f)
{
}

void DerivativeTransformationPlugin::transform()
{
    auto points = getInputDataset<Points>();

    if (!points.isValid())
        return;

    const auto numDimensions = static_cast<int>(points->getNumDimensions());
    const auto numPoints = points->getNumPoints();

    if (numDimensions < 2) {
        qWarning() << "Derivative Transformation: input dataset needs at least 2 dimensions (spectral samples), got" << numDimensions;
        return;
    }

    auto kernelDescription = getKernelName(_kernel);

    if (_kernel == Kernel::SavitzkyGolay)
        kernelDescription += QString(" w%1 p%2").arg(_sgWindowSize).arg(_sgPolynomialOrder);

    if (_kernel == Kernel::Gaussian)
        kernelDescription += QString(" sigma %1").arg(_sigma);

    auto& task = points->getTask();

    task.setName("Derivative Transformation");
    task.setRunning();
    task.setProgressDescription(QString("1st derivative, %1").arg(kernelDescription));

    const auto weightTable = derivative::buildWeightTable(_kernel, numDimensions, _sgWindowSize, _sgPolynomialOrder, _sigma);

    // Gather the input into a flat point-major buffer on the GUI thread
    // (dataset access is GUI-thread only). visitData is subset-aware; do NOT
    // use constVisitFromBeginToEnd here — it walks the raw source buffer,
    // which is larger than numPoints * numDimensions when the input is a
    // subset (buffer overflow).
    std::vector<float> input(static_cast<std::size_t>(numPoints) * numDimensions);

    points->visitData([&](auto pointData) {
        std::size_t i = 0;
        for (const auto point : pointData)
            for (int d = 0; d < numDimensions; ++d)
                input[i++] = static_cast<float>(point[d]);
    });

    std::vector<float> derivativeValues(input.size());

    // Parallelize over points: workers only touch the flat buffers, never the
    // dataset or any UI. Chunks are sized for ~4 chunks per core so progress
    // stays granular without scheduling overhead.
    const auto threadCount = std::max(1, QThread::idealThreadCount());
    const auto chunkSize = std::max<std::uint32_t>(1024, (numPoints + threadCount * 4 - 1) / (threadCount * 4));

    std::vector<std::pair<std::uint32_t, std::uint32_t>> chunks;
    for (std::uint32_t begin = 0; begin < numPoints; begin += chunkSize)
        chunks.emplace_back(begin, std::min<std::uint32_t>(numPoints, begin + chunkSize));

    std::atomic<std::uint32_t> pointsProcessed{ 0 };

    auto future = QtConcurrent::map(chunks, [&](const std::pair<std::uint32_t, std::uint32_t>& chunk) {
        derivative::convolveRange(input.data(), derivativeValues.data(), chunk.first, chunk.second, numDimensions, weightTable);
        pointsProcessed.fetch_add(chunk.second - chunk.first, std::memory_order_relaxed);
    });

    // Wait on the GUI thread, reporting progress from the atomic counter
    {
        QEventLoop waitLoop;

        QFutureWatcher<void> watcher;
        QObject::connect(&watcher, &QFutureWatcher<void>::finished, &waitLoop, &QEventLoop::quit);
        watcher.setFuture(future);

        QTimer progressTimer;
        progressTimer.setInterval(50);
        QObject::connect(&progressTimer, &QTimer::timeout, &progressTimer, [&]() {
            task.setProgress(static_cast<float>(pointsProcessed.load(std::memory_order_relaxed)) / numPoints);
        });
        progressTimer.start();

        if (!future.isFinished())
            waitLoop.exec();
    }

    input.clear();
    input.shrink_to_fit();

    auto derived = mv::data().createDerivedDataset<Points>(
        points->getGuiName() + QString(" (1st derivative, %1)").arg(kernelDescription), points);

    derived->setData(std::move(derivativeValues), numDimensions);

    auto dimensionNames = points->getDimensionNames();

    if (static_cast<int>(dimensionNames.size()) == numDimensions) {
        for (auto& name : dimensionNames)
            name = QString("d/dλ %1").arg(name);
    } else {
        dimensionNames.resize(numDimensions);
        for (int d = 0; d < numDimensions; ++d)
            dimensionNames[d] = QString("d/dλ %1").arg(d);
    }

    derived->setDimensionNames(dimensionNames);

    events().notifyDatasetDataChanged(derived);

    task.setFinished();
}

// -----------------------------------------------------------------------------
// Factory
// -----------------------------------------------------------------------------

DerivativeTransformationPluginFactory::DerivativeTransformationPluginFactory()
{
    getPluginMetadata().setDescription("First derivative of spectral response functions");
    getPluginMetadata().setSummary("Computes the first derivative of per-point spectra via selectable derivative kernels (finite differences, Savitzky-Golay, Gaussian).");
    getPluginMetadata().setCopyrightHolder({ "REPLACE ME" });
    getPluginMetadata().setAuthors({ { "REPLACE ME", { "Developer" }, { "REPLACE ME org" } } });
    getPluginMetadata().setLicenseText("REPLACE ME (e.g. LGPL v3.0)");
}

DerivativeTransformationPlugin* DerivativeTransformationPluginFactory::produce()
{
    return new DerivativeTransformationPlugin(this);
}

mv::DataTypes DerivativeTransformationPluginFactory::supportedDataTypes() const
{
    return { PointType };
}

mv::gui::PluginTriggerActions DerivativeTransformationPluginFactory::getPluginTriggerActions(const mv::Datasets& datasets) const
{
    mv::gui::PluginTriggerActions pluginTriggerActions;

    if (datasets.count() >= 1 && PluginFactory::areAllDatasetsOfTheSameType(datasets, PointType)) {
        const auto addTriggerAction = [this, &pluginTriggerActions, datasets](const Kernel& kernel) {
            const auto kernelName = DerivativeTransformationPlugin::getKernelName(kernel);
            const auto isParameterized = kernel == Kernel::SavitzkyGolay || kernel == Kernel::Gaussian;
            const auto menuName = isParameterized ? kernelName + "..." : kernelName;

            auto pluginTriggerAction = new mv::gui::PluginTriggerAction(const_cast<DerivativeTransformationPluginFactory*>(this), this,
                QString("Derivative Transformation/%1").arg(menuName), QString("Compute the first derivative (%1)").arg(kernelName), icon(),
                [this, datasets, kernel, isParameterized](mv::gui::PluginTriggerAction& pluginTriggerAction) -> void {
                    int sgWindowSize = 7;
                    int sgPolynomialOrder = 2;
                    float sigma = 1.0f;

                    if (isParameterized && !promptKernelParameters(kernel, sgWindowSize, sgPolynomialOrder, sigma))
                        return;

                    for (const auto& dataset : datasets) {
                        auto pluginInstance = dynamic_cast<DerivativeTransformationPlugin*>(plugins().requestPlugin(getKind()));
                        if (pluginInstance) {
                            pluginInstance->setInputDataset(dataset);
                            pluginInstance->setKernel(kernel);
                            pluginInstance->setSavitzkyGolayParameters(sgWindowSize, sgPolynomialOrder);
                            pluginInstance->setGaussianSigma(sigma);
                            pluginInstance->transform();
                        }
                    }
                });

            pluginTriggerActions << pluginTriggerAction;
        };

        addTriggerAction(Kernel::Forward);
        addTriggerAction(Kernel::Central);
        addTriggerAction(Kernel::Central5);
        addTriggerAction(Kernel::SavitzkyGolay);
        addTriggerAction(Kernel::Gaussian);
    }

    return pluginTriggerActions;
}
