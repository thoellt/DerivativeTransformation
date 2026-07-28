#include "DerivativeTransformationPlugin.h"

#include <PointData/PointData.h>

#include <actions/DecimalAction.h>
#include <actions/IntegralAction.h>
#include <actions/WidgetActionViewWidget.h>

#include <QDebug>
#include <QDialog>
#include <QDialogButtonBox>
#include <QEventLoop>
#include <QFutureWatcher>
#include <QGridLayout>
#include <QGroupBox>
#include <QThread>
#include <QTimer>
#include <QVBoxLayout>
#include <QtConcurrent/QtConcurrent>

#include <algorithm>
#include <atomic>
#include <type_traits>
#include <utility>
#include <vector>

Q_PLUGIN_METADATA(IID "studio.manivault.DerivativeTransformationPlugin")

using namespace mv;

using Kernel = DerivativeTransformationPlugin::Kernel;
using Output = DerivativeTransformationPlugin::Output;

const QMap<Kernel, QString> DerivativeTransformationPlugin::kernels = QMap<Kernel, QString>({
    { Kernel::Forward, "Forward differences" },
    { Kernel::Central, "Central differences" },
    { Kernel::Central5, "Central differences (5-point)" },
    { Kernel::SavitzkyGolay, "Savitzky-Golay" },
    { Kernel::Gaussian, "Gaussian derivative" }
});

const QMap<Output, QString> DerivativeTransformationPlugin::outputs = QMap<Output, QString>({
    { Output::InPlace, "Replace existing dataset" },
    { Output::Derived, "Add derived dataset" }
});

namespace {

/**
 * Highest Savitzky-Golay polynomial order offered, regardless of the window size. Fitting a
 * high order polynomial over a handful of samples is numerically poor long before the window
 * runs out of room for it.
 */
constexpr int maxPolynomialOrder = 10;

} // namespace

DerivativeTransformationPlugin::DerivativeTransformationPlugin(const PluginFactory* factory) :
    TransformationPlugin(factory),
    _kernel(Kernel::SavitzkyGolay),
    _output(Output::Derived),
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
    task.setProgressDescription(QString("1st derivative, %1 (%2)").arg(kernelDescription, getOutputName(_output)));

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

    // The samples the result holds are derivatives, so they are named for that. The names
    // of the input only carry over when there is one per sample, which a dataset that was
    // never given dimension names does not have.
    const auto derivativeDimensionNames = [&points, numDimensions]() -> std::vector<QString> {
        auto dimensionNames = points->getDimensionNames();

        if (static_cast<int>(dimensionNames.size()) != numDimensions) {
            dimensionNames.resize(numDimensions);

            for (int d = 0; d < numDimensions; ++d)
                dimensionNames[d] = QString::number(d);
        }

        for (auto& name : dimensionNames)
            name = QString("d/dλ %1").arg(name);

        return dimensionNames;
    };

    switch (_output)
    {
        case Output::InPlace:
        {
            points->setLocked(true);

            // Written back in the same point-major order it was gathered in; visitData walks
            // a subset in its own order, so gather and scatter line up only because both go
            // through it. Note that a dataset whose elements are not floating point rounds
            // the derivative to its own type here, which the derived output does not do.
            points->visitData([&derivativeValues, numDimensions](auto pointData) {
                std::size_t i = 0;

                for (auto point : pointData) {
                    for (int d = 0; d < numDimensions; ++d) {
                        using ValueType = std::remove_reference_t<decltype(point[d])>;

                        point[d] = static_cast<ValueType>(derivativeValues[i++]);
                    }
                }
            });

            points->setLocked(false);

            // Dimension names sit on the raw data that a subset shares with its parent, so
            // renaming them from a subset would relabel points that were left alone
            if (points->isFull())
                points->setDimensionNames(derivativeDimensionNames());

            events().notifyDatasetDataChanged(points);

            break;
        }

        case Output::Derived:
        {
            auto derived = mv::data().createDerivedDataset<Points>(
                points->getGuiName() + QString(" (1st derivative, %1)").arg(kernelDescription), points);

            derived->setData(std::move(derivativeValues), numDimensions);
            derived->setDimensionNames(derivativeDimensionNames());

            events().notifyDatasetDataChanged(derived);

            break;
        }
    }

    task.setFinished();
}

// -----------------------------------------------------------------------------
// Factory
// -----------------------------------------------------------------------------

DerivativeTransformationPluginFactory::Settings::Settings(QObject* parent, const QString& title, const Output& defaultOutput) :
    mv::gui::GroupAction(parent, title),
    _kernelAction(this, "Kernel", DerivativeTransformationPlugin::kernels.values(), DerivativeTransformationPlugin::getKernelName(Kernel::SavitzkyGolay)),
    _outputAction(this, "Output", { DerivativeTransformationPlugin::getOutputName(Output::InPlace), DerivativeTransformationPlugin::getOutputName(Output::Derived) }, DerivativeTransformationPlugin::getOutputName(defaultOutput)),
    _sgWindowSizeAction(this, "Window size", 3, 101, 7),
    _sgPolynomialOrderAction(this, "Polynomial order", 1, 10, 2),
    _sigmaAction(this, "Sigma", 0.1f, 25.0f, 1.0f, 2),
    _sgWindowSizeLast(_sgWindowSizeAction.getValue())
{
    _kernelAction.setToolTip("Derivative kernel every spectrum is convolved with");
    _outputAction.setToolTip("Whether the derivative overwrites the input dataset or is written to a derived one");
    _sgWindowSizeAction.setToolTip("Sliding window size in samples (odd; stepping moves it by two)");
    _sgPolynomialOrderAction.setToolTip("Order of the fitted polynomial (kept below the window size)");
    _sigmaAction.setToolTip("Standard deviation of the Gaussian in samples");

    // The kernel is only defined for some combinations of the two, so rather than correct
    // them behind the user's back when the transformation runs, the actions are held to
    // those combinations as they are edited — whichever UI is doing the editing
    connect(&_sgWindowSizeAction, &mv::gui::IntegralAction::valueChanged, this, [this]() -> void {
        constrainSavitzkyGolayParameters();
    });

    connect(&_kernelAction, &mv::gui::OptionAction::currentTextChanged, this, [this]() -> void {
        constrainParametersToKernel();
    });

    constrainSavitzkyGolayParameters();
    constrainParametersToKernel();

    setToolTip("Derivative settings");
    setLabelSizingType(mv::gui::GroupAction::LabelSizingType::Auto);

    addAction(&_kernelAction);
    addAction(&_outputAction);
    addAction(&_sgWindowSizeAction);
    addAction(&_sgPolynomialOrderAction);
    addAction(&_sigmaAction);
}

Kernel DerivativeTransformationPluginFactory::Settings::getKernel() const
{
    return DerivativeTransformationPlugin::kernels.key(_kernelAction.getCurrentText(), Kernel::SavitzkyGolay);
}

DerivativeTransformationPlugin::Output DerivativeTransformationPluginFactory::Settings::getOutput() const
{
    return DerivativeTransformationPlugin::outputs.key(_outputAction.getCurrentText(), Output::Derived);
}

void DerivativeTransformationPluginFactory::Settings::getSavitzkyGolayParameters(int& windowSize, int& polynomialOrder) const
{
    windowSize      = _sgWindowSizeAction.getValue();
    polynomialOrder = _sgPolynomialOrderAction.getValue();
}

float DerivativeTransformationPluginFactory::Settings::getGaussianSigma() const
{
    return _sigmaAction.getValue();
}

void DerivativeTransformationPluginFactory::Settings::constrainParametersToKernel()
{
    const auto kernel = getKernel();

    _sgWindowSizeAction.setForceDisabled(kernel != Kernel::SavitzkyGolay);
    _sgPolynomialOrderAction.setForceDisabled(kernel != Kernel::SavitzkyGolay);
    _sigmaAction.setForceDisabled(kernel != Kernel::Gaussian);
}

void DerivativeTransformationPluginFactory::Settings::constrainSavitzkyGolayParameters()
{
    const auto windowSize = _sgWindowSizeAction.getValue();

    // An even window has no centre sample to differentiate about. Which odd neighbour it
    // becomes follows the direction it was moved in, so that stepping the spin box up from 7
    // reaches 9 rather than being pulled back to 7 and looking broken.
    if (windowSize % 2 == 0) {
        _sgWindowSizeAction.setValue(windowSize > _sgWindowSizeLast ? windowSize + 1 : windowSize - 1);

        // Setting it re-enters here with an odd value, which does the rest of the work
        return;
    }

    _sgWindowSizeLast = windowSize;

    // A polynomial of the window's own order fits it exactly and so has no smoothing left in
    // it; the action's maximum keeps the UI from offering that in the first place
    const auto highestOrder = std::min(maxPolynomialOrder, windowSize - 1);

    _sgPolynomialOrderAction.setMaximum(highestOrder);

    // Lowering the maximum leaves the value where it was, so it has to be brought down by hand
    _sgPolynomialOrderAction.setValue(std::min(_sgPolynomialOrderAction.getValue(), highestOrder));
}

QWidget* DerivativeTransformationPluginFactory::Settings::getWidget(QWidget* parent, const std::int32_t& widgetFlags)
{
    auto settingsWidget = createSettingsWidget(parent);

    if (!(widgetFlags & mv::gui::WidgetActionViewWidget::PopupLayout))
        return settingsWidget;

    // In a popup the layout this override replaces would have boxed the whole group under its
    // own title, so that box has to be put back by hand here
    auto popupWidget = new QWidget(parent);
    auto popupLayout = new QVBoxLayout(popupWidget);
    auto groupBox    = new QGroupBox(text(), popupWidget);
    auto groupLayout = new QVBoxLayout(groupBox);

    popupLayout->setContentsMargins(4, 4, 4, 4);
    groupLayout->addWidget(settingsWidget);
    popupLayout->addWidget(groupBox);

    return popupWidget;
}

QWidget* DerivativeTransformationPluginFactory::Settings::createSettingsWidget(QWidget* parent)
{
    auto widget = new QWidget(parent);
    auto layout = new QVBoxLayout(widget);

    layout->setContentsMargins(0, 0, 0, 0);

    // What is computed is a separate decision from where the result goes, so the two are boxed
    // apart rather than laid out as one run of rows
    const auto addGroupBox = [widget, layout](const QString& title, const QVector<mv::gui::WidgetAction*>& actions, bool showLabels) -> void {
        auto groupBox   = new QGroupBox(title, widget);
        auto gridLayout = new QGridLayout(groupBox);

        for (auto action : actions) {
            const auto row = gridLayout->rowCount();

            if (showLabels)
                gridLayout->addWidget(action->createLabelWidget(groupBox), row, 0);

            gridLayout->addWidget(action->createWidget(groupBox), row, showLabels ? 1 : 0, 1, showLabels ? 1 : 2);
        }

        gridLayout->setColumnStretch(1, 1);

        layout->addWidget(groupBox);
    };

    addGroupBox("Parameters", { &_kernelAction, &_sgWindowSizeAction, &_sgPolynomialOrderAction, &_sigmaAction }, true);

    // The box already says what the drop-down is, so it carries no label of its own
    addGroupBox("Output", { &_outputAction }, false);

    return widget;
}

bool DerivativeTransformationPluginFactory::Settings::edit()
{
    // The very actions that are shown here are the ones kept afterwards, so cancelling has to
    // put back what they held before the dialog was opened
    const auto kernel           = _kernelAction.getCurrentText();
    const auto output           = _outputAction.getCurrentText();
    const auto windowSize       = _sgWindowSizeAction.getValue();
    const auto polynomialOrder  = _sgPolynomialOrderAction.getValue();
    const auto sigma            = _sigmaAction.getValue();

    QDialog dialog;

    dialog.setWindowTitle("Derivative Transformation");

    auto layout = new QVBoxLayout(&dialog);

    layout->addWidget(createSettingsWidget(&dialog));

    auto buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);

    QObject::connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);

    layout->addWidget(buttons);

    if (dialog.exec() == QDialog::Accepted)
        return true;

    _kernelAction.setCurrentText(kernel);
    _outputAction.setCurrentText(output);
    _sgWindowSizeAction.setValue(windowSize);
    _sgPolynomialOrderAction.setValue(polynomialOrder);
    _sigmaAction.setValue(sigma);

    return false;
}

DerivativeTransformationPluginFactory::DerivativeTransformationPluginFactory() :
    _hierarchySettings(this, "Derivative", Output::Derived),
    _hostSettings(this, "Derivative", Output::InPlace)
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

DerivativeTransformationPluginFactory::Settings& DerivativeTransformationPluginFactory::getSettings(bool configurable)
{
    return configurable ? _hostSettings : _hierarchySettings;
}

mv::DataTypes DerivativeTransformationPluginFactory::supportedDataTypes() const
{
    return { PointType };
}

mv::gui::PluginTriggerActions DerivativeTransformationPluginFactory::createTriggerActions(const mv::Datasets& datasets, bool configurable) const
{
    mv::gui::PluginTriggerActions pluginTriggerActions;

    auto factory = const_cast<DerivativeTransformationPluginFactory*>(this);

    // Where the settings are not shown alongside the entry they are asked for in a dialog,
    // which the ellipsis promises
    const auto menuName = configurable ? QString("Derivative Transformation") : QString("Derivative Transformation...");

    auto& settings = factory->getSettings(configurable);

    auto pluginTriggerAction = new mv::gui::PluginTriggerAction(factory, this,
        menuName, "Compute the first derivative of each point's spectrum", icon(),
        [this, &settings, datasets, configurable](mv::gui::PluginTriggerAction& pluginTriggerAction) -> void {
            // Where the caller offers no way of showing the configuration action, the
            // settings have to be asked for here; it edits the same actions either way
            if (!configurable && !settings.edit())
                return;

            // An action built from data types alone is bound to its datasets by whoever
            // triggers it, which may be long after the action itself was created
            const auto targets = datasets.isEmpty() ? pluginTriggerAction.getDatasets() : datasets;

            int sgWindowSize        = 0;
            int sgPolynomialOrder   = 0;

            settings.getSavitzkyGolayParameters(sgWindowSize, sgPolynomialOrder);

            for (const auto& dataset : targets) {
                auto pluginInstance = dynamic_cast<DerivativeTransformationPlugin*>(plugins().requestPlugin(getKind()));
                if (pluginInstance) {
                    pluginInstance->setInputDataset(dataset);
                    pluginInstance->setKernel(settings.getKernel());
                    pluginInstance->setOutput(settings.getOutput());
                    pluginInstance->setSavitzkyGolayParameters(sgWindowSize, sgPolynomialOrder);
                    pluginInstance->setGaussianSigma(settings.getGaussianSigma());
                    pluginInstance->transform();
                }
            }
        });

    pluginTriggerAction->setConfigurationAction(&settings);

    pluginTriggerActions << pluginTriggerAction;

    return pluginTriggerActions;
}

mv::gui::PluginTriggerActions DerivativeTransformationPluginFactory::getPluginTriggerActions(const mv::Datasets& datasets) const
{
    if (datasets.count() < 1 || !PluginFactory::areAllDatasetsOfTheSameType(datasets, PointType))
        return {};

    // Reached from the dataset right-click menu, which triggers straight away and shows
    // nothing of the configuration action
    return createTriggerActions(datasets, false);
}

mv::gui::PluginTriggerActions DerivativeTransformationPluginFactory::getPluginTriggerActions(const mv::DataTypes& dataTypes) const
{
    if (dataTypes.isEmpty() || dataTypes.count(PointType) != dataTypes.count())
        return {};

    // Asked for without any dataset in hand, so by a caller picking a transformation up front
    // — one that can show the configuration action and assigns the datasets before triggering
    return createTriggerActions({}, true);
}
