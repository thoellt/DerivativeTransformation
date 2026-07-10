#include "DerivativeTransformationPlugin.h"

#include <PointData/PointData.h>

#include <actions/DecimalAction.h>
#include <actions/IntegralAction.h>

#include <QApplication>
#include <QDebug>
#include <QDialog>
#include <QDialogButtonBox>
#include <QGridLayout>

#include <algorithm>
#include <cmath>
#include <map>
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

/** Weights of one output sample: out[d] = sum_j weights[j] * spectrum[start + j] */
struct WeightRow {
    int                 start = 0;
    std::vector<double> weights;
};

/**
 * Savitzky-Golay first-derivative weights: least-squares fit of a polynomial of
 * order @p polynomialOrder to @p windowSize samples, derivative evaluated at
 * window position @p evaluationIndex (asymmetric positions handle boundaries).
 * Solves the normal equations (A^T A) y = e1 with A_{ji} = (j - t)^i, then w = A y.
 */
std::vector<double> savitzkyGolayWeights(int windowSize, int polynomialOrder, int evaluationIndex)
{
    const int n = polynomialOrder + 1;

    std::vector<double> moments(2 * polynomialOrder + 1, 0.0);

    for (int j = 0; j < windowSize; ++j) {
        const double u = j - evaluationIndex;
        double power = 1.0;
        for (int k = 0; k < static_cast<int>(moments.size()); ++k) {
            moments[k] += power;
            power *= u;
        }
    }

    // Augmented system [M | e1], M_{ab} = moments[a + b]
    std::vector<std::vector<double>> system(n, std::vector<double>(n + 1, 0.0));

    for (int a = 0; a < n; ++a) {
        for (int b = 0; b < n; ++b)
            system[a][b] = moments[a + b];
        system[a][n] = (a == 1) ? 1.0 : 0.0;
    }

    // Gaussian elimination with partial pivoting
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        for (int row = col + 1; row < n; ++row)
            if (std::abs(system[row][col]) > std::abs(system[pivot][col]))
                pivot = row;
        std::swap(system[col], system[pivot]);

        for (int row = col + 1; row < n; ++row) {
            const double factor = system[row][col] / system[col][col];
            for (int k = col; k <= n; ++k)
                system[row][k] -= factor * system[col][k];
        }
    }

    std::vector<double> y(n, 0.0);
    for (int row = n - 1; row >= 0; --row) {
        double sum = system[row][n];
        for (int k = row + 1; k < n; ++k)
            sum -= system[row][k] * y[k];
        y[row] = sum / system[row][row];
    }

    std::vector<double> weights(windowSize, 0.0);
    for (int j = 0; j < windowSize; ++j) {
        const double u = j - evaluationIndex;
        double power = 1.0;
        for (int i = 0; i < n; ++i) {
            weights[j] += y[i] * power;
            power *= u;
        }
    }

    return weights;
}

/**
 * Per-output-sample weight rows for the whole spectrum (unit sample spacing).
 * Interior samples get the nominal kernel; boundary samples fall back to
 * one-sided/asymmetric variants that stay exact for linear signals.
 */
std::vector<WeightRow> buildWeightTable(Kernel kernel, int numSamples, int sgWindowSize, int sgPolynomialOrder, double sigma)
{
    const auto forwardAt = [numSamples](int d) -> WeightRow {
        return { std::min(d, numSamples - 2), { -1.0, 1.0 } };
    };

    const auto centralAt = [&](int d) -> WeightRow {
        if (d == 0 || d == numSamples - 1)
            return forwardAt(d);
        return { d - 1, { -0.5, 0.0, 0.5 } };
    };

    std::vector<WeightRow> rows(numSamples);

    switch (kernel)
    {
        case Kernel::Forward:
            for (int d = 0; d < numSamples; ++d)
                rows[d] = forwardAt(d);
            break;

        case Kernel::Central:
            for (int d = 0; d < numSamples; ++d)
                rows[d] = centralAt(d);
            break;

        case Kernel::Central5:
            for (int d = 0; d < numSamples; ++d) {
                if (d >= 2 && d <= numSamples - 3)
                    rows[d] = { d - 2, { 1.0 / 12.0, -8.0 / 12.0, 0.0, 8.0 / 12.0, -1.0 / 12.0 } };
                else
                    rows[d] = centralAt(d);
            }
            break;

        case Kernel::SavitzkyGolay:
        {
            int windowSize = std::min(sgWindowSize, numSamples);
            if (windowSize % 2 == 0)
                windowSize -= 1;

            if (windowSize < 3) {
                for (int d = 0; d < numSamples; ++d)
                    rows[d] = centralAt(d);
                break;
            }

            const int polynomialOrder = std::clamp(sgPolynomialOrder, 1, windowSize - 1);
            const int halfWindow = (windowSize - 1) / 2;

            std::map<int, std::vector<double>> weightsByPosition;

            for (int d = 0; d < numSamples; ++d) {
                const int start = std::clamp(d - halfWindow, 0, numSamples - windowSize);
                const int evaluationIndex = d - start;

                auto cached = weightsByPosition.find(evaluationIndex);
                if (cached == weightsByPosition.end())
                    cached = weightsByPosition.emplace(evaluationIndex, savitzkyGolayWeights(windowSize, polynomialOrder, evaluationIndex)).first;

                rows[d] = { start, cached->second };
            }
            break;
        }

        case Kernel::Gaussian:
        {
            const int radius = std::max(1, static_cast<int>(std::ceil(3.0 * sigma)));

            for (int d = 0; d < numSamples; ++d) {
                const int start = std::max(0, d - radius);
                const int end = std::min(numSamples - 1, d + radius);
                const int count = end - start + 1;

                std::vector<double> weights(count);
                double mean = 0.0;

                for (int j = 0; j < count; ++j) {
                    const double k = start + j - d;
                    weights[j] = k * std::exp(-k * k / (2.0 * sigma * sigma));
                    mean += weights[j];
                }
                mean /= count;

                // Enforce exactness for constant (sum w = 0) and linear (sum w*k = 1)
                // signals; required where the truncated kernel loses its symmetry.
                double scale = 0.0;
                for (int j = 0; j < count; ++j) {
                    weights[j] -= mean;
                    scale += weights[j] * (start + j - d);
                }

                if (std::abs(scale) < 1e-12) {
                    rows[d] = centralAt(d);
                    continue;
                }

                for (auto& weight : weights)
                    weight /= scale;

                rows[d] = { start, std::move(weights) };
            }
            break;
        }
    }

    return rows;
}

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

    const auto weightTable = buildWeightTable(_kernel, numDimensions, _sgWindowSize, _sgPolynomialOrder, _sigma);

    auto derived = mv::data().createDerivedDataset<Points>(
        points->getGuiName() + QString(" (1st derivative, %1)").arg(kernelDescription), points);

    std::vector<float> derivative(static_cast<std::size_t>(numPoints) * numDimensions);

    // visitData is subset-aware; do NOT use constVisitFromBeginToEnd here —
    // it walks the raw source buffer, which is larger than
    // numPoints * numDimensions when the input is a subset (buffer overflow).
    // (Non-const overload used read-only: the const one does not compile in core 1.5.)
    points->visitData([&](auto pointData) {
        std::vector<double> spectrum(numDimensions);
        std::size_t out = 0;
        std::uint32_t pointsProcessed = 0;

        for (const auto point : pointData) {
            for (int d = 0; d < numDimensions; ++d)
                spectrum[d] = static_cast<double>(point[d]);

            for (int d = 0; d < numDimensions; ++d) {
                const auto& row = weightTable[d];

                double sum = 0.0;
                for (std::size_t j = 0; j < row.weights.size(); ++j)
                    sum += row.weights[j] * spectrum[row.start + j];

                derivative[out++] = static_cast<float>(sum);
            }

            if (++pointsProcessed % 1000 == 0) {
                task.setProgress(static_cast<float>(pointsProcessed) / numPoints);
                QApplication::processEvents();
            }
        }
    });

    derived->setData(std::move(derivative), numDimensions);

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
