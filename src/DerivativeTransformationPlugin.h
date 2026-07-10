#pragma once

#include <TransformationPlugin.h>

#include "DerivativeKernels.h"

#include <QMap>
#include <QString>

using namespace mv::plugin;

/**
 * DerivativeTransformation transformation plugin
 *
 * Computes the first derivative of spectral response functions stored in a
 * Points dataset: every data point is one function, the dimensions (in order)
 * are its discrete samples. Each spectrum is convolved with a selectable 1D
 * derivative kernel (unit sample spacing); the result is written to a derived
 * dataset with the same number of samples (one-sided/asymmetric kernels are
 * used near the boundaries).
 */
class DerivativeTransformationPlugin : public TransformationPlugin
{
    Q_OBJECT

public:
    /** Available derivative kernels (math lives in DerivativeKernels.{h,cpp}) */
    using Kernel = derivative::Kernel;

    static const QMap<Kernel, QString> kernels;

public:
    DerivativeTransformationPlugin(const PluginFactory* factory);
    ~DerivativeTransformationPlugin() override = default;

    void init() override {};

    /** Called by the core after the input dataset was set */
    void transform() override;

    Kernel getKernel() const { return _kernel; }
    void setKernel(const Kernel& kernel) { _kernel = kernel; }
    static QString getKernelName(const Kernel& kernel) { return kernels[kernel]; }

    /**
     * Set the Savitzky-Golay parameters
     * @param windowSize Sliding window size in samples (odd, >= 3)
     * @param polynomialOrder Order of the fitted polynomial (>= 1, < windowSize)
     */
    void setSavitzkyGolayParameters(int windowSize, int polynomialOrder) {
        _sgWindowSize = windowSize;
        _sgPolynomialOrder = polynomialOrder;
    }

    /**
     * Set the derivative-of-Gaussian parameter
     * @param sigma Standard deviation of the Gaussian in samples
     */
    void setGaussianSigma(float sigma) { _sigma = sigma; }

private:
    Kernel  _kernel;                /** Selected derivative kernel */
    int     _sgWindowSize;          /** Savitzky-Golay window size (samples, odd) */
    int     _sgPolynomialOrder;     /** Savitzky-Golay polynomial order */
    float   _sigma;                 /** Derivative-of-Gaussian standard deviation (samples) */
};

class DerivativeTransformationPluginFactory : public TransformationPluginFactory
{
    Q_INTERFACES(mv::plugin::TransformationPluginFactory mv::plugin::PluginFactory)
    Q_OBJECT
    Q_PLUGIN_METADATA(IID   "studio.manivault.DerivativeTransformationPlugin"
                      FILE  "PluginInfo.json")

public:
    DerivativeTransformationPluginFactory();

    DerivativeTransformationPlugin* produce() override;

    mv::DataTypes supportedDataTypes() const override;

    /** One right-click "Transform" entry per kernel */
    mv::gui::PluginTriggerActions getPluginTriggerActions(const mv::Datasets& datasets) const override;
};
