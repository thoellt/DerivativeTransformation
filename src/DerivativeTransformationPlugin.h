#pragma once

#include <TransformationPlugin.h>

#include "DerivativeKernels.h"

#include <actions/DecimalAction.h>
#include <actions/GroupAction.h>
#include <actions/IntegralAction.h>
#include <actions/OptionAction.h>

#include <QMap>
#include <QString>

using namespace mv::plugin;

/**
 * DerivativeTransformation transformation plugin
 *
 * Computes the first derivative of spectral response functions stored in a
 * Points dataset: every data point is one function, the dimensions (in order)
 * are its discrete samples. Each spectrum is convolved with a selectable 1D
 * derivative kernel (unit sample spacing); the result has the same number of
 * samples (one-sided/asymmetric kernels are used near the boundaries) and is
 * written either over the input dataset or into a derived one.
 */
class DerivativeTransformationPlugin : public TransformationPlugin
{
    Q_OBJECT

public:
    /** Available derivative kernels (math lives in DerivativeKernels.{h,cpp}) */
    using Kernel = derivative::Kernel;

    /** Where the computed derivative is written */
    enum class Output {
        InPlace,        /** Overwrite the values of the input dataset */
        Derived         /** Write the values to a dataset derived from the input */
    };

    static const QMap<Kernel, QString> kernels;
    static const QMap<Output, QString> outputs;

public:
    DerivativeTransformationPlugin(const PluginFactory* factory);
    ~DerivativeTransformationPlugin() override = default;

    void init() override {};

    /** Called by the core after the input dataset was set */
    void transform() override;

    Kernel getKernel() const { return _kernel; }
    void setKernel(const Kernel& kernel) { _kernel = kernel; }
    static QString getKernelName(const Kernel& kernel) { return kernels[kernel]; }

    Output getOutput() const { return _output; }
    void setOutput(const Output& output) { _output = output; }
    static QString getOutputName(const Output& output) { return outputs[output]; }

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
    Output  _output;                /** Where the derivative is written */
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
    /** Available derivative kernels, spelled as the plugin spells them */
    using Kernel = DerivativeTransformationPlugin::Kernel;

    DerivativeTransformationPluginFactory();

    DerivativeTransformationPlugin* produce() override;

    mv::DataTypes supportedDataTypes() const override;

    /** One right-click "Transform" entry per output mode and kernel, bound to \p datasets */
    mv::gui::PluginTriggerActions getPluginTriggerActions(const mv::Datasets& datasets) const override;

    /**
     * Get plugin trigger actions for \p dataTypes
     *
     * The same entries, but for a caller that has no dataset yet and assigns one to the
     * action it picked before triggering it.
     *
     * @param dataTypes Vector of input data types
     * @return Vector of plugin trigger actions
     */
    mv::gui::PluginTriggerActions getPluginTriggerActions(const mv::DataTypes& dataTypes) const override;

    /**
     * Get the action that configures the transformation
     *
     * Handed to the trigger action so that a caller offering a configuration UI, such as the
     * gear button of a plugin trigger picker, can present it before triggering. Carries the
     * kernel, the output mode and the parameters of the kernels that have any.
     *
     * @return Pointer to the configuration action, never null
     */
    mv::gui::WidgetAction* getConfigurationAction();

    /**
     * Get the derivative kernel to convolve with
     * @return Derivative kernel
     */
    Kernel getKernel() const;

    /**
     * Get where the derivative should be written
     * @return Output mode
     */
    DerivativeTransformationPlugin::Output getOutput() const;

    /**
     * Get the Savitzky-Golay parameters
     *
     * Always a usable combination: the actions themselves are kept to one, so there is
     * nothing left to correct here.
     *
     * @param windowSize Sliding window size in samples, odd
     * @param polynomialOrder Order of the fitted polynomial, smaller than the window size
     */
    void getSavitzkyGolayParameters(int& windowSize, int& polynomialOrder) const;

    /**
     * Get the derivative-of-Gaussian standard deviation in samples
     * @return Standard deviation
     */
    float getGaussianSigma() const;

private:

    /**
     * Build the trigger action
     *
     * @param datasets Datasets to bind, empty to have the action resolve its own when triggered
     * @param configurable Whether the caller shows the configuration action itself, in which
     *                     case the settings are not asked for again when triggering
     * @return Vector of plugin trigger actions
     */
    mv::gui::PluginTriggerActions createTriggerActions(const mv::Datasets& datasets, bool configurable) const;

    /**
     * Ask for the settings in a modal dialog
     *
     * For callers that trigger the action straight away and so never get to show the
     * configuration action themselves, such as the dataset right-click menu. Shows the very
     * same action, and puts its values back when the dialog is cancelled.
     *
     * @return Boolean determining whether the transformation should go ahead
     */
    bool editSettings();

    /**
     * Grey out the parameters that do not belong to the selected kernel
     *
     * The configuration action is handed out once and cannot be swapped afterwards, so every
     * parameter lives in the one group and the ones that do not apply are disabled rather
     * than taken out of it.
     */
    void constrainParametersToKernel();

    /**
     * Keep the Savitzky-Golay actions to combinations the kernel is defined for
     *
     * Rounds the window size to odd in whichever direction it was moved, and caps the
     * polynomial order below it. Called whenever the window size changes, and once at
     * construction to bring the initial values under the same rule.
     */
    void constrainSavitzkyGolayParameters();


    mv::gui::OptionAction   _kernelAction;              /** Derivative kernel to convolve with */
    mv::gui::OptionAction   _outputAction;              /** Where the derivative is written */
    mv::gui::IntegralAction _sgWindowSizeAction;        /** Savitzky-Golay window size (samples, odd) */
    mv::gui::IntegralAction _sgPolynomialOrderAction;   /** Savitzky-Golay polynomial order */
    mv::gui::DecimalAction  _sigmaAction;               /** Derivative-of-Gaussian standard deviation (samples) */
    mv::gui::GroupAction    _groupAction;               /** Configuration UI holding all of the above */
    int                     _sgWindowSizeLast;          /** Last window size seen, to tell which way it is being moved */
};
