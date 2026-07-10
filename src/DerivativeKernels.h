#pragma once

#include <cstdint>
#include <vector>

/**
 * Kernel-based first-derivative machinery for uniformly sampled functions
 * (unit sample spacing). Pure math without Qt/ManiVault dependencies, so it
 * can be compiled and tested standalone.
 */
namespace derivative {

/** Available derivative kernels */
enum class Kernel {
    Forward,        /** Forward differences, first order accurate */
    Central,        /** Central differences, second order accurate */
    Central5,       /** Five-point central differences, fourth order accurate */
    SavitzkyGolay,  /** Savitzky-Golay least-squares derivative (window size, polynomial order) */
    Gaussian        /** Derivative-of-Gaussian kernel (sigma) */
};

/** Weights of one output sample: out[d] = sum_j weights[j] * f[start + j] */
struct WeightRow {
    int                 start = 0;
    std::vector<float>  weights;
};

/**
 * Convolution plan for one spectrum length. Boundary samples use their
 * per-sample rows; the uniform interior [interiorBegin, interiorEnd) uses the
 * single contiguous interiorKernel (input window starts at d - interiorOffset),
 * which keeps the hot loop auto-vectorizable.
 */
struct WeightTable {
    std::vector<WeightRow>  rows;               /** Per-sample rows (consulted outside the interior) */
    std::vector<float>      interiorKernel;     /** Contiguous kernel of the uniform interior */
    int                     interiorBegin = 0;  /** First interior sample */
    int                     interiorEnd = 0;    /** One past the last interior sample (empty: begin == end) */
    int                     interiorOffset = 0; /** Interior window for sample d starts at d - interiorOffset */
};

/**
 * Savitzky-Golay first-derivative weights: least-squares fit of a polynomial of
 * order @p polynomialOrder to @p windowSize samples, derivative evaluated at
 * window position @p evaluationIndex (asymmetric positions handle boundaries).
 */
std::vector<double> savitzkyGolayWeights(int windowSize, int polynomialOrder, int evaluationIndex);

/**
 * Per-output-sample weight rows plus interior descriptor for a spectrum of
 * @p numSamples samples. Boundary samples fall back to one-sided/asymmetric
 * variants that stay exact for linear signals.
 */
WeightTable buildWeightTable(Kernel kernel, int numSamples, int sgWindowSize, int sgPolynomialOrder, double sigma);

/**
 * Differentiate the points [firstPoint, lastPoint): @p input and @p output are
 * point-major buffers of numPoints x @p numDimensions floats. Thread-safe for
 * disjoint point ranges on shared buffers.
 */
void convolveRange(const float* input, float* output, std::uint32_t firstPoint, std::uint32_t lastPoint, int numDimensions, const WeightTable& table);

}
