#include "DerivativeKernels.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace derivative {

namespace {

WeightRow makeRow(int start, const std::vector<double>& weights)
{
    WeightRow row;
    row.start = start;
    row.weights.assign(weights.begin(), weights.end());
    return row;
}

} // namespace

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

    // Augmented system [M | e1], M_{ab} = moments[a + b]; solving (A^T A) y = e1
    // with A_{ji} = (j - t)^i yields the derivative weights as w = A y
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

WeightTable buildWeightTable(Kernel kernel, int numSamples, int sgWindowSize, int sgPolynomialOrder, double sigma)
{
    const auto forwardAt = [numSamples](int d) -> WeightRow {
        return { std::min(d, numSamples - 2), { -1.0f, 1.0f } };
    };

    const auto centralAt = [&](int d) -> WeightRow {
        if (d == 0 || d == numSamples - 1)
            return forwardAt(d);
        return { d - 1, { -0.5f, 0.0f, 0.5f } };
    };

    WeightTable table;
    auto& rows = table.rows;
    rows.resize(numSamples);

    switch (kernel)
    {
        case Kernel::Forward:
            for (int d = 0; d < numSamples; ++d)
                rows[d] = forwardAt(d);

            table.interiorBegin = 0;
            table.interiorEnd = numSamples - 1;
            table.interiorOffset = 0;
            break;

        case Kernel::Central:
            for (int d = 0; d < numSamples; ++d)
                rows[d] = centralAt(d);

            table.interiorBegin = 1;
            table.interiorEnd = numSamples - 1;
            table.interiorOffset = 1;
            break;

        case Kernel::Central5:
            for (int d = 0; d < numSamples; ++d) {
                if (d >= 2 && d <= numSamples - 3)
                    rows[d] = { d - 2, { 1.0f / 12.0f, -8.0f / 12.0f, 0.0f, 8.0f / 12.0f, -1.0f / 12.0f } };
                else
                    rows[d] = centralAt(d);
            }

            if (numSamples >= 5) {
                table.interiorBegin = 2;
                table.interiorEnd = numSamples - 2;
                table.interiorOffset = 2;
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

                table.interiorBegin = 1;
                table.interiorEnd = numSamples - 1;
                table.interiorOffset = 1;
                break;
            }

            const int polynomialOrder = std::clamp(sgPolynomialOrder, 1, windowSize - 1);
            const int halfWindow = (windowSize - 1) / 2;

            std::map<int, WeightRow> rowsByPosition;

            for (int d = 0; d < numSamples; ++d) {
                const int start = std::clamp(d - halfWindow, 0, numSamples - windowSize);
                const int evaluationIndex = d - start;

                auto cached = rowsByPosition.find(evaluationIndex);
                if (cached == rowsByPosition.end())
                    cached = rowsByPosition.emplace(evaluationIndex, makeRow(0, savitzkyGolayWeights(windowSize, polynomialOrder, evaluationIndex))).first;

                rows[d] = cached->second;
                rows[d].start = start;
            }

            table.interiorBegin = halfWindow;
            table.interiorEnd = numSamples - halfWindow;
            table.interiorOffset = halfWindow;
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

                rows[d] = makeRow(start, weights);
            }

            if (numSamples >= 2 * radius + 1) {
                table.interiorBegin = radius;
                table.interiorEnd = numSamples - radius;
                table.interiorOffset = radius;
            }
            break;
        }
    }

    if (table.interiorBegin < table.interiorEnd)
        table.interiorKernel = rows[table.interiorBegin].weights;
    else
        table.interiorBegin = table.interiorEnd = 0;

    return table;
}

void convolveRange(const float* input, float* output, std::uint32_t firstPoint, std::uint32_t lastPoint, int numDimensions, const WeightTable& table)
{
    const auto kernelSize = static_cast<int>(table.interiorKernel.size());
    const float* kernel = table.interiorKernel.data();

    for (std::uint32_t p = firstPoint; p < lastPoint; ++p) {
        const float* spectrum = input + static_cast<std::size_t>(p) * numDimensions;
        float* out = output + static_cast<std::size_t>(p) * numDimensions;

        const auto edgeSample = [&](int d) {
            const auto& row = table.rows[d];

            float sum = 0.0f;
            for (std::size_t j = 0; j < row.weights.size(); ++j)
                sum += row.weights[j] * spectrum[row.start + j];

            out[d] = sum;
        };

        for (int d = 0; d < table.interiorBegin; ++d)
            edgeSample(d);

        // uniform interior: contiguous fixed kernel, auto-vectorizable
        for (int d = table.interiorBegin; d < table.interiorEnd; ++d) {
            const float* window = spectrum + d - table.interiorOffset;

            float sum = 0.0f;
            for (int j = 0; j < kernelSize; ++j)
                sum += kernel[j] * window[j];

            out[d] = sum;
        }

        for (int d = std::max(table.interiorEnd, table.interiorBegin); d < numDimensions; ++d)
            edgeSample(d);
    }
}

} // namespace derivative
