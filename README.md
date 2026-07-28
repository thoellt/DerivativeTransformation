# Derivative Transformation

A [ManiVault Studio](https://github.com/ManiVaultStudio) **Transformation** plugin that computes the first derivative of spectral response functions stored in a Points dataset.

## What it does

Each point in the input dataset is treated as one function whose dimensions are ordered samples with unit spacing (e.g. a spectral response curve). The plugin convolves every function with a precomputed per-sample weight table to produce its first derivative. Output dimensions match the input, and dimension names are prefixed with `d/dλ`.

It is triggered from the dataset right-click menu:

> **Transform → Derivative Transformation**

A single entry; everything it does is a setting.

### Settings

**Kernel** picks the derivative kernel (see below); it defaults to Savitzky–Golay.

**Output** decides where the result goes, and applies to every kernel:

| Mode | Description |
|---|---|
| **Add derived dataset** (default) | The derivative is written to a **new derived dataset** named after the input and the kernel; the input is left untouched. |
| **Replace existing dataset** | The derivative overwrites the values of the input dataset. Applied to a subset, only that subset's points are written and the dimension names are left alone, since they are shared with the parent. |

### Kernels

| Kernel | Description |
|---|---|
| **Forward differences** | First-order accurate |
| **Central differences** | Second-order accurate |
| **5-point central differences** | Fourth-order accurate |
| **Savitzky–Golay** (default) | Least-squares polynomial derivative; adds **Window size** and **Polynomial order** |
| **Derivative-of-Gaussian** | Gaussian-smoothed derivative; adds **Sigma** |

The parameters of all kernels sit in the same panel; the ones that do not belong to the
selected kernel are greyed out rather than removed, so what a kernel offers stays visible.

The settings are exposed through a configuration action on the trigger, so a host that offers one — the gear button next to a plugin trigger picker, for instance — can present them in place. The right-click entries have nowhere to show it and are marked with `...`; they ask for the same settings in a modal dialog when triggered. Both routes read and write the same values, and cancelling the dialog restores what was set before it opened.

The settings only offer combinations the kernel is defined for: the window size stays odd (stepping it moves it by two, in the direction it was stepped), and the polynomial order's upper bound follows the window size — capped at 10, above which fitting is numerically poor regardless of how much room the window has.

Boundary samples fall back to one-sided/asymmetric kernels that stay exact for linear signals, so the output keeps the input's dimensionality.

## Requirements

- ManiVault Studio core **1.5**, installed and available via CMake (`ManiVault_DIR`).
- Qt **6** (built and tested against 6.9.3) with the Widgets, Network, and Concurrent
  components.
- A C++20 compiler and CMake ≥ 3.22.

## Build

Configure (here uses Xcode generator), build, then launch ManiVault Studio (which loads the freshly
built plugin from its install directory). Adjust the paths to match your environment.

```bash
# Configure (out-of-source build directory next to the repo)
cmake -S . -B build -G Xcode \
  -DCMAKE_PREFIX_PATH=Qt6_DIR \
  -DManiVault_DIR=ManiVault_DIR

# Build
cmake --build build --config Debug --parallel
```

A post-build step installs the plugin straight into the ManiVault install directory, so no
manual copy is required. On Windows/Linux, drop the `-G Xcode` generator and use the
platform-appropriate way.

## Project layout

| Path | Purpose |
| --- | --- |
| `CMakeLists.txt` | Single shared-library target; links `ManiVault::Core` + `ManiVault::PointData`; auto-installs on build. |
| `PluginInfo.json` | Plugin metadata (name, version, type, dependencies) |
| `src/DerivativeTransformationPlugin.{h,cpp}` | Plugin class and factory. |

## License

LGPL v3.0.