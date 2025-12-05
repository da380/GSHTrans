# 🌐 GSHTrans: Generalized Spherical Harmonics Transformation Library

## Overview

**GSHTrans** is a modern C++ library designed for working with fields defined on the sphere, specifically focusing on **Canonical Component Fields** and their expansion into **Generalized Spherical Harmonics (GSH)**.

The library employs generic programming, C++20 concepts, and the Curiously Recurring Template Pattern (CRTP) to provide a type-safe, performant, and flexible framework suitable for high-performance numerical simulations in fields like geophysics, computational fluid dynamics, and spherical analysis.

## Key Features

* **Canonical Component Fields ($\Psi_n$):** Provides a robust implementation of canonical component fields, allowing for vector and tensor fields to be analyzed and synthesized in the GSH basis. The upper index $n$ is tracked as a template parameter (`std::ptrdiff_t _N`).
* **Expression Templates (Zero-Overhead Operations):** Uses expression templates to define lazy-evaluated arithmetic operations (addition, subtraction, multiplication, and division) on fields. This minimizes temporary object creation and enables compiler optimizations for efficient computation.
* **Wigner $d$-Functions (Numerical Core):** Implements the Wigner $d$-functions, $d^l_{mn}(\theta)$, a critical component for spherical transformations. The class `Wigner` is template-parameterized by floating-point type, normalization, and indexing options for flexibility.
    * The values are computed using stable recurrence relations.
    * The computation is parallelized using OpenMP (`#pragma omp parallel for`).
* **Flexible Grid Support:** Supports various computational grids, with `GridBase` serving as the generic interface. The provided `GaussLegendreGrid` implements a spectral-space-compatible grid with highly efficient FFT-based forward and inverse transformations.
* **Modern C++ Design:** Built with C++23, leveraging features like `std::ranges`, `std::complex`, and C++ concepts (`RealFloatingPoint`, `ComplexFloatingPoint`, etc.) for strong type checking and cleaner code.
* **Efficient Data Handling:** Uses `FFTWpp::vector` for memory management, which is optimized for use with the FFTW library, facilitating highly efficient Fast Fourier Transforms (FFTs) for longitudinal processing.

---

## Expression Template Structure (Field Algebra)

The library is designed for natural, math-like syntax:

```cpp
// Create two fields
CanonicalComponentField<N, Grid, ComplexValued> u0(grid);
CanonicalComponentField<N, Grid, ComplexValued> u1(grid);

// Expression template for addition (u_sum is a temporary view, not a new vector)
auto u_sum = u0 + u1;

// Expression template for multiplication (updates the upper index N)
auto u_product = u0 * u1; // Resulting type has upper index N0 + N1
```

The arithmetic operators (`+`, `-`, `*`, `/`) are overloaded to return lightweight **view classes** (e.g., `CanonicalComponentFieldAdd`, `CanonicalComponentFieldMultiply`) instead of calculating the result immediately. The actual computation occurs only when the result is assigned to a concrete, non-view field object.

| Operation | Resulting Type (Upper Index) | Class |
| :--- | :--- | :--- |
| `u0 + u1` | $N_0$ (requires $N_0 = N_1$) | `CanonicalComponentFieldAdd` |
| `u0 - u1` | $N_0$ (requires $N_0 = N_1$) | `CanonicalComponentFieldSubtract` |
| `u0 * u1` | $N_0 + N_1$ | `CanonicalComponentFieldMultiply` |
| `u0 / u1` | $N_0$ (requires $N_1 = 0$) | `CanonicalComponentFieldDivide` |
| `-u` | $N$ | `CanonicalComponentFieldUnary` |
| `conj(u)` | $N$ | `CanonicalComponentFieldConj` |
| `real(u)` | $N$ (Value becomes `RealValued`) | `CanonicalComponentFieldReal` |
| `imag(u)` | $N$ (Value becomes `RealValued`) | `CanonicalComponentFieldImag` |

---

## Building the Library

The project uses CMake for configuration.

### Prerequisites

* A C++23 compliant compiler (e.g., GCC 13+ or Clang 16+).
* CMake (minimum version 3.20).
* FFTW (required by the `FFTWpp` dependency).

### External Dependencies

All major external dependencies are handled automatically using CMake's `FetchContent`:

* `Eigen3` (version 3.4.1)
* `GaussQuad`
* `FFTWpp`

* `NumericConcepts`
* OpenMP (for parallel computation in Wigner $d$-function calculation)

### CMake Options

The following options can be set to control the build process:

* `MY_PROJECT_BUILD_EXAMPLES`: Builds example executables (default: `ON`).
* `MY_PROJECT_BUILD_TESTS`: Builds the test suite and enables testing (default: `ON`).

---

## License

This project is released under the **BSD 3-Clause License**.
Copyright (c) 2025, David Al-Attar.
