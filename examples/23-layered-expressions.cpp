// 23 -- Products over a stack, on the caller's own arrays
//
// What this shows. The pseudospectral step of a three-dimensional code:
// coefficients at every radius go to the grid, are multiplied there by a
// lateral weight, and come back as coefficients. Here the arrays belong to
// the caller -- they are the code's own workspaces, kept from step to step --
// and the library works on them in place through layered views, without
// copying them into its own types.
//
// Two dimensions stay the primitive. A layered expression such as
// `weight * strain` is a map from a radius to the two-dimensional expression
// `weight.Slice(i) * strain.Slice(i)`, so the upper-index and reality rules
// are example 03's and 04's: a real weight at N = 0 times a strain component
// at N = -2 is a complex field at N = -2, and the sum of the two would not
// compile. A depth-independent weight enters through AtEveryRadius, which
// stores nothing per radius.
//
// Expanding an expression evaluates it a block of radii at a time into
// scratch and transforms each block in one batched call, so the product is
// never stored whole -- and the result is the bits of expanding the stored
// product, which the example checks.
//
// Read first. 05 (views), 07 (batches and threads), 16 (layered fields).
//
// Introduced. LayeredSpinFieldView, LayeredSpinExpansionView and their const
// forms; layered expressions; AtEveryRadius; Assign and Materialise; the
// writing forms of Expand and Evaluate.
//
// The output. The sizes involved, and the largest difference between the
// expression route and the stored-product route, which is zero.

#include <GSHTrans/GSHTrans.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <span>
#include <vector>

int main() {
  using namespace GSHTrans;
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  constexpr auto lMax = std::ptrdiff_t{16};
  constexpr auto nR = std::ptrdiff_t{40};
  const auto grid = Grid::ForBand(lMax, 2, 1.5, FFTWpp::Estimate);
  auto radii = std::vector<Real>(nR);
  for (auto i = std::ptrdiff_t{0}; i < nR; ++i) radii[i] = 0.55 + 0.01 * i;
  const auto radial = RadialGrid<Real>(radii);
  const auto policy = Execution::Parallel();

  // The caller's arrays: coefficients of one strain component (N = -2) at
  // every radius, its values on the grid, a real lateral weight, and the
  // analysed product. Radius-major, as such codes lay them out.
  const auto coefficients =
      static_cast<std::size_t>(nR * grid.CoefficientSize(lMax, -2));
  const auto points = static_cast<std::size_t>(nR * grid.FieldSize());
  auto spectral = std::vector<Complex>(coefficients);
  auto values = std::vector<Complex>(points);
  auto weights = std::vector<Real>(points);
  auto analysed = std::vector<Complex>(coefficients);
  for (auto k = std::size_t{0}; k < coefficients; ++k) {
    spectral[k] = Complex{std::sin(0.01 * k), std::cos(0.013 * k)};
  }
  for (auto k = std::size_t{0}; k < points; ++k) {
    weights[k] = 1.0 + 0.2 * std::sin(0.001 * k);
  }

  // Views over them. Nothing is copied.
  const auto strainCoefficients = ConstLayeredSpinExpansionView<-2, Grid>(
      radial, grid, lMax, std::span<const Complex>(spectral));
  auto strain = LayeredSpinFieldView<-2, Grid>(radial, grid, std::span(values));
  const auto weight = ConstLayeredSpinFieldView<0, Grid, RealValued>(
      radial, grid, std::span<const Real>(weights));
  auto product = LayeredSpinExpansionView<-2, Grid>(radial, grid, lMax,
                                                    std::span(analysed));

  // The step: to the grid, weight, and back -- into the caller's arrays.
  Evaluate(strainCoefficients, strain, policy);
  Expand(weight * strain, product, policy);

  // The same through a stored product, for comparison.
  const auto stored = Materialise(weight * strain, policy);
  const auto expected = Expand(stored, lMax, policy);
  auto largest = Real{0};
  for (auto k = std::size_t{0}; k < coefficients; ++k) {
    largest = std::max(largest, std::abs(analysed[k] - expected.Data()[k]));
  }

  // A depth-independent weight, explicitly the same at every radius.
  auto flat = SpinField<0, Grid, RealValued>(grid);
  for (auto& w : flat.Data()) w = 1.25;
  const auto scaled =
      Expand(AtEveryRadius(radial, flat) * strain, lMax, policy);
  auto largestScaled = Real{0};
  for (auto k = std::size_t{0}; k < coefficients; ++k) {
    largestScaled = std::max(largestScaled,
                             std::abs(scaled.Data()[k] - 1.25 * spectral[k]));
  }

  std::cout << "radii " << nR << ", grid " << grid.NumberOfCoLatitudes()
            << " x " << grid.NumberOfLongitudes() << ", degree " << lMax << "\n"
            << "expression route against stored product: largest difference "
            << largest << "\n"
            << "constant weight 1.25 against scaling the coefficients: "
            << "largest difference " << largestScaled << "\n";

  // The first is an identity, bit for bit; the second holds to round-off,
  // the grid being exact for a band-limited field times a constant.
  return (largest == 0 && largestScaled < 1e-12) ? 0 : 1;
}
