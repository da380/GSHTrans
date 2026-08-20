// 06 -- The transform
//
// Between the spatial and spectral representations. The one thing worth
// understanding before using it: transforming an arbitrary field is a
// *projection*, not a round trip.

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <iostream>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  constexpr auto lMax = std::ptrdiff_t{16};
  auto grid = Grid(lMax, 2);

  auto f = SpinField<2, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::sin(theta) * std::cos(phi), std::cos(theta)};
  });

  // Expand takes any spin-weighted node -- a field, a view, or a lazy
  // expression -- so nothing needs materialising first.
  auto e = Expand(f, lMax);
  static_assert(decltype(e)::UpperIndex == 2);

  // The coefficients start at degree |N|. Below that no harmonic of that
  // upper index exists, so there is nothing to hold: d^l_{mN} vanishes
  // identically for l < |N|.
  std::cout << "degrees " << e.MinDegree() << " .. " << e.MaxDegree()
            << ",  " << e.Size() << " coefficients\n"
            << "f^2_{2,0} = " << (e[2, 0]) << "\n\n";

  // Now the projection. The field above has content at degrees 0 and 1 that
  // a spin-2 field cannot carry, and the transform discards it -- so
  // evaluating the coefficients does not give back what we started with.
  auto projected = Evaluate(e);
  auto drift = Real{0};
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      drift = std::max(drift, std::abs(projected[iTheta, iPhi] - f[iTheta, iPhi]));
    }
  }
  std::cout << "arbitrary field, once through: " << drift
            << "   (a projection, not an error)\n";

  // Once the field is band-limited, the round trip is an identity.
  auto again = Expand(projected, lMax);
  auto twice = Evaluate(again);
  drift = 0;
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      drift = std::max(
          drift, std::abs(twice[iTheta, iPhi] - projected[iTheta, iPhi]));
    }
  }
  std::cout << "band-limited field, again:     " << drift << "\n\n";

  // A real scalar field uses the reduced m >= 0 storage, since its negative
  // orders follow from f_{l,-m} = (-1)^m conj(f_{lm}). Half the coefficients
  // for the same information.
  auto g = SpinField<0, Grid, RealValued>(
      grid, [](auto theta, auto phi) { return std::cos(theta) * std::sin(phi); });
  auto real = Expand(g, lMax);
  auto complexOne = Expand(Materialise(g * Complex{1.0, 0.0}), lMax);
  std::cout << "real scalar    " << real.Size() << " coefficients\n"
            << "complex scalar " << complexOne.Size() << "\n";

  FFTWpp::CleanUp();
}
