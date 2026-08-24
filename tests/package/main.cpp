// Consume an installed GSHTrans through find_package, exercising enough of it
// that a package which resolves but does not work would still fail here.
//
// It reaches through the umbrella header, uses all four dependencies
// indirectly -- the quadrature from GaussQuad, the FFT from FFTWpp, the
// concepts from NumericConcepts, the threading from OpenMP -- and checks a
// round trip, because a transform that returns the wrong numbers is a broken
// package just as surely as one that fails to link.

#include <GSHTrans/All>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

int main() {
  using namespace GSHTrans;

  static_assert(GSHTRANS_VERSION >= GSHTRANS_VERSION_NUMBER(1, 0, 0),
                "The version header did not come through the package");

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  constexpr auto lMax = std::ptrdiff_t{8};
  auto grid = Grid(lMax, 0);

  // A band-limited field, so that the round trip is an identity rather than a
  // projection: Y_1^0 up to normalisation.
  auto f = SpinField<0, Grid>(
      grid, [](auto theta, auto) { return Complex{std::cos(theta), 0}; });

  auto twice = Evaluate(Expand(f, lMax));

  auto drift = Real{0};
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      drift = std::max(drift, std::abs(twice[iTheta, iPhi] - f[iTheta, iPhi]));
    }
  }

  std::cout << "GSHTrans " << VersionString << ": round-trip drift " << drift
            << '\n';

  if (!(drift < 1e-12)) {
    std::cerr << "The round trip did not reproduce the field\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
