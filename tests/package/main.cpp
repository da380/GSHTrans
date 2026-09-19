// Consume an installed GSHTrans through find_package, exercising enough of it
// that a package which resolves but does not work would still fail here.
//
// It reaches through the umbrella header, uses all four dependencies
// indirectly -- the quadrature from GaussQuad, the FFT from FFTWpp, the
// concepts from NumericConcepts, the threading from OpenMP -- and checks a
// round trip, because a transform that returns the wrong numbers is a broken
// package just as surely as one that fails to link.
//
// The two optional dependencies are exercised when the installed package says
// it has them, and only then. Each is reached through a call that needs its
// symbols -- a GEMM from the BLAS, a spline from Interpolation -- because a
// dependency the package forgot to carry shows up at configure or link time
// only if something asks for it. The BLAS was forgotten once, and this file
// did not notice because it never asked.

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

  // A band-limited field, so that the round trip is an identity rather than a
  // projection: Y_1^0 up to normalisation.
  const auto shape = [](auto theta, auto) {
    return Complex{std::cos(theta), 0};
  };

  const auto roundTrip = [&](const Grid& grid, const char* name) {
    auto f = SpinField<0, Grid>(grid, shape);
    auto twice = Evaluate(Expand(f, lMax));

    auto drift = Real{0};
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        drift =
            std::max(drift, std::abs(twice[iTheta, iPhi] - f[iTheta, iPhi]));
      }
    }
    std::cout << "GSHTrans " << VersionString << ", " << name
              << ": round-trip drift " << drift << '\n';
    return drift < 1e-12;
  };

  auto grid = Grid(lMax, 0);
  if (!roundTrip(grid, "loop kernel")) {
    std::cerr << "The round trip did not reproduce the field\n";
    return EXIT_FAILURE;
  }

#ifdef GSHTRANS_HAVE_BLAS
  auto matrix = Grid(lMax, 0, FFTWpp::Estimate, Chunking::Automatic(),
                     WignerValues::Stored(), TransformKernel::Matrix());
  if (!roundTrip(matrix, "matrix kernel")) {
    std::cerr << "The matrix kernel did not reproduce the field\n";
    return EXIT_FAILURE;
  }
#else
  std::cout << "This package was built without a BLAS\n";
#endif

#ifdef GSHTRANS_HAVE_INTERPOLATION
  {
    // Fourth order on a grid this coarse is worth two or three digits, which
    // is enough to tell a spline from a link error.
    auto f = SpinField<0, Grid>(grid, shape);
    const auto bicubic = Interpolate(f, Scheme::Bicubic());
    const auto error = std::abs(bicubic(1.1, 2.2) - shape(1.1, 2.2));
    std::cout << "bicubic interpolation error " << error << '\n';
    if (!(error < 1e-2)) {
      std::cerr << "The bicubic interpolant is wrong\n";
      return EXIT_FAILURE;
    }
  }
#else
  std::cout << "This package was built without Interpolation\n";
#endif

  return EXIT_SUCCESS;
}
