// 21 -- Interpolating a field, as a callable of the two angles
//
// What this shows. A field is samples on a fixed grid. Interpolate turns one
// -- or an expansion, or a lazy expression -- into a function you can
// evaluate at any (theta, phi), in one of three schemes:
//
//   Spectral -- sum the expansion. Exact for a band-limited field, and the
//               reference the other two are measured against. O(lMax^2) a
//               point.
//   Bilinear -- Interpolation's rectilinear linear scheme over the samples.
//   Bicubic  -- the same, tensor-product cubic. Both are O(1) a point and
//               need the Interpolation dependency.
//
// The callable models ScalarFunctionS2, which is what a field constructor
// takes -- so remeshing onto another grid is one line, and is probably the
// commonest use of the whole feature.
//
// Two things about the sphere that a rectilinear scheme does not know, and
// which the library fixes by handing it a padded grid rather than the field's
// own:
//
//   -- the longitudes stop one step short of 2 pi, so the last cell would
//      have nothing to interpolate against. A wrap column is added, exactly:
//      phi = 2 pi is phi = 0.
//   -- neither pole is a grid point, so every local scheme would extrapolate
//      there. Two polar rows are added, and they are computed from the
//      expansion rather than guessed -- which is why building a local
//      interpolant costs a forward transform.
//
// At a pole only one order survives (m = N at theta = 0, m = -N at
// theta = pi; example 14), so the value there is c exp(+-i N phi): a *row*,
// not a constant. That is not a defect. A spin-weighted field at a coordinate
// pole is genuinely not single-valued, because the frame e_+- depends on the
// azimuth of approach.
//
// Read first. 06 (the transform), 12 (expansions), 14 (the d-functions and
// their behaviour at the poles). docs/gshtrans-reference.tex,
// "Interpolation", has the polar-row formulae.
//
// Introduced. Interpolate(field, scheme) and Interpolate(expansion);
// Scheme::Spectral(), Scheme::Bilinear(), Scheme::Bicubic(); constructing a
// field from an interpolant; Grid::ForBand for oversampling.
//
// The output. Spectral interpolation agreeing with the transform at the
// grid points; the phase rotation of a spin-2 field at the north pole; the
// refusal of a colatitude outside [0, pi] and the wrapping of longitude; a
// remesh onto a finer grid; and, with the Interpolation dependency, the
// errors of the local schemes at the band limit and on an oversampled grid.
//

#include <GSHTrans/GSHTrans.hpp>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  using Int = std::ptrdiff_t;

  constexpr auto pi = std::numbers::pi_v<Real>;
  constexpr auto N = Int{2};
  constexpr auto band = Int{12};

  // A band-limited field of upper index 2, built from its coefficients so
  // that we know exactly what it is: every (l, m) set, with an amplitude
  // falling like 1/(l+1)^2.
  auto grid = Grid(band, 2);
  auto expansion = SpinExpansion<N, Grid>(grid, band);
  for (auto l : expansion.Degrees()) {
    for (auto m : expansion.Orders(l)) {
      const auto scale = 1 / static_cast<Real>((l + 1) * (l + 1));
      expansion[l, m] = Complex(scale, scale / 2);
    }
  }
  const auto field = Evaluate(expansion);

  std::cout << std::scientific << std::setprecision(3);

  //------------------------------------------------------------------------//
  //  The reference, and that it agrees with the transform where both apply  //
  //------------------------------------------------------------------------//

  // The spectral interpolant expands the field once and sums the expansion
  // at each query. At the grid's own points it must reproduce the samples,
  // which is what the loop below checks.
  const auto exact = Interpolate(field, Scheme::Spectral());

  {
    auto worst = Real{0};
    auto iTheta = Int{0};
    for (auto theta : grid.CoLatitudes()) {
      auto iPhi = Int{0};
      for (auto phi : grid.Longitudes()) {
        worst =
            std::max(worst, std::abs(exact(theta, phi) - field[iTheta, iPhi]));
        iPhi++;
      }
      iTheta++;
    }
    std::cout << "Spectral against the transform, at the grid's own points: "
              << worst << "\n";
    if (!(worst < 1e-12)) return 1;
  }

  //------------------------------------------------------------------------//
  //                     The poles, and the frame phase                      //
  //------------------------------------------------------------------------//

  // At theta = 0 only m = N = 2 survives, so f(0, phi) = c exp(2 i phi):
  // same modulus at every phi, and a phase that turns by N * (pi / 2N) =
  // pi/2 between phi = 0 and phi = pi/2N.
  {
    const auto atZero = exact(0.0, 0.0);
    const auto atQuarter = exact(0.0, pi / (2 * N));
    std::cout << "\nAt the north pole the value turns with the frame:\n"
              << "  |f(0, 0)|            = " << std::abs(atZero) << "\n"
              << "  |f(0, pi/2N)|        = " << std::abs(atQuarter)
              << "   (the same)\n"
              << "  arg difference       = "
              << std::abs(std::arg(atQuarter / atZero)) << "   (N pi / 2N)\n";
    if (!(std::abs(std::abs(atZero) - std::abs(atQuarter)) < 1e-12)) return 1;
  }

  // Off the sphere in theta is refused with std::invalid_argument; phi is
  // periodic and is reduced modulo 2 pi, which is exact.
  try {
    exact(-0.1, 0.0);
    std::cout << "\nA colatitude below zero should have thrown\n";
    return 1;
  } catch (const std::invalid_argument&) {
    std::cout << "\nA colatitude outside [0, pi] is refused; a longitude is "
                 "wrapped:\n"
              << "  |f(1, 0.3) - f(1, 0.3 + 2 pi)| = "
              << std::abs(exact(1.0, 0.3) - exact(1.0, 0.3 + 2 * pi)) << "\n";
  }

  //------------------------------------------------------------------------//
  //                        Remeshing is one line                            //
  //------------------------------------------------------------------------//

  // The interpolant is a callable of (theta, phi), which is what a field
  // constructor takes (example 02), so building the field on another grid
  // from it is a remesh. Interpolate(field) with no scheme is spectral.
  {
    auto finer = Grid(band + 6, 2);
    const auto moved = SpinField<N, Grid>(finer, Interpolate(field));

    // The same field, evaluated on the finer grid directly from its
    // coefficients, as a check.
    auto there = SpinExpansion<N, Grid>(finer, band);
    for (auto l : expansion.Degrees())
      for (auto m : expansion.Orders(l)) there[l, m] = expansion[l, m];
    const auto want = Evaluate(there);

    auto worst = Real{0};
    for (auto iTheta : finer.CoLatitudeIndices())
      for (auto iPhi : finer.LongitudeIndices())
        worst =
            std::max(worst, std::abs(moved[iTheta, iPhi] - want[iTheta, iPhi]));
    std::cout << "\nRemeshed onto a grid of degree " << finer.MaxDegree()
              << ", against the exact answer there: " << worst << "\n";
    if (!(worst < 1e-12)) return 1;
  }

  //------------------------------------------------------------------------//
  //                          The cheap schemes                              //
  //------------------------------------------------------------------------//

#ifdef GSHTRANS_HAVE_INTERPOLATION
  {
    const auto bilinear = Interpolate(field, Scheme::Bilinear());
    const auto bicubic = Interpolate(field, Scheme::Bicubic());

    std::cout << "\nAgainst the reference, at a point off the grid:\n"
              << "  bilinear error = "
              << std::abs(bilinear(1.1, 2.2) - exact(1.1, 2.2)) << "\n"
              << "  bicubic error  = "
              << std::abs(bicubic(1.1, 2.2) - exact(1.1, 2.2)) << "\n";

    // A local scheme wants an oversampled grid. At the band limit it has
    // barely enough samples to see the field at all; give it room and the
    // error falls at its own order -- second for bilinear, fourth for
    // bicubic. ForBand(band, nMax, oversampling) is how you ask for that
    // room: a grid of degree 8 * band for a field of band `band`.
    //
    // The third argument of Interpolate is the degree at which the field is
    // expanded -- for the polar rows of a local scheme, or for the whole sum
    // of the spectral one. It defaults to the grid's own degree; here the
    // field's band is enough.
    auto roomy = Grid::ForBand(band, 2, 8);
    auto same = SpinExpansion<N, Grid>(roomy, band);
    for (auto l : expansion.Degrees())
      for (auto m : expansion.Orders(l)) same[l, m] = expansion[l, m];
    const auto sampled = Evaluate(same);

    const auto better = Interpolate(sampled, Scheme::Bicubic(), band);
    const auto reference = Interpolate(sampled, Scheme::Spectral(), band);
    std::cout << "  bicubic on an 8x oversampled grid = "
              << std::abs(better(1.1, 2.2) - reference(1.1, 2.2)) << "\n";
  }
#else
  std::cout << "\nBuilt without Interpolation, so Scheme::Bilinear() and\n"
               "Scheme::Bicubic() do not exist. Scheme::Spectral() always\n"
               "does, and nothing above is withdrawn.\n";
#endif

  std::cout
      << "\nWhich to use: Spectral for a handful of points and as the\n"
         "reference; a local scheme on an oversampled grid for many.\n"
         "One whole-grid remesh costs about what fifty spectral point\n"
         "evaluations cost, so past a few dozen scattered points it is\n"
         "cheaper to transform onto a finer grid and interpolate there.\n";

  return 0;
}
