// 06 -- The transform
//
// Between the spatial and spectral representations. The one thing worth
// understanding before using it: transforming an arbitrary field is a
// *projection*, not a round trip.
//
// The mathematics, briefly. A field f of upper index N expands as
//
//   f(theta, phi) = sum_{l >= |N|} sum_{m = -l..l} f^N_lm Y^N_lm(theta, phi),
//   Y^N_lm = sqrt((2l + 1) / (4 pi)) d^l_Nm(theta) exp(i m phi),
//
// with d^l_Nm the Wigner d-function, upper index first (Dahlen & Tromp's
// generalised Legendre function P^N_lm). The Y^N_lm are orthonormal, so
// f^N_lm is the integral of conj(Y^N_lm) f over the sphere. The forward
// transform computes that integral with the grid's quadrature -- an FFT in
// longitude at each colatitude, then a Gauss-Legendre sum against the stored
// d-functions -- for every l up to the requested lMax; the inverse transform
// sums the series back onto the grid. See docs/gshtrans-reference.tex,
// sections "Generalized spherical harmonics" and "The transform".
//
// What this shows
//   Expand and Evaluate, the coefficient range, why a round trip of an
//   arbitrary field is not an identity, and the reduced storage of a real
//   scalar.
//
// Assumes
//   Examples 02-04.
//
// Introduced
//   Expand(field, lMax), Evaluate(expansion), SpinExpansion (MinDegree,
//   MaxDegree, Size, operator[](l, m)).
//
// Output
//   The degree range and coefficient count of a spin-2 expansion, the error
//   of a first round trip (order one: a projection) and of a second (rounding
//   level), and the coefficient counts of a real and a complex scalar.
//
// The spectral side has more to it than this -- indexing by degree and order
// and the operators that act there -- which example 12 and those after it
// take up.

#include <GSHTrans/GSHTrans.hpp>
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
  // expression -- so nothing needs materialising first. The degree may be at
  // most the grid's, and an optional Execution policy (example 07) is the
  // third argument. The result is a SpinExpansion carrying the same upper
  // index in its type.
  auto e = Expand(f, lMax);
  static_assert(decltype(e)::UpperIndex == 2);

  // The coefficients start at degree |N|. Below that no harmonic of that
  // upper index exists, so there is nothing to hold: d^l_Nm vanishes
  // identically for l < |N|. They are stored degree by degree, each degree
  // holding its orders m = -l..l, and are read as e[l, m]: 285 here, the
  // sum of 2l + 1 over l = 2..16.
  std::cout << "degrees " << e.MinDegree() << " .. " << e.MaxDegree() << ",  "
            << e.Size() << " coefficients\n"
            << "f^2_{2,0} = " << (e[2, 0]) << "\n\n";

  // Now the projection. The field above has content at degrees 0 and 1 that
  // a spin-2 field cannot carry, and the transform discards it -- so
  // evaluating the coefficients does not give back what we started with. It
  // gives the projection onto the span of the Y^2_lm, l = 2..lMax.
  auto projected = Evaluate(e);
  auto drift = Real{0};
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      drift =
          std::max(drift, std::abs(projected[iTheta, iPhi] - f[iTheta, iPhi]));
    }
  }
  std::cout << "arbitrary field, once through: " << drift
            << "   (a projection, not an error)\n";

  // Once the field is band-limited, the round trip is an identity. A field
  // that came out of Evaluate is band-limited by construction.
  auto again = Expand(projected, lMax);
  auto twice = Evaluate(again);
  drift = 0;
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      drift = std::max(drift,
                       std::abs(twice[iTheta, iPhi] - projected[iTheta, iPhi]));
    }
  }
  std::cout << "band-limited field, again:     " << drift << "\n\n";

  // A real scalar field uses the reduced m >= 0 storage, since its negative
  // orders follow from f_{l,-m} = (-1)^m conj(f_{lm}). Half the coefficients
  // for the same information. Expand chooses the real path from the field's
  // value kind; the coefficients themselves are always complex.
  //
  // For comparison, multiplying by a complex scalar promotes the field to
  // ComplexValued, and Materialise stores it as a complex SpinField<0>,
  // which takes every order.
  auto g = SpinField<0, Grid, RealValued>(grid, [](auto theta, auto phi) {
    return std::cos(theta) * std::sin(phi);
  });
  auto real = Expand(g, lMax);
  auto complexOne = Expand(Materialise(g * Complex{1.0, 0.0}), lMax);
  std::cout << "real scalar    " << real.Size() << " coefficients\n"
            << "complex scalar " << complexOne.Size() << "\n";
}
