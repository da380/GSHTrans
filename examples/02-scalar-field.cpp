// 02 -- A scalar field, and integration
//
// The unit of the field algebra is a single field of definite upper index.
// At upper index zero that is an ordinary scalar field, and it may be real.
//
// What this shows
//   How to sample a function of (theta, phi) onto a grid as a SpinField, how
//   to integrate over the sphere with the grid's quadrature, and how to read a
//   single sample.
//
// Assumes
//   Example 01 (the grid).
//
// Introduced
//   SpinField<N, Grid, Value>, RealValued, Integrate, operator[](iTheta, iPhi).
//
// Output
//   Integrals of two degree-1 functions (zero to rounding), of their product
//   (zero: they are orthogonal) and of f^2 next to its exact value 4 pi / 3.

#include <GSHTrans/GSHTrans.hpp>
#include <cmath>
#include <iostream>
#include <numbers>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Grid = GaussLegendreGrid<Real, All, All>;

  auto grid = Grid(24, 2);

  // Built from a function of (theta, phi), theta the colatitude in [0, pi]
  // and phi the longitude in [0, 2 pi), called once at every grid point.
  //
  // The template arguments are the upper index N (fixed at compile time),
  // the grid type, and the value kind, which defaults to ComplexValued.
  // RealValued is available only at upper index zero: real-valuedness is not
  // preserved by a rotation of the local frame, so no component of any tensor
  // can have it elsewhere (example 03 says more).
  auto f = SpinField<0, Grid, RealValued>(
      grid, [](auto theta, auto) { return std::cos(theta); });

  auto g = SpinField<0, Grid, RealValued>(grid, [](auto theta, auto phi) {
    return std::sin(theta) * std::cos(phi);
  });

  // Integrate is the only reduction this layer offers, and it exists only at
  // upper index zero: the integral of a field over the sphere vanishes
  // identically unless N = 0. It is the grid's quadrature, the sum over
  // points of weight times value, so it is exact for an integrand of band up
  // to 2 lMax -- a product of two band-lMax fields, for instance. It returns
  // Real for a
  // real field and std::complex<Real> for a complex one.
  std::cout << "integral of cos(theta)        " << Integrate(f) << "\n"
            << "integral of sin cos phi       " << Integrate(g) << "\n";

  // The two are orthogonal, and each has a norm the quadrature gets exactly.
  // These are Y_10 and a combination of Y_11, Y_1-1 up to normalisation.
  // f * g is not computed here: it is a lazy expression that Integrate
  // evaluates point by point (example 04).
  std::cout << "integral of the product       " << Integrate(f * g) << "\n"
            << "integral of f^2               " << Integrate(f * f)
            << "   4 pi / 3 = " << 4 * std::numbers::pi / 3 << "\n\n";

  // Point access is by (iTheta, iPhi), longitude fastest, using C++23's
  // multidimensional subscript. Reading is by value on every node in the
  // library, so an expression and a field are interchangeable; writing needs
  // storage, so it is offered on fields and views alone. Size() is the number
  // of samples, the grid's FieldSize.
  std::cout << "f at the first point " << (f[0, 0]) << "\n"
            << "field size           " << f.Size() << "\n";
}
