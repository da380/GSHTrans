// 02 -- A scalar field, and integration
//
// The unit of the field algebra is a single field of definite upper index.
// At upper index zero that is an ordinary scalar field, and it may be real.

#include <GSHTrans/All>
#include <cmath>
#include <iostream>
#include <numbers>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Grid = GaussLegendreGrid<Real, All, All>;

  auto grid = Grid(24, 2);

  // Built from a function of (theta, phi). RealValued is available only at
  // upper index zero: real-valuedness is not preserved by a rotation of the
  // local frame, so no component of any tensor can have it elsewhere.
  auto f = SpinField<0, Grid, RealValued>(
      grid, [](auto theta, auto) { return std::cos(theta); });

  auto g = SpinField<0, Grid, RealValued>(grid, [](auto theta, auto phi) {
    return std::sin(theta) * std::cos(phi);
  });

  // Integrate is the only reduction this layer offers, and it exists only at
  // upper index zero: the integral of a field over the sphere vanishes
  // identically unless N = 0.
  std::cout << "integral of cos(theta)        " << Integrate(f) << "\n"
            << "integral of sin cos phi       " << Integrate(g) << "\n";

  // The two are orthogonal, and each has a norm the quadrature gets exactly.
  // These are Y_10 and a combination of Y_11, Y_1-1 up to normalisation.
  std::cout << "integral of the product       " << Integrate(f * g) << "\n"
            << "integral of f^2               " << Integrate(f * f)
            << "   4 pi / 3 = " << 4 * std::numbers::pi / 3 << "\n\n";

  // Point access is by (iTheta, iPhi), longitude fastest. Reading is by
  // value on every node in the library, so an expression and a field are
  // interchangeable; writing needs storage, so it is offered on fields and
  // views alone.
  std::cout << "f at the first point " << (f[0, 0]) << "\n"
            << "field size           " << f.Size() << "\n";

  FFTWpp::CleanUp();
}
