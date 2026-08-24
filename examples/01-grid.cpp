// 01 -- Grids
//
// The grid is where everything starts: it fixes the point set, the quadrature
// and the range of upper indices a field can carry.

#include <GSHTrans/All>
#include <iostream>
#include <numbers>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Grid = GaussLegendreGrid<Real, All, All>;

  // Maximum degree, and the largest upper index the grid will be asked for.
  // A rank-p tensor needs upper indices up to p, so nMax = 2 covers rank 2.
  auto grid = Grid(16, 2);

  // Gauss-Legendre in colatitude, equally spaced in longitude. nPhi is the
  // smallest fast FFT length that resolves every order |m| <= lMax, which
  // needs 2 lMax + 1 points and not 2 lMax -- at 2 lMax the orders m = +lMax
  // and m = -lMax are the same discrete mode.
  std::cout << "lMax   " << grid.MaxDegree() << "\n"
            << "nTheta " << grid.NumberOfCoLatitudes() << "  (lMax + 1)\n"
            << "nPhi   " << grid.NumberOfLongitudes()
            << "  (least fast FFT length >= 2 lMax + 1)\n"
            << "points " << grid.FieldSize() << "\n\n";

  // The quadrature integrates a band-limited function exactly. Summing the
  // weights integrates the constant 1, which is the area of the sphere.
  auto area = Real{0};
  for (auto w : grid.Weights()) area += w;
  std::cout << "sum of weights " << area << "   4 pi = " << 4 * std::numbers::pi
            << "\n\n";

  // A product of two band-limited fields has twice the band, and integrating
  // |f|^2 integrates a degree-2l quantity. ForBand asks for the headroom
  // rather than leaving a caller to compute a degree and hope: the argument
  // is the band of the *fields*, and the factor is how much room to leave.
  auto exact = Grid::ForBand(16, 2, 2.0);
  auto dealiased = Grid::ForBand(16, 2, 1.5);  // the 3/2 rule
  std::cout << "ForBand(16, 2, 2.0)  -> lMax " << exact.MaxDegree() << "\n"
            << "ForBand(16, 2, 1.5)  -> lMax " << dealiased.MaxDegree()
            << "\n\n";

  // Grids are value-semantic handles over shared, immutable state, so copying
  // one is a pointer copy and not a copy of its Wigner table -- which at
  // lMax = 256 would be 648 MB. Two grids are the same grid when their
  // identities agree; equal parameters are not enough.
  auto copy = grid;
  auto twin = Grid(16, 2);
  std::cout << "copy shares the original: "
            << (copy.Identity() == grid.Identity() ? "yes" : "no") << "\n"
            << "an identical grid does not: "
            << (twin.Identity() == grid.Identity() ? "yes" : "no") << "\n";
}
