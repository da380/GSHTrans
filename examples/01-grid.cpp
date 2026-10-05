// 01 -- Grids
//
// The grid is where everything starts: it fixes the point set, the quadrature
// and the range of upper indices a field can carry.
//
// What this shows
//   How to build a GaussLegendreGrid, what its point set looks like, that its
//   quadrature integrates over the sphere, how ForBand asks for a grid with
//   headroom above the band of the fields, and that a grid is a cheap
//   value-semantic handle with an identity.
//
// Assumes
//   Nothing from the library. Some familiarity with spherical harmonics
//   Y_lm(theta, phi) of degree l and order m, |m| <= l.
//
// Introduced
//   GaussLegendreGrid<Real, MRange, NRange>, Grid(lMax, nMax), Grid::ForBand,
//   MaxDegree, NumberOfCoLatitudes, NumberOfLongitudes, FieldSize, Weights,
//   Identity.
//
// Output
//   The grid's sizes, the sum of its quadrature weights next to 4 pi, the
//   degrees ForBand chose, and which of two grids share an identity.
//
// The quadrature and the choice of nPhi are set out in
// docs/gshtrans-reference.tex, section "Quadrature".

#include <GSHTrans/GSHTrans.hpp>
#include <iostream>
#include <numbers>

int main() {
  using namespace GSHTrans;

  using Real = double;
  // The three template parameters are the precision (float, double or long
  // double), which orders m are stored and which upper indices N are covered.
  // All for both is the general choice: every order -l..l, and every upper
  // index -nMax..nMax. The narrower options (NonNegative orders, which serve
  // real scalars alone; NonNegative or Single upper indices) are specialised
  // and not needed in this series.
  using Grid = GaussLegendreGrid<Real, All, All>;

  // Maximum degree, and the largest upper index the grid will be asked for.
  // A rank-p tensor needs upper indices up to p, so nMax = 2 covers rank 2.
  // Construction builds the table of Wigner values for every degree and every
  // covered upper index, which is the expensive part of a grid.
  auto grid = Grid(16, 2);

  // Gauss-Legendre in colatitude, equally spaced in longitude. nPhi is the
  // smallest fast FFT length that resolves every order |m| <= lMax, which
  // needs 2 lMax + 1 points and not 2 lMax -- at 2 lMax the orders m = +lMax
  // and m = -lMax are the same discrete mode.
  //
  // The lMax + 1 colatitudes lie strictly inside (0, pi), symmetric about the
  // equator, so neither pole is a grid point. Samples are stored
  // colatitude-major with longitude fastest: point (iTheta, iPhi) sits at
  // flat index iTheta * nPhi + iPhi, and FieldSize is nTheta * nPhi.
  std::cout << "lMax   " << grid.MaxDegree() << "\n"
            << "nTheta " << grid.NumberOfCoLatitudes() << "  (lMax + 1)\n"
            << "nPhi   " << grid.NumberOfLongitudes()
            << "  (least fast FFT length >= 2 lMax + 1)\n"
            << "points " << grid.FieldSize() << "\n\n";

  // The quadrature integrates a band-limited function exactly. Summing the
  // weights integrates the constant 1, which is the area of the sphere.
  // Weights() gives one weight per point, in storage order: the Gauss-Legendre
  // weight of its colatitude times the uniform longitude step 2 pi / nPhi.
  auto area = Real{0};
  for (auto w : grid.Weights()) area += w;
  std::cout << "sum of weights " << area << "   4 pi = " << 4 * std::numbers::pi
            << "\n\n";

  // A product of two band-limited fields has twice the band, so a grid that
  // only just resolves the fields cannot hold, or transform exactly, their
  // product. ForBand asks for the headroom rather than leaving a caller to
  // compute a degree and hope: the argument is the band of the *fields*, and
  // the factor is how much room to leave.
  // The grid resolves degrees up to ceil(oversampling * lBand); 2.0 is exact
  // for a single product, and the local interpolation schemes of example 21
  // want more still.
  auto exact = Grid::ForBand(16, 2, 2.0);
  auto dealiased = Grid::ForBand(16, 2, 1.5);  // the 3/2 rule
  std::cout << "ForBand(16, 2, 2.0)  -> lMax " << exact.MaxDegree() << "\n"
            << "ForBand(16, 2, 1.5)  -> lMax " << dealiased.MaxDegree()
            << "\n\n";

  // Grids are value-semantic handles over shared, immutable state, so copying
  // one is a pointer copy and not a copy of its Wigner table -- which at
  // lMax = 256, nMax = 2 would be about 648 MB. Two grids are the same grid
  // when their identities agree; equal parameters are not enough. Fields
  // remember their grid, and combining fields whose grids have different
  // identities throws std::invalid_argument, so build a grid once and pass
  // it (or copies of it) around.
  auto copy = grid;
  auto twin = Grid(16, 2);
  std::cout << "copy shares the original: "
            << (copy.Identity() == grid.Identity() ? "yes" : "no") << "\n"
            << "an identical grid does not: "
            << (twin.Identity() == grid.Identity() ? "yes" : "no") << "\n";
}
