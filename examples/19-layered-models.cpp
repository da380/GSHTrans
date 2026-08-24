// 19 -- Layered models: elements, interfaces, and remeshing
//
// A radial grid has always been able to hold a repeated radius. What it could
// not do was say what one *meant*: at the core-mantle boundary a repetition is
// a material interface with two sides, and in a typo it is a mistake, and
// nothing could tell them apart.
//
// `RadialGrid::WithElements` is the difference. It carries which radii belong
// to which element -- the smallest fact that distinguishes a discretisation
// from a list of numbers, and the one thing here that more than one facility
// needs and none can infer.

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <span>
#include <vector>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  using Int = std::ptrdiff_t;

  constexpr auto lMax = Int{8};
  auto grid = Grid(lMax, 2);

  std::cout << std::scientific << std::setprecision(3);

  //------------------------------------------------------------------------//
  // A model with an interface
  //------------------------------------------------------------------------//

  // Two elements meeting at r = 0.8, which is held twice: once as the top of
  // the lower element and once as the bottom of the upper one. The blocks are
  // disjoint, so every radius belongs to exactly one element -- which is what
  // makes a derivative at the interface well defined without anyone having to
  // choose between the two sides.
  const auto radii = std::vector<Real>{0.40, 0.60, 0.80, 0.80, 1.00, 1.20};
  const auto starts = std::vector<Int>{0, 3, 6};
  const auto mesh = RadialGrid<Real>::WithElements(radii, starts);

  std::cout << "a two-element mesh\n"
            << "  radii     " << mesh.NumberOfRadii() << '\n'
            << "  elements  " << mesh.ElementCount() << '\n'
            << "  breaks    " << mesh.Breakpoint(0) << ", "
            << mesh.Breakpoint(1) << ", " << mesh.Breakpoint(2) << '\n'
            << "  index 2 is in element " << mesh.ElementOf(2)
            << ", index 3 in element " << mesh.ElementOf(3) << "\n\n";

  // The partition has to *be* one, and each way of not being one has its own
  // message: not covering the radii, an element of a single node, a repeated
  // radius inside an element -- which is a mistake rather than an interface --
  // or a gap where the field would be undefined.
  try {
    (void)RadialGrid<Real>::WithElements(
        std::vector<Real>{0.4, 0.6, 0.6, 0.8, 1.0, 1.2}, {0, 3, 6});
  } catch (const std::invalid_argument& e) {
    std::cout << "a repeat inside an element is refused:\n  " << e.what()
              << "\n\n";
  }

  //------------------------------------------------------------------------//
  // The derivative is block-diagonal, and two-valued at the interface
  //------------------------------------------------------------------------//

  const auto element = ElementDerivative<Real>(mesh);

  // A field with a jump: slope +1 below the interface, slope -3 above it.
  auto profile = std::vector<Real>{};
  for (std::size_t i = 0; i < radii.size(); i++) {
    profile.push_back(i < 3 ? radii[i] : 10.0 - 3.0 * radii[i]);
  }

  auto slope = std::vector<Real>(radii.size());
  element(std::span<const Real>(profile), std::span<Real>(slope));

  std::cout << "d/dr across an interface\n"
            << "     r        value      slope\n";
  for (std::size_t i = 0; i < radii.size(); i++) {
    std::cout << "  " << radii[i] << "   " << profile[i] << "   " << slope[i]
              << (i == 2 ? "   <- from below\n"
                         : (i == 3 ? "   <- from above\n" : "\n"));
  }
  std::cout
      << "  both sides are reported, each at its own index. That is what\n"
         "  a discontinuity is, and it is why the blocks are disjoint.\n\n";

  // Nothing in one element can move the answer in another, which is what
  // block-diagonal means and what a global operator on the same radii would
  // not give.
  auto disturbed = profile;
  disturbed[4] += 1.0;
  auto after = std::vector<Real>(radii.size());
  element(std::span<const Real>(disturbed), std::span<Real>(after));
  std::cout << "disturbing the upper element moves the lower one by "
            << std::abs(after[0] - slope[0]) << "\n\n";

#ifdef GSHTRANS_HAVE_INTERPOLATION
  //------------------------------------------------------------------------//
  // Splines stop refusing, and fit per element
  //------------------------------------------------------------------------//

  // Without a partition a spline through a repeated radius is refused: it does
  // not know whether the repetition is meaningful. With one, it fits each
  // element separately and never spans the interface.
  const auto spline = SplineDerivative<Real>(mesh);
  std::cout << "the spline derivative on this mesh fits " << spline.PieceCount()
            << " pieces\n\n";

  //------------------------------------------------------------------------//
  // Remeshing, which must not cross an interface either
  //------------------------------------------------------------------------//

  auto f = LayeredSpinField<0, Grid, ComplexValued>(mesh, grid);
  for (auto i : mesh.RadiusIndices()) {
    const auto r = mesh.Radius(i);
    const auto value = i < 3 ? Complex{r, 0} : Complex{10.0 - 3.0 * r, 0};
    auto layer = f.Slice(i);
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) layer[iTheta, iPhi] = value;
    }
  }

  // Targets on both sides of the break, none of them a source node. Each is
  // answered from the piece that owns it, so the two straight lines are
  // reproduced exactly rather than smeared into one curve through the jump.
  const auto onto =
      RadialGrid<Real>(std::vector<Real>{0.5, 0.7, 0.8, 0.9, 1.1});
  const auto moved = Resample(f, onto, RadialInterpolation::CubicSpline());

  std::cout << "resampling across the interface\n"
            << "     r        value      exact\n";
  for (auto i : onto.RadiusIndices()) {
    const auto r = onto.Radius(i);
    const auto exact = r < 0.8 ? r : 10.0 - 3.0 * r;
    std::cout << "  " << r << "   " << (moved.Slice(i)[1, 2]).real() << "   "
              << exact << (i == 2 ? "   <- on the break\n" : "\n");
  }
  std::cout << "  a radius landing exactly on a breakpoint is answered from\n"
               "  above, right-continuously -- the convention Interpolation's\n"
               "  Piecewise uses, so the two libraries cannot disagree about\n"
               "  which side of the core-mantle boundary a query belongs to.\n";
#endif

  return 0;
}
