#include <gtest/gtest.h>

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

// Preconditions, in every build mode.
//
// Indexing.h sets the policy: the index classes assert, because they sit in
// inner loops, and "the layers above validate and throw in every build mode".
// This file is where that second half is held to. Each case here once did
// something other than throw -- an assert that a Release build compiles out,
// leaving a heap overwrite or a table of the wrong values, or an assert that a
// Debug build reaches before the documented exception. Benchmarks and
// production both run with NDEBUG, so an assert is not a check that a user of
// this library ever sees.
//
// The suite is built in both modes, and these tests must pass in both: that
// is the point of them.

using namespace GSHTrans;

namespace {

using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

// A grid of someone else's making, which is what the base class is for. The
// two-point Gauss-Legendre rule, with as many longitudes as it is asked for.
class TwoPointGrid : public SphericalGrid<Real, All, All> {
 public:
  TwoPointGrid(Int lMax, Int nPhi)
      : SphericalGrid<Real, All, All>(
            lMax, 0,
            std::vector<Real>{std::acos(1 / std::sqrt(Real{3})),
                              std::acos(-1 / std::sqrt(Real{3}))},
            std::vector<Real>{1, 1}, FFTWpp::Estimate, Chunking::Automatic(),
            WignerValues::Stored(), TransformKernel::Loop(), nPhi) {}
};

}  // namespace

//--------------------------------------------------------------------------//
//                                   Grids                                   //
//--------------------------------------------------------------------------//

TEST(Preconditions, AGridRefusesDegreesAndUpperIndicesThatMakeNoSense) {
  EXPECT_THROW(Grid(-1, 0, FFTWpp::Estimate), std::invalid_argument);
  EXPECT_THROW(Grid(2, 3, FFTWpp::Estimate), std::invalid_argument);
  EXPECT_THROW(Grid(2, -3, FFTWpp::Estimate), std::invalid_argument);
  EXPECT_NO_THROW(Grid(2, 2, FFTWpp::Estimate));
  EXPECT_NO_THROW(Grid(0, 0, FFTWpp::Estimate));
}

TEST(Preconditions, AGridRefusesAPlannerFlagThatCannotPlan) {
  // WisdomOnly makes FFTW return no plan at all unless wisdom is already
  // there, and the plans are made lazily, inside the transform.
  EXPECT_THROW(Grid(4, 0, FFTWpp::WisdomOnly), std::invalid_argument);
  const auto grid = Grid(4, 0, FFTWpp::Estimate);
  EXPECT_THROW(grid.With(FFTWpp::WisdomOnly), std::invalid_argument);
}

TEST(Preconditions, ADerivedGridNeedsEnoughLongitudesForItsOrders) {
  // Both kernels index the Fourier output at every order up to lMax, so fewer
  // than 2 lMax + 1 longitudes is a read past its end.
  EXPECT_THROW(TwoPointGrid(1, 2), std::invalid_argument);
  EXPECT_NO_THROW(TwoPointGrid(1, 3));
}

TEST(Preconditions, ADegreeZeroGridNeedNotBeTheOnePointGrid) {
  // lMax = 0 used to mean "one sample, weight two" to the transform, which is
  // what GaussLegendreGrid builds and not what the base class requires. Here
  // there are eight samples, and every one of them counts.
  const auto grid = TwoPointGrid(0, 4);
  ASSERT_EQ(grid.FieldSize(), 8u);
  const auto y00 = std::numbers::inv_sqrtpi_v<Real> / 2;

  auto coefficient = std::vector<Complex>{Complex{3, -1}};
  auto field = std::vector<Complex>(8, Complex{-7, -7});
  grid.InverseTransformation(0, 0, coefficient, field);
  for (auto value : field) {
    EXPECT_NEAR(std::abs(value - Complex{3, -1} * y00), 0, 1e-14);
  }

  auto back = std::vector<Complex>(1);
  grid.ForwardTransformation(0, 0, field, back);
  EXPECT_NEAR(std::abs(back[0] - Complex{3, -1}), 0, 1e-14);
}

//--------------------------------------------------------------------------//
//                                   Wigner                                  //
//--------------------------------------------------------------------------//

TEST(Preconditions, AWignerTableRefusesWhatItCannotHold) {
  using Table = Wigner<Real, All, All, Single>;
  EXPECT_THROW(Table(-1, 0, 0, 1.0), std::invalid_argument);
  EXPECT_THROW(Table(2, -1, 0, 1.0), std::invalid_argument);
  EXPECT_THROW(Table(2, 2, 3, 1.0), std::invalid_argument);
  // A negative upper index means something only when there is exactly one.
  EXPECT_THROW(Table(4, 4, -2, 1.0), std::invalid_argument);
  EXPECT_NO_THROW((Wigner<Real, All, Single, Single>(4, 4, -2, 1.0)));

  EXPECT_THROW(Table(2, 2, 1, -0.1), std::invalid_argument);
  EXPECT_THROW(Table(2, 2, 1, 3.2), std::invalid_argument);
  EXPECT_THROW(Table(2, 2, 1, std::numeric_limits<Real>::quiet_NaN()),
               std::invalid_argument);
  EXPECT_NO_THROW(Table(2, 2, 1, 0.0));
  EXPECT_NO_THROW(Table(2, 2, 1, std::numbers::pi_v<Real>));
}

TEST(Preconditions, AWignerTableRecomputesOnlyAtAsManyAngles) {
  auto table =
      Wigner<Real, All, All, Multiple>(3, 3, 1, std::vector<Real>{0.4, 1.1});
  EXPECT_NO_THROW(table.ReCompute(std::vector<Real>{0.5, 1.2}));
  EXPECT_THROW(table.ReCompute(std::vector<Real>{0.5}), std::invalid_argument);
}

TEST(Preconditions, ANegativeSingleUpperIndexHasItsSquareRoots) {
  // The recursion's square-root tables run to lMax + max(mMax, |n|). Sized
  // with n rather than |n| they were short whenever mMax < -n, and the
  // recursion read past the end. The values must agree with a table wide
  // enough never to have been affected.
  constexpr auto lMax = Int{10};
  constexpr auto n = Int{-2};
  const auto theta = Real{1.0};
  const auto narrow = Wigner<Real, All, Single, Single>(lMax, 0, n, theta);
  const auto wide = Wigner<Real, All, Single, Single>(lMax, lMax, n, theta);
  for (auto l = Int{2}; l <= lMax; l++) {
    EXPECT_EQ(narrow[l][0], wide[l][0]) << "degree " << l;
  }
}

TEST(Preconditions, AnEmptyWignerTableIsEmpty) {
  const auto table = Wigner<Real, All, All, Multiple>();
  EXPECT_EQ(table.MaxDegree(), 0);
  EXPECT_EQ(table.NumberOfAngles(), 0u);
}

//--------------------------------------------------------------------------//
//                                 3-j symbols                               //
//--------------------------------------------------------------------------//

TEST(Preconditions, ThreeJStorageMustBeLargeEnough) {
  // The rows are scattered by index into the caller's storage, so storage
  // that is too small is a write past its end.
  auto right = std::vector<Real>(5 * 7);
  auto small = std::vector<Real>(5 * 7 - 1);
  EXPECT_NO_THROW(FillWigner3jMatrix(2, 3, 3, right));
  EXPECT_THROW(FillWigner3jMatrix(2, 3, 3, small), std::invalid_argument);
  EXPECT_NO_THROW(FillCouplingMatrix(2, 3, 3, right));
  EXPECT_THROW(FillCouplingMatrix(2, 3, 3, small), std::invalid_argument);
}

TEST(Preconditions, ThreeJDegreesAreNonNegative) {
  auto storage = std::vector<Real>(64);
  EXPECT_THROW(FillWigner3jMatrix(-1, 1, 1, storage), std::invalid_argument);
  EXPECT_THROW(FillCouplingMatrix(1, -1, 1, storage), std::invalid_argument);
  EXPECT_THROW(Wigner3jMatrix<Real>(1, 1, -1), std::invalid_argument);
  EXPECT_THROW(Wigner3jStack<Real>(-1, 1), std::invalid_argument);
  EXPECT_THROW(Wigner3jStack<Real>(1, 1, 2, 1), std::invalid_argument);
  // Degrees that fail the triangle condition are not an error: the symbols
  // there are zero, and zero is what is returned.
  EXPECT_NO_THROW(Wigner3jMatrix<Real>(1, 5, 1));
}

//--------------------------------------------------------------------------//
//                                 Expansions                                //
//--------------------------------------------------------------------------//

TEST(Preconditions, AnExpansionViewThrowsBelowItsUpperIndex) {
  // The documented exception, which a Debug build used never to reach: the
  // index block is a member, and its own assert came first.
  auto grid = Grid(4, 2, FFTWpp::Estimate);
  auto storage = std::vector<Complex>(64);
  EXPECT_THROW(
      (SpinExpansionView<2, Grid>(grid, 1, std::span<Complex>(storage))),
      std::invalid_argument);
  EXPECT_THROW((SpinExpansion<2, Grid>(grid, 1)), std::invalid_argument);
}

TEST(Preconditions, ATensorExpansionThrowsBelowItsRank) {
  auto grid = Grid(4, 2, FFTWpp::Estimate);
  using E = TensorExpansion<2, NoSymmetry<2>, ComplexTensor, Grid>;
  EXPECT_THROW(E(grid, 1), std::invalid_argument);
  EXPECT_NO_THROW(E(grid, 2));

  using L = LayeredTensorExpansion<2, NoSymmetry<2>, ComplexTensor, Grid>;
  const auto radial = RadialGrid<Real>(std::vector<Real>{0.5, 1.0});
  EXPECT_THROW(L(radial, grid, 1), std::invalid_argument);
  EXPECT_NO_THROW(L(radial, grid, 2));
}

TEST(Preconditions, TheGradientOfAConstantIsZeroAndNotAnError) {
  // A scalar expansion of degree zero is a perfectly good expansion, and its
  // gradient is a perfectly good vector: zero. It has components at upper
  // index one, so it cannot be held at degree zero, and the result is built at
  // the lowest degree that can hold it instead of being refused.
  auto grid = Grid(4, 2, FFTWpp::Estimate);
  using Scalar = TensorExpansion<0, NoSymmetry<0>, ComplexTensor, Grid>;
  auto constant = Scalar(grid, 0);
  constant.Data()[0] = Complex{2, 1};

  const auto gradient = SurfaceGradient(constant);
  EXPECT_EQ(gradient.MaxDegree(), 1);
  for (auto value : gradient.Data()) EXPECT_EQ(value, Complex{});
}

//--------------------------------------------------------------------------//
//                                Radial grids                               //
//--------------------------------------------------------------------------//

TEST(Preconditions, ARadialGridRefusesARadiusThatIsNotANumber) {
  // A NaN compares false with everything, so it passes a sortedness check and
  // a sign check alike.
  const auto nan = std::numeric_limits<Real>::quiet_NaN();
  EXPECT_THROW(RadialGrid<Real>(std::vector<Real>{0.5, nan, 1.0}),
               std::invalid_argument);
  EXPECT_THROW(RadialGrid<Real>(std::vector<Real>{nan}), std::invalid_argument);
  EXPECT_THROW(RadialGrid<Real>(std::vector<Real>{
                   0.5, std::numeric_limits<Real>::infinity()}),
               std::invalid_argument);
}

TEST(Preconditions, ElementOfAnswersOnlyForARadiusInAnElement) {
  const auto mesh =
      RadialGrid<Real>::WithElements({0.4, 0.6, 0.8, 0.8, 1.0, 1.2}, {0, 3, 6});
  EXPECT_EQ(mesh.ElementOf(0), 0);
  EXPECT_EQ(mesh.ElementOf(2), 0);
  EXPECT_EQ(mesh.ElementOf(3), 1);
  EXPECT_EQ(mesh.ElementOf(5), 1);
  EXPECT_THROW(mesh.ElementOf(-1), std::invalid_argument);
  EXPECT_THROW(mesh.ElementOf(6), std::invalid_argument);

  const auto plain = RadialGrid<Real>(std::vector<Real>{0.4, 0.6, 0.8});
  EXPECT_THROW(plain.ElementOf(1), std::invalid_argument);
}

TEST(Preconditions, ALayeredExpansionRefusesARadiusItDoesNotHave) {
  // The layered field throws for this; its expansion read out of bounds.
  auto grid = Grid(4, 0, FFTWpp::Estimate);
  const auto radial = RadialGrid<Real>(std::vector<Real>{0.5, 1.0});
  auto e = LayeredSpinExpansion<0, Grid>(radial, grid, 4);
  EXPECT_NO_THROW((e[1, 0, 0]));
  EXPECT_THROW((e[2, 0, 0]), std::invalid_argument);
  EXPECT_THROW((e[-1, 0, 0]), std::invalid_argument);
}

//--------------------------------------------------------------------------//
//                                  Indexing                                 //
//--------------------------------------------------------------------------//

TEST(Preconditions, TheIndicesOfATemporaryBlockOutliveIt) {
  // Views are made, used within one expression and let go, so the range of a
  // temporary block is the natural thing to loop over -- and until range-for
  // extends the lifetime of such a temporary, the range must not refer to it.
  // Under the address sanitiser this was a stack-use-after-scope.
  auto count = Int{0};
  auto sumOfOrders = Int{0};
  for (auto [l, m] : GSHIndices<All>(3, 3, 1).Indices()) {
    count++;
    sumOfOrders += m;
    EXPECT_GE(l, 1);
    EXPECT_LE(std::abs(m), l);
  }
  EXPECT_EQ(count, GSHIndices<All>(3, 3, 1).Size());
  EXPECT_EQ(sumOfOrders, 0);
}

//--------------------------------------------------------------------------//
//                    What a grid keeps, and what it lends                   //
//--------------------------------------------------------------------------//

namespace {

template <typename G>
concept LendsItsNodesFromATemporary =
    requires { std::declval<G>().CoLatitudes(); };
template <typename G>
concept LendsItsPointsFromATemporary = requires { std::declval<G>().Points(); };

}  // namespace

TEST(Preconditions, AGridsNodesAreNotLentFromATemporary) {
  // They are views into storage the grid shares with its other handles and
  // frees with the last of them, which for a temporary is the end of the
  // statement. Longitudes are computed and carry what they need.
  static_assert(!LendsItsNodesFromATemporary<Grid>);
  static_assert(!LendsItsPointsFromATemporary<Grid>);
  static_assert(LendsItsNodesFromATemporary<const Grid&>);
  static_assert(requires { std::declval<Grid>().Longitudes(); });

  const auto grid = Grid(4, 0, FFTWpp::Estimate);
  EXPECT_EQ(std::ranges::distance(grid.CoLatitudes()), 5);
}

TEST(Preconditions, AThreadsCachesCanBeGivenBack) {
  // Plans and buffers are kept per thread for the life of the thread, which
  // for an OpenMP worker is the life of the process. Without a way to release
  // them FFTWpp::CleanUp could never be called again after the first
  // transform, since it refuses while a plan is alive.
  const auto grid = Grid(8, 0, FFTWpp::Estimate);
  const auto size = static_cast<std::size_t>(grid.CoefficientSize(8, 0));
  auto given = std::vector<Complex>(size);
  for (std::size_t j = 0; j < size; j++) {
    given[j] = Complex{std::sin(0.3 * static_cast<Real>(j)), 0.25};
  }
  auto field = std::vector<Complex>(grid.FieldSize());
  auto back = std::vector<Complex>(size);

  const auto before = FFTWpp::LivePlanCount();
  grid.InverseTransformation(8, 0, given, field);
  grid.ForwardTransformation(8, 0, field, back);
  const auto during = FFTWpp::LivePlanCount();
  EXPECT_GT(during, 0);

  Grid::ReleaseThreadCaches();
  const auto after = FFTWpp::LivePlanCount();
  EXPECT_LT(after, during);
  EXPECT_LE(after, before);

  // And the next transform simply makes them again.
  auto again = std::vector<Complex>(size);
  grid.InverseTransformation(8, 0, given, field);
  grid.ForwardTransformation(8, 0, field, again);
  for (std::size_t j = 0; j < size; j++) EXPECT_EQ(again[j], back[j]);
}
