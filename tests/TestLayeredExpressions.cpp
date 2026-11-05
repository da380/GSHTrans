// Layered views over a caller's storage, layered expressions, and Expand of a
// layered expression.
//
// What is checked: a view transforms to the owning type's bits and its slices
// alias the caller's array; layered expressions follow the two-dimensional
// index, reality and grid rules slice by slice, evaluate the same on any
// team, and may alias their target; and Expand of an expression, made a
// block of radii at a time, gives the bits of Expand of its stored value.

#include <gtest/gtest.h>

#include <GSHTrans/GSHTrans.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace {

using namespace GSHTrans;
using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

constexpr auto lMax = Int{8};

auto TestGrid() { return Grid(lMax, 2, FFTWpp::Estimate); }

auto TestRadii(Int nR = 5) {
  auto r = std::vector<Real>(static_cast<std::size_t>(nR));
  for (auto i = Int{0}; i < nR; ++i) {
    r[static_cast<std::size_t>(i)] =
        0.5 + 0.5 * static_cast<Real>(i) / static_cast<Real>(nR - 1);
  }
  return RadialGrid<Real>(r);
}

template <typename Span>
void Fill(Span data, Real tag = 0) {
  for (auto i = std::size_t{0}; i < data.size(); ++i) {
    const auto x = static_cast<Real>(i);
    if constexpr (std::same_as<std::remove_const_t<typename Span::element_type>,
                               Real>) {
      data[i] = 1.5 + std::sin(0.23 * x + tag);
    } else {
      data[i] = Complex{std::sin(0.31 * x + tag), std::cos(0.17 * x - tag)};
    }
  }
}

template <typename A, typename B>
void ExpectSameBits(const A& a, const B& b) {
  ASSERT_EQ(a.size(), b.size());
  for (auto i = std::size_t{0}; i < a.size(); ++i) {
    ASSERT_EQ(a[i], b[i]) << "element " << i;
  }
}

}  // namespace

//--------------------------------------------------------------------------//
//                                   Views                                   //
//--------------------------------------------------------------------------//

TEST(LayeredViews, TransformLikeTheOwningTypes) {
  auto grid = TestGrid();
  auto radial = TestRadii();
  auto owning = LayeredSpinField<-2, Grid>(radial, grid);
  Fill(owning.Data());

  // The caller's own arrays, holding the same values.
  auto samples =
      std::vector<Complex>(owning.Data().begin(), owning.Data().end());
  auto view = LayeredSpinFieldView<-2, Grid>(radial, grid, std::span(samples));
  const auto expected = Expand(owning, lMax);

  auto coefficients = std::vector<Complex>(expected.Data().size());
  Expand(view, LayeredSpinExpansionView<-2, Grid>(radial, grid, lMax,
                                                  std::span(coefficients)));
  ExpectSameBits(coefficients, expected.Data());
  ExpectSameBits(Expand(view, lMax).Data(), expected.Data());

  // Back into the caller's array, through a const view of the coefficients.
  std::ranges::fill(samples, Complex{});
  Evaluate(ConstLayeredSpinExpansionView<-2, Grid>(
               radial, grid, lMax, std::span<const Complex>(coefficients)),
           view);
  ExpectSameBits(samples, Evaluate(expected).Data());
}

TEST(LayeredViews, ARealComponentViewTransformsLikeTheOwningType) {
  auto grid = TestGrid();
  auto radial = TestRadii();
  auto values = std::vector<Real>(
      static_cast<std::size_t>(radial.NumberOfRadii() * grid.FieldSize()));
  Fill(std::span(values));
  auto owning = LayeredSpinField<0, Grid, RealValued>(radial, grid);
  std::ranges::copy(values, owning.Data().begin());

  const auto view = ConstLayeredSpinFieldView<0, Grid, RealValued>(
      radial, grid, std::span<const Real>(values));
  ExpectSameBits(Expand(view, lMax).Data(), Expand(owning, lMax).Data());
}

TEST(LayeredViews, SlicesAliasTheCallersStorage) {
  auto grid = TestGrid();
  auto radial = TestRadii(3);
  auto samples =
      std::vector<Complex>(static_cast<std::size_t>(3 * grid.FieldSize()));
  auto view = LayeredSpinFieldView<1, Grid>(radial, grid, std::span(samples));

  view.Slice(2)[1, 3] = Complex{4.0, -2.0};
  EXPECT_EQ(samples[static_cast<std::size_t>(
                2 * grid.FieldSize() + 1 * grid.NumberOfLongitudes() + 3)],
            (Complex{4.0, -2.0}));

  auto coefficients = std::vector<Complex>(
      static_cast<std::size_t>(3 * grid.CoefficientSize(lMax, 1)));
  auto expansion = LayeredSpinExpansionView<1, Grid>(radial, grid, lMax,
                                                     std::span(coefficients));
  expansion[1, 4, -2] = Complex{1.0, 1.0};
  EXPECT_EQ((expansion.Slice(1)[4, -2]), (Complex{1.0, 1.0}));
}

TEST(LayeredViews, RefuseStorageOfTheWrongSize) {
  auto grid = TestGrid();
  auto radial = TestRadii(3);
  auto tooShort =
      std::vector<Complex>(static_cast<std::size_t>(3 * grid.FieldSize() - 1));
  EXPECT_THROW(
      (LayeredSpinFieldView<1, Grid>(radial, grid, std::span(tooShort))),
      std::invalid_argument);
  EXPECT_THROW((LayeredSpinExpansionView<1, Grid>(radial, grid, lMax,
                                                  std::span(tooShort))),
               std::invalid_argument);
}

//--------------------------------------------------------------------------//
//                                Expressions                                //
//--------------------------------------------------------------------------//

namespace {

// Whether `a op b` is offered, so refused combinations can be tested.
template <typename A, typename B>
concept CanAdd = requires(A a, B b) { a + b; };
template <typename A, typename B>
concept CanMultiply = requires(A a, B b) { a* b; };

using Strain = LayeredSpinField<-2, Grid>;
using Weight = LayeredSpinField<0, Grid, RealValued>;
using Flat = SpinField<0, Grid, RealValued>;

}  // namespace

// The two-dimensional rules, at every radius: indices add under a product and
// must agree under a sum, a real factor keeps a complex field complex, and a
// two-dimensional field does not enter except through AtEveryRadius.
TEST(LayeredExpressions, FollowTheTwoDimensionalRules) {
  using Product =
      decltype(std::declval<const Weight&>() * std::declval<const Strain&>());
  static_assert(Product::UpperIndex == -2);
  static_assert(std::same_as<Product::Value, ComplexValued>);
  static_assert(CanMultiply<const Weight&, const Strain&>);
  static_assert(CanAdd<const Strain&, const Strain&>);
  static_assert(!CanAdd<const Strain&, const Weight&>);
  static_assert(!CanMultiply<const Strain&, const Flat&>);
  static_assert(!CanAdd<const Weight&, const Flat&>);
  static_assert(CanMultiply<const Strain&, double>);
  static_assert(!CanMultiply<const Strain&, float>);
}

TEST(LayeredExpressions, EvaluateLikeTheirSlices) {
  auto grid = TestGrid();
  auto radial = TestRadii();
  auto strain = Strain(radial, grid);
  auto weight = Weight(radial, grid);
  Fill(strain.Data(), 1.0);
  Fill(weight.Data(), 2.0);

  const auto product = Materialise(2.0 * (weight * strain) - strain);
  for (auto i : radial.RadiusIndices()) {
    auto expected = SpinField<-2, Grid>(grid);
    expected = 2.0 * (weight.Slice(i) * strain.Slice(i)) - strain.Slice(i);
    ExpectSameBits(product.Slice(i).Data(), expected.Data());
  }

  // A depth-independent weight, made explicit and stored once.
  auto flat = Flat(grid);
  Fill(flat.Data(), 3.0);
  const auto weighted = Materialise(AtEveryRadius(radial, flat) * strain);
  for (auto i : radial.RadiusIndices()) {
    auto expected = SpinField<-2, Grid>(grid);
    expected = flat * strain.Slice(i);
    ExpectSameBits(weighted.Slice(i).Data(), expected.Data());
  }
}

TEST(LayeredExpressions, EvaluateTheSameOnAnyTeam) {
  auto grid = TestGrid();
  auto radial = TestRadii(7);
  auto strain = Strain(radial, grid);
  auto weight = Weight(radial, grid);
  Fill(strain.Data(), 1.0);
  Fill(weight.Data(), 2.0);

  const auto reference = Materialise(conj(weight * strain));
  for (auto threads : {2, 3, 8}) {
    auto threaded = LayeredSpinField<2, Grid>(radial, grid);
    Assign(threaded, conj(weight * strain), Execution::Parallel(threads));
    ExpectSameBits(threaded.Data(), reference.Data());
  }
}

// The target may appear in the expression, as in two dimensions; and a view
// over the caller's storage is a target like any other.
TEST(LayeredExpressions, MayAliasTheirTarget) {
  auto grid = TestGrid();
  auto radial = TestRadii();
  auto strain = Strain(radial, grid);
  Fill(strain.Data());
  const auto before = Materialise(1.0 * strain);

  strain = 3.0 * strain + strain;
  const auto expected = Materialise(4.0 * before);
  for (auto k = std::size_t{0}; k < strain.Data().size(); ++k) {
    EXPECT_NEAR(std::abs(strain.Data()[k] - expected.Data()[k]), 0.0, 1e-14);
  }

  auto storage =
      std::vector<Complex>(before.Data().begin(), before.Data().end());
  auto view = LayeredSpinFieldView<-2, Grid>(radial, grid, std::span(storage));
  Assign(view, 2.0 * view, Execution::Parallel(4));
  ExpectSameBits(storage, Materialise(2.0 * before).Data());
}

TEST(LayeredExpressions, RefuseOperandsOnDifferentGrids) {
  auto grid = TestGrid();
  auto radial = TestRadii();
  auto strain = Strain(radial, grid);
  auto otherRadii = Strain(TestRadii(), grid);
  auto otherGrid = Strain(radial, TestGrid());
  EXPECT_THROW(strain + otherRadii, std::invalid_argument);
  EXPECT_THROW(strain + otherGrid, std::invalid_argument);
  auto target = Strain(TestRadii(), grid);
  EXPECT_THROW(Assign(target, 2.0 * strain), std::invalid_argument);
}

//--------------------------------------------------------------------------//
//                       Expand of a layered expression                      //
//--------------------------------------------------------------------------//

// Made a block of radii at a time, and still the bits of one call over the
// stored value: the blocks are multiples of the transform's chunk, so the
// chunks fall where that call's would. Enough radii for several blocks, under
// the default chunking and a small fixed one, sequentially and threaded.
TEST(LayeredExpressions, ExpandMatchesExpandOfTheStoredValue) {
  auto radial = TestRadii(150);
  for (auto chunking : {Chunking::Automatic(), Chunking::Fixed(5)}) {
    auto grid = Grid(lMax, 2, FFTWpp::Estimate, chunking);
    auto strain = Strain(radial, grid);
    auto weight = Weight(radial, grid);
    Fill(strain.Data(), 1.0);
    Fill(weight.Data(), 2.0);

    for (auto policy : {Execution::Sequential(), Execution::Parallel(3)}) {
      const auto product = weight * strain;
      const auto expected = Expand(Materialise(product), lMax, policy);
      auto out = LayeredSpinExpansion<-2, Grid>(radial, grid, lMax);
      Expand(product, out, policy);
      ExpectSameBits(out.Data(), expected.Data());
      ExpectSameBits(Expand(product, lMax, policy).Data(), expected.Data());

      // A real expression goes through the real transform.
      const auto squared = weight * weight;
      ExpectSameBits(Expand(squared, lMax, policy).Data(),
                     Expand(Materialise(squared), lMax, policy).Data());
    }
  }
}

// What sl3dpp's pass is, on its own arrays: evaluate the strain into the
// caller's grid array, weight it, and expand the product into the caller's
// coefficient array, with nothing stored but those arrays.
TEST(LayeredExpressions, APassOnTheCallersArrays) {
  auto grid = TestGrid();
  auto radial = TestRadii(9);
  const auto nR = static_cast<std::size_t>(radial.NumberOfRadii());
  const auto coefficientSize =
      static_cast<std::size_t>(grid.CoefficientSize(lMax, -2));

  auto spectral = std::vector<Complex>(nR * coefficientSize);
  auto gridValues =
      std::vector<Complex>(nR * static_cast<std::size_t>(grid.FieldSize()));
  auto weights =
      std::vector<Real>(nR * static_cast<std::size_t>(grid.FieldSize()));
  auto analysed = std::vector<Complex>(nR * coefficientSize);
  Fill(std::span(spectral), 1.0);
  Fill(std::span(weights), 2.0);

  const auto policy = Execution::Parallel(4);
  auto strain =
      LayeredSpinFieldView<-2, Grid>(radial, grid, std::span(gridValues));
  const auto weight = ConstLayeredSpinFieldView<0, Grid, RealValued>(
      radial, grid, std::span<const Real>(weights));
  Evaluate(ConstLayeredSpinExpansionView<-2, Grid>(
               radial, grid, lMax, std::span<const Complex>(spectral)),
           strain, policy);
  Expand(weight * strain,
         LayeredSpinExpansionView<-2, Grid>(radial, grid, lMax,
                                            std::span(analysed)),
         policy);

  // The same pass through owning stacks.
  auto coefficients = LayeredSpinExpansion<-2, Grid>(radial, grid, lMax);
  std::ranges::copy(spectral, coefficients.Data().begin());
  auto weightStack = Weight(radial, grid);
  std::ranges::copy(weights, weightStack.Data().begin());
  const auto field = Evaluate(coefficients, policy);
  const auto expected = Expand(Materialise(weightStack * field), lMax, policy);
  ExpectSameBits(analysed, expected.Data());
}
