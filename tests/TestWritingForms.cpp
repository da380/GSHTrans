// The writing forms of Expand, Evaluate, ExpandToLines and EvaluateLines: the
// result is written into an object the caller holds, or a view over the
// caller's storage, and nothing is allocated.
//
// What is checked, for the flat and layered spin fields, the flat and layered
// tensors and radial lines: that the writing form gives the returning form's
// bits; that writing twice into the same result overwrites rather than
// accumulates; that results on another grid are refused; and, for the flat
// spin field, that views are read and written in place at their stride and
// leave the caller's other storage alone.

#include <gtest/gtest.h>

#include <GSHTrans/GSHTrans.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
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

auto TestRadii() { return RadialGrid<Real>(std::vector<Real>{0.5, 0.7, 1.0}); }

Complex Sample(Int i, Real tag = 0) {
  const auto x = static_cast<Real>(i);
  return {std::sin(0.31 * x + tag), std::cos(0.17 * x - tag)};
}

template <typename Span>
void Fill(Span data, Real tag = 0) {
  for (auto i = std::size_t{0}; i < data.size(); ++i) {
    if constexpr (std::same_as<typename Span::value_type, Real>) {
      data[i] = Sample(static_cast<Int>(i), tag).real();
    } else {
      data[i] = Sample(static_cast<Int>(i), tag);
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

TEST(WritingForms, FlatSpinFieldMatchesTheReturningForm) {
  auto grid = TestGrid();
  auto u = SpinField<1, Grid>(grid);
  Fill(u.Data());

  const auto expected = Expand(u, lMax);
  auto out = SpinExpansion<1, Grid>(grid, lMax);
  Expand(u, out);
  ExpectSameBits(out.Data(), expected.Data());

  // Writing again overwrites: the transform owns its output.
  Expand(u, out);
  ExpectSameBits(out.Data(), expected.Data());

  // A lazy expression goes through scratch and gives the same bits as its
  // materialised value does.
  auto twice = SpinField<1, Grid>(grid);
  twice = 2.0 * u;
  auto fromNode = SpinExpansion<1, Grid>(grid, lMax);
  Expand(2.0 * u, fromNode);
  ExpectSameBits(fromNode.Data(), Expand(twice, lMax).Data());

  const auto field = Evaluate(expected);
  auto back = SpinField<1, Grid>(grid);
  Evaluate(expected, back);
  ExpectSameBits(back.Data(), field.Data());
}

TEST(WritingForms, FlatRealFieldMatchesTheReturningForm) {
  auto grid = TestGrid();
  auto u = SpinField<0, Grid, RealValued>(grid);
  Fill(u.Data());

  const auto expected = Expand(u, lMax);
  auto out = SpinExpansion<0, Grid, RealValued>(grid, lMax);
  Expand(u, out);
  ExpectSameBits(out.Data(), expected.Data());

  auto back = SpinField<0, Grid, RealValued>(grid);
  Evaluate(expected, back);
  ExpectSameBits(back.Data(), Evaluate(expected).Data());
}

// The caller's own storage, interleaved with something else: a view at stride
// two is read in place and written in place, and the other half is left as it
// was.
TEST(WritingForms, ViewsAreReadAndWrittenInPlaceAtTheirStride) {
  auto grid = TestGrid();
  const auto size = static_cast<std::size_t>(grid.FieldSize());
  auto u = SpinField<1, Grid>(grid);
  Fill(u.Data());
  const auto expected = Expand(u, lMax);

  auto interleaved = std::vector<Complex>(2 * size);
  for (auto i = std::size_t{0}; i < size; ++i) {
    interleaved[2 * i] = u.Data()[i];
    interleaved[2 * i + 1] = Complex{-7.0, 7.0};
  }
  auto view = SpinFieldView<1, Grid>(grid, std::span(interleaved), 2);

  auto coefficients = std::vector<Complex>(expected.Data().size());
  Expand(view, SpinExpansionView<1, Grid>(grid, lMax, std::span(coefficients)));
  ExpectSameBits(coefficients, expected.Data());

  // Evaluate back into the same interleaved storage, after clearing the
  // samples the view owns.
  for (auto i = std::size_t{0}; i < size; ++i) interleaved[2 * i] = Complex{};
  Evaluate(expected, view);
  const auto field = Evaluate(expected);
  for (auto i = std::size_t{0}; i < size; ++i) {
    ASSERT_EQ(interleaved[2 * i], field.Data()[i]) << "sample " << i;
    ASSERT_EQ(interleaved[2 * i + 1], (Complex{-7.0, 7.0})) << "gap " << i;
  }
}

TEST(WritingForms, ResultsOnAnotherGridAreRefused) {
  auto grid = TestGrid();
  auto other = TestGrid();  // equal parameters, a different grid
  auto u = SpinField<1, Grid>(grid);
  auto onOther = SpinExpansion<1, Grid>(other, lMax);
  EXPECT_THROW(Expand(u, onOther), std::invalid_argument);
  auto fieldOnOther = SpinField<1, Grid>(other);
  EXPECT_THROW(Evaluate(Expand(u, lMax), fieldOnOther), std::invalid_argument);

  auto radial = TestRadii();
  auto stack = LayeredSpinField<1, Grid>(radial, grid);
  auto stackOnOther = LayeredSpinExpansion<1, Grid>(radial, other, lMax);
  EXPECT_THROW(Expand(stack, stackOnOther), std::invalid_argument);
  auto otherRadii = LayeredSpinExpansion<1, Grid>(TestRadii(), grid, lMax);
  EXPECT_THROW(Expand(stack, otherRadii), std::invalid_argument);
}

TEST(WritingForms, LayeredSpinFieldMatchesTheReturningForm) {
  auto grid = TestGrid();
  auto radial = TestRadii();
  auto stack = LayeredSpinField<2, Grid>(radial, grid);
  Fill(stack.Data());

  const auto expected = Expand(stack, lMax, Execution::Parallel(3));
  auto out = LayeredSpinExpansion<2, Grid>(radial, grid, lMax);
  Expand(stack, out, Execution::Parallel(3));
  Expand(stack, out, Execution::Parallel(3));
  ExpectSameBits(out.Data(), expected.Data());

  auto back = LayeredSpinField<2, Grid>(radial, grid);
  Evaluate(expected, back);
  ExpectSameBits(back.Data(), Evaluate(expected).Data());
}

TEST(WritingForms, RadialLinesMatchTheReturningForm) {
  auto grid = TestGrid();
  auto radial = TestRadii();
  auto stack = LayeredSpinField<0, Grid, RealValued>(radial, grid);
  Fill(stack.Data());

  const auto expected = ExpandToLines(stack, lMax);
  auto lines = expected.SameShape();
  ExpandToLines(stack, lMax, lines);
  ExpectSameBits(lines.Data(), expected.Data());

  auto back = LayeredSpinField<0, Grid, RealValued>(radial, grid);
  EvaluateLines(expected, lMax, back);
  ExpectSameBits(back.Data(), EvaluateLines(expected, grid, lMax).Data());

  // Lines for another degree are refused, not misread.
  EXPECT_THROW(ExpandToLines(stack, lMax - 1, lines), std::invalid_argument);
  EXPECT_THROW(EvaluateLines(expected, lMax - 1, back), std::invalid_argument);
}

TEST(WritingForms, FlatTensorMatchesTheReturningForm) {
  using T = TensorField<1, NoSymmetry<1>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  Fill(t.Data());

  const auto expected = Expand(t, lMax);
  auto out = TensorExpansion<1, NoSymmetry<1>, ComplexTensor, Grid>(grid, lMax);
  Expand(t, out);
  ExpectSameBits(out.Data(), expected.Data());

  auto back = T(grid);
  Evaluate(expected, back);
  ExpectSameBits(back.Data(), Evaluate(expected).Data());
}

TEST(WritingForms, LayeredTensorMatchesTheReturningForm) {
  using T = LayeredVectorField<Grid, ComplexTensor>;
  auto grid = TestGrid();
  auto radial = TestRadii();
  auto t = T(radial, grid);
  Fill(t.ComponentStack<-1>().Data(), 1.0);
  Fill(t.ComponentStack<0>().Data(), 2.0);
  Fill(t.ComponentStack<1>().Data(), 3.0);

  const auto expected = Expand(t, lMax);
  auto out = LayeredVectorExpansion<Grid, ComplexTensor>(radial, grid, lMax);
  Expand(t, out);
  ExpectSameBits(out.ComponentStack<-1>().Data(),
                 expected.ComponentStack<-1>().Data());
  ExpectSameBits(out.ComponentStack<0>().Data(),
                 expected.ComponentStack<0>().Data());
  ExpectSameBits(out.ComponentStack<1>().Data(),
                 expected.ComponentStack<1>().Data());

  auto back = T(radial, grid);
  Evaluate(expected, back);
  const auto field = Evaluate(expected);
  ExpectSameBits(back.ComponentStack<1>().Data(),
                 field.ComponentStack<1>().Data());

  auto onOther =
      LayeredVectorExpansion<Grid, ComplexTensor>(TestRadii(), grid, lMax);
  EXPECT_THROW(Expand(t, onOther), std::invalid_argument);
}
