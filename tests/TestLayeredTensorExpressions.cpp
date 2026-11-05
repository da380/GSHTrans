// Tensor expressions over layered tensors.
//
// What is checked: a layered tensor's components are layered expressions,
// derived ones included; the flat tensor operations -- Permute, Transpose,
// TensorProduct, Contract, Symmetrise, and the scale by a scalar field --
// applied to layered tensors give at every radius exactly what they give on
// that radius's flat tensors; evaluation is the same on any team; and Expand
// of a layered tensor expression gives the bits of expanding its stored value,
// a real symmetric target included.

#include <gtest/gtest.h>

#include <GSHTrans/GSHTrans.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <utility>
#include <vector>

namespace {

using namespace GSHTrans;
using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

constexpr auto lMax = Int{8};

auto TestGrid() { return Grid(lMax, 2, FFTWpp::Estimate); }

auto TestRadii(Int nR = 4) {
  auto r = std::vector<Real>(static_cast<std::size_t>(nR));
  for (auto i = Int{0}; i < nR; ++i) {
    r[static_cast<std::size_t>(i)] = 0.5 + 0.1 * static_cast<Real>(i);
  }
  return RadialGrid<Real>(r);
}

using Matrix = LayeredTensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
using Vector = LayeredTensorField<1, NoSymmetry<1>, ComplexTensor, Grid>;
using RealSymmetric = LayeredTensorField<2, Symmetric<2>, RealTensor, Grid>;
using Weight = LayeredSpinField<0, Grid, RealValued>;

// Every stored component of a layered tensor, filled with distinct values.
template <typename Layered>
void Fill(Layered& t, Real tag) {
  auto slot = Real{0};
  LayeredTensorDetails::ForEachStoredSlot(
      t, [&]<auto indices, ComponentConstraint constraint>(auto& stack) {
        auto data = stack.Data();
        for (auto k = std::size_t{0}; k < data.size(); ++k) {
          const auto x = static_cast<Real>(k) + 100 * slot + tag;
          if constexpr (std::same_as<std::remove_cvref_t<decltype(data[k])>,
                                     Real>) {
            data[k] = std::sin(0.13 * x);
          } else {
            data[k] = Complex{std::sin(0.13 * x), std::cos(0.07 * x)};
          }
        }
        slot += 1;
      });
}

template <auto Indices, typename Flat, typename Layered, std::size_t... I>
void CopyComponent(Flat& flat, const Layered& t, Int i,
                   std::index_sequence<I...>) {
  auto target = flat.template Component<Indices[I]...>();
  const auto source = t.template Component<Indices[I]...>(i);
  for (auto iTheta : flat.Grid().CoLatitudeIndices()) {
    for (auto iPhi : flat.Grid().LongitudeIndices()) {
      target[iTheta, iPhi] = source[iTheta, iPhi];
    }
  }
}

// The flat tensor at radius i of a layered one.
template <typename Layered>
auto FlatSlice(const Layered& t, Int i) {
  using L = Layered;
  auto flat = typename L::Flat(t.Grid());
  [&]<std::size_t... S>(std::index_sequence<S...>) {
    (
        [&] {
          constexpr auto f = L::ComponentLayout.flatOfSlot[S];
          constexpr auto indices =
              MultiIndex<L::Rank, typename L::SlotSet>::FromFlat(f).Slots();
          CopyComponent<indices>(
              flat, t, i,
              std::make_index_sequence<static_cast<std::size_t>(L::Rank)>{});
        }(),
        ...);
  }(std::make_index_sequence<static_cast<std::size_t>(L::StoredComponents)>{});
  return flat;
}

template <auto Indices, typename Layered, typename Flat, std::size_t... I>
void ExpectSameComponent(const Layered& layered, Int i, const Flat& flat,
                         std::index_sequence<I...>) {
  const auto a = layered.template Component<Indices[I]...>(i);
  const auto b = flat.template Component<Indices[I]...>();
  for (auto iTheta : flat.Grid().CoLatitudeIndices()) {
    for (auto iPhi : flat.Grid().LongitudeIndices()) {
      ASSERT_EQ((a[iTheta, iPhi]), (b[iTheta, iPhi]))
          << "radius " << i << " at " << iTheta << ", " << iPhi;
    }
  }
}

// Every stored component of a layered result, at every radius, against the
// flat result computed on that radius's slices.
template <typename Layered, typename FlatAt>
void ExpectSameAtEveryRadius(const Layered& layered, FlatAt flatAt) {
  using L = Layered;
  for (auto i : layered.RadiusIndices()) {
    const auto flat = flatAt(i);
    [&]<std::size_t... S>(std::index_sequence<S...>) {
      (
          [&] {
            constexpr auto f = L::ComponentLayout.flatOfSlot[S];
            constexpr auto indices =
                MultiIndex<L::Rank, typename L::SlotSet>::FromFlat(f).Slots();
            ExpectSameComponent<indices>(
                layered, i, flat,
                std::make_index_sequence<static_cast<std::size_t>(L::Rank)>{});
          }(),
          ...);
    }(std::make_index_sequence<static_cast<std::size_t>(
          L::StoredComponents)>{});
  }
}

}  // namespace

TEST(LayeredTensorExpressions, AreToldApartFromFlatOnes) {
  static_assert(LayeredTensorExpr<Matrix>);
  static_assert(!FlatTensorExpr<Matrix>);
  static_assert(FlatTensorExpr<Matrix::Flat>);
  static_assert(!LayeredTensorExpr<Matrix::Flat>);
  static_assert(LayeredTensorExpr<decltype(Contract<0, 1>(
                    std::declval<const Matrix&>()))>);
}

// A derived component -- here the (+, -) of a real symmetric tensor, the
// conjugate of a stored one -- is its representative's, as the slice
// accessor already says.
TEST(LayeredTensorExpressions, ComponentsAreLayeredExpressions) {
  auto grid = TestGrid();
  auto radial = TestRadii();
  auto s = RealSymmetric(radial, grid);
  Fill(s, 1.0);
  const auto plusMinus = Materialise(s.Component<1, -1>());
  const auto zeroZero = Materialise(s.Component<0, 0>());
  for (auto i : radial.RadiusIndices()) {
    const auto expected = s.Component<1, -1>(i);
    const auto real = s.Component<0, 0>(i);
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        ASSERT_EQ((plusMinus.Slice(i)[iTheta, iPhi]), (expected[iTheta, iPhi]));
        ASSERT_EQ((zeroZero.Slice(i)[iTheta, iPhi]), (real[iTheta, iPhi]));
      }
    }
  }
}

TEST(LayeredTensorExpressions, OperationsAgreeWithTheFlatOnesAtEveryRadius) {
  auto grid = TestGrid();
  auto radial = TestRadii();
  auto a = Matrix(radial, grid);
  auto v = Vector(radial, grid);
  auto w = Weight(radial, grid);
  Fill(a, 1.0);
  Fill(v, 2.0);
  for (auto& x : w.Data()) x = 1.5;

  // A contraction of a product: the shape of a kernel.
  ExpectSameAtEveryRadius(
      Materialise(Contract<1, 2>(TensorProduct(a, v))), [&](Int i) {
        const auto fa = FlatSlice(a, i);
        const auto fv = FlatSlice(v, i);
        return Materialise(Contract<1, 2>(TensorProduct(fa, fv)));
      });

  ExpectSameAtEveryRadius(Materialise(Transpose(a)), [&](Int i) {
    return Materialise(Transpose(FlatSlice(a, i)));
  });

  ExpectSameAtEveryRadius(
      Materialise<Symmetric<2>>(Symmetrise<Symmetric<2>>(a)), [&](Int i) {
        return Materialise<Symmetric<2>>(
            Symmetrise<Symmetric<2>>(FlatSlice(a, i)));
      });

  ExpectSameAtEveryRadius(Materialise(w * a), [&](Int i) {
    const auto fa = FlatSlice(a, i);
    return Materialise(w.Slice(i) * fa);
  });
}

TEST(LayeredTensorExpressions, EvaluateTheSameOnAnyTeam) {
  auto grid = TestGrid();
  auto radial = TestRadii(7);
  auto a = Matrix(radial, grid);
  auto v = Vector(radial, grid);
  Fill(a, 1.0);
  Fill(v, 2.0);
  const auto expr = Contract<0, 1>(TensorProduct(a, v));
  const auto reference = Materialise(expr);
  for (auto threads : {2, 3, 8}) {
    const auto threaded = Materialise(expr, Execution::Parallel(threads));
    ExpectSameAtEveryRadius(threaded,
                            [&](Int i) { return FlatSlice(reference, i); });
  }
}

// Expanded a stored component at a time, each through the blocked layered
// Expand: the bits of expanding the stored value. A real symmetric target
// takes the surviving part of each pinned component.
TEST(LayeredTensorExpressions, ExpandMatchesExpandOfTheStoredValue) {
  auto grid = TestGrid();
  auto radial = TestRadii(70);
  auto s = RealSymmetric(radial, grid);
  auto w = Weight(radial, grid);
  Fill(s, 1.0);
  for (auto k = std::size_t{0}; k < w.Data().size(); ++k) {
    w.Data()[k] = 1.0 + 0.1 * std::sin(0.01 * static_cast<Real>(k));
  }

  for (auto policy : {Execution::Sequential(), Execution::Parallel(3)}) {
    const auto expr = w * s;
    const auto stored = Materialise<Symmetric<2>, RealTensor>(expr, policy);
    const auto expected = Expand(stored, lMax, policy);
    const auto got = Expand<Symmetric<2>, RealTensor>(expr, lMax, policy);
    LayeredTensorDetails::ForEachStoredSlot(
        expected, [&]<auto indices, ComponentConstraint>(const auto& stack) {
          const auto& mine = LayeredTensorDetails::StackOf<indices>(
              got, std::make_index_sequence<2>{});
          ASSERT_EQ(stack.Data().size(), mine.Data().size());
          for (auto k = std::size_t{0}; k < stack.Data().size(); ++k) {
            ASSERT_EQ(stack.Data()[k], mine.Data()[k]) << "element " << k;
          }
        });
  }
}

// The flat scale product, for completeness: the weight times each component.
TEST(LayeredTensorExpressions, AFlatTensorScalesByAFlatField) {
  auto grid = TestGrid();
  auto t = Matrix::Flat(grid);
  auto w = SpinField<0, Grid, RealValued>(grid);
  for (auto k = std::size_t{0}; k < t.Data().size(); ++k) {
    t.Data()[k] = Complex{std::sin(0.1 * k), std::cos(0.2 * k)};
  }
  for (auto& x : w.Data()) x = 0.75;
  const auto scaled = Materialise(w * t);
  const auto back = Materialise(t * w);
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      ASSERT_EQ((scaled.Component<-1, 1>()[iTheta, iPhi]),
                0.75 * (t.Component<-1, 1>()[iTheta, iPhi]));
      ASSERT_EQ((back.Component<0, -1>()[iTheta, iPhi]),
                (scaled.Component<0, -1>()[iTheta, iPhi]));
    }
  }
}

// A layered tensor over the caller's arrays, one per stored component: it
// transforms, and is assigned to, exactly as the owning type.
TEST(LayeredTensorExpressions, ViewsWorkOnTheCallersArrays) {
  auto grid = TestGrid();
  auto radial = TestRadii(5);
  auto owning = Vector(radial, grid);
  Fill(owning, 3.0);

  const auto points =
      static_cast<std::size_t>(radial.NumberOfRadii() * grid.FieldSize());
  auto a = std::vector<Complex>(points);
  auto b = std::vector<Complex>(points);
  auto c = std::vector<Complex>(points);
  auto arrays = std::array{&a, &b, &c};
  auto slot = std::size_t{0};
  LayeredTensorDetails::ForEachStoredSlot(
      owning, [&]<auto indices, ComponentConstraint>(auto& stack) {
        std::ranges::copy(stack.Data(), arrays[slot++]->begin());
      });
  auto view = LayeredTensorFieldView<1, NoSymmetry<1>, ComplexTensor, Grid>(
      radial, grid, std::span(a), std::span(b), std::span(c));

  // Expanding the view, and expanding into a view of the caller's arrays.
  const auto expected = Expand(owning, lMax);
  // Each component's coefficient count is its own upper index's.
  const auto coefficients = [&](std::size_t slot) {
    const auto n = Vector::ComponentLayout.upperIndexOfSlot[slot];
    return static_cast<std::size_t>(radial.NumberOfRadii() *
                                    grid.CoefficientSize(lMax, n));
  };
  auto ca = std::vector<Complex>(coefficients(0));
  auto cb = std::vector<Complex>(coefficients(1));
  auto cc = std::vector<Complex>(coefficients(2));
  auto out = LayeredTensorExpansionView<1, NoSymmetry<1>, ComplexTensor, Grid>(
      radial, grid, lMax, std::span(ca), std::span(cb), std::span(cc));
  Expand(view, out);
  auto outSlot = std::size_t{0};
  auto outArrays = std::array{&ca, &cb, &cc};
  LayeredTensorDetails::ForEachStoredSlot(
      expected, [&]<auto indices, ComponentConstraint>(const auto& stack) {
        const auto& mine = *outArrays[outSlot++];
        for (auto k = std::size_t{0}; k < stack.Data().size(); ++k) {
          ASSERT_EQ(stack.Data()[k], mine[k]) << "element " << k;
        }
      });

  // Back into the view, and an expression assigned into it.
  // (Against the owning route's Evaluate of the same coefficients: the
  // samples were not band-limited, so this is not the identity.)
  std::ranges::fill(a, Complex{});
  Evaluate(out, view);
  ExpectSameAtEveryRadius(Evaluate(expected),
                          [&](Int i) { return FlatSlice(view, i); });
  auto w = Weight(radial, grid);
  for (auto& x : w.Data()) x = 2.0;
  Assign(view, w * owning);
  ExpectSameAtEveryRadius(Materialise(w * owning),
                          [&](Int i) { return FlatSlice(view, i); });
}
