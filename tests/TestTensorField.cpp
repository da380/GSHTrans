#include <gtest/gtest.h>

#include <GSHTrans/All>

#include <array>
#include <complex>
#include <cstddef>
#include <type_traits>
#include <vector>

namespace {

using namespace GSHTrans;
using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

auto TestGrid(Int lMax = 4, Int nMax = 4) {
  return Grid(lMax, nMax, FFTWpp::Estimate);
}

// Whether a component is readable, and whether it is writable, asked as a
// concept rather than as a bare requires-expression: inside a non-template
// context GCC diagnoses "no matching function" eagerly instead of reporting
// the requirement as unsatisfied, so the negative cases have to be asked from
// somewhere that substitutes.
template <typename T, std::ptrdiff_t... Alphas>
concept Readable = requires(const T& t) { t.template Component<Alphas...>(); };

template <typename T, std::ptrdiff_t... Alphas>
concept Assignable = requires(T& t) { t.template Component<Alphas...>(); };

// A value that identifies the component it was written into, so that a later
// read can say which one it actually reached.
Complex Marker(Int component, Int point) {
  return Complex{static_cast<Real>(100 * component + point),
                 static_cast<Real>(-component)};
}

}  // namespace

//--------------------------------------------------------------------------//
//                            Storage and layout                             //
//--------------------------------------------------------------------------//

TEST(TensorField, StoresOneComponentPerOrbit) {
  using NoSym = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  using Sym = TensorField<2, Symmetric<2>, ComplexTensor, Grid>;
  using Elastic = TensorField<4, ElasticSymmetry, ComplexTensor, Grid>;

  static_assert(NoSym::Components == 9);
  static_assert(NoSym::StoredComponents == 9);
  static_assert(Sym::StoredComponents == 6);
  static_assert(Elastic::StoredComponents == 21);

  auto grid = TestGrid();
  const auto fieldSize = static_cast<Int>(grid.FieldSize());

  auto t = Sym(grid);
  EXPECT_EQ(t.Size(), Sym::StoredComponents * fieldSize);
  EXPECT_EQ(t.FieldSize(), fieldSize);

  // One buffer, not six allocations: the components are contiguous slices of
  // it, which is what lets a batch of them be transformed together.
  EXPECT_EQ(t.Data().size(), static_cast<std::size_t>(t.Size()));
}

TEST(TensorField, ComponentsCarryTheUpperIndexOfTheirMultiIndex) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);

  // Theory note table 1, read off the returned nodes' types.
  static_assert(decltype(t.Component<-1, -1>())::UpperIndex == -2);
  static_assert(decltype(t.Component<-1, 0>())::UpperIndex == -1);
  static_assert(decltype(t.Component<0, -1>())::UpperIndex == -1);
  static_assert(decltype(t.Component<-1, 1>())::UpperIndex == 0);
  static_assert(decltype(t.Component<0, 0>())::UpperIndex == 0);
  static_assert(decltype(t.Component<1, -1>())::UpperIndex == 0);
  static_assert(decltype(t.Component<0, 1>())::UpperIndex == 1);
  static_assert(decltype(t.Component<1, 0>())::UpperIndex == 1);
  static_assert(decltype(t.Component<1, 1>())::UpperIndex == 2);

  // Three distinct components share N = 0, which is the fact the multi-index
  // exists to keep. Their *types* are identical -- a node is labelled by its
  // upper index and nothing else, which is exactly why a collection labelled
  // only by N does not determine a tensor -- so what distinguishes them is the
  // storage they name.
  static_assert(std::same_as<decltype(t.Component<-1, 1>()),
                             decltype(t.Component<0, 0>())>);

  t.Component<-1, 1>()[0, 0] = Complex{1.0, 0.0};
  t.Component<0, 0>()[0, 0] = Complex{2.0, 0.0};
  t.Component<1, -1>()[0, 0] = Complex{3.0, 0.0};
  EXPECT_EQ((t.Component<-1, 1>()[0, 0]), (Complex{1.0, 0.0}));
  EXPECT_EQ((t.Component<0, 0>()[0, 0]), (Complex{2.0, 0.0}));
  EXPECT_EQ((t.Component<1, -1>()[0, 0]), (Complex{3.0, 0.0}));
}

TEST(TensorField, ComponentsAreSpinWeightedNodes) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  const auto& constT = t;

  static_assert(SpinWeighted<decltype(t.Component<1, 0>())>);
  static_assert(SpinWeighted<decltype(constT.Component<1, 0>())>);

  // And they participate in the phase-1 algebra like any other node. The
  // product of two components lands at the sum of their upper indices, and
  // conj reverses -- so this is integrable, which a component pair at
  // unequal upper index would not be.
  auto pairing = Integrate(conj(constT.Component<1, 1>()) *
                           constT.Component<1, 1>());
  static_assert(std::same_as<decltype(pairing), Complex>);
  SUCCEED();
}

//--------------------------------------------------------------------------//
//                          Reading and writing                              //
//--------------------------------------------------------------------------//

TEST(TensorField, WritingAComponentIsVisibleThroughTheBuffer) {
  using T = TensorField<1, NoSymmetry<1>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  const auto fieldSize = t.FieldSize();

  auto u = t.Component<1>();
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      u[iTheta, iPhi] = Marker(1, iTheta * grid.NumberOfLongitudes() + iPhi);
    }
  }

  // At rank 1 the multi-index and the upper index coincide, so ordering the
  // buffer by upper index puts (+1) last of the three.
  const auto data = t.Data();
  for (auto i = Int{0}; i < fieldSize; i++) {
    EXPECT_EQ(data[2 * fieldSize + i], Marker(1, i)) << "sample " << i;
  }

  // The others are untouched.
  for (auto i = Int{0}; i < fieldSize; i++) {
    EXPECT_EQ(data[i], Complex{});
  }
}

TEST(TensorField, SymmetricComponentsShareStorage) {
  using T = TensorField<2, Symmetric<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);

  t.Component<0, 1>()[0, 0] = Complex{3.0, -1.0};

  // The transposed component is the same field, not a copy of it.
  EXPECT_EQ((t.Component<1, 0>()[0, 0]), (Complex{3.0, -1.0}));
  EXPECT_EQ((t.Component<0, 1>()[0, 0]), (Complex{3.0, -1.0}));

  // Written the other way round, it is still one place.
  t.Component<1, 0>()[0, 0] = Complex{-2.0, 5.0};
  EXPECT_EQ((t.Component<0, 1>()[0, 0]), (Complex{-2.0, 5.0}));
}

TEST(TensorField, AntisymmetricComponentsAreTheNegativeOfTheirRepresentative) {
  using T = TensorField<2, Antisymmetric<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);

  t.Component<-1, 0>()[0, 0] = Complex{2.0, -3.0};

  // The transposed component is a phase-1 expression, not a view, and it
  // evaluates to the negative.
  const auto& constT = t;
  EXPECT_EQ((constT.Component<0, -1>()[0, 0]), (Complex{-2.0, 3.0}));
  static_assert(!std::same_as<decltype(constT.Component<-1, 0>()),
                              decltype(constT.Component<0, -1>())>);

  // It is still a node, so it still composes.
  static_assert(SpinWeighted<decltype(constT.Component<0, -1>())>);
  static_assert(decltype(constT.Component<0, -1>())::UpperIndex == -1);
}

//--------------------------------------------------------------------------//
//                        What cannot be expressed                           //
//--------------------------------------------------------------------------//

TEST(TensorField, VanishingComponentsAreNotRepresentable) {
  using T = TensorField<2, Antisymmetric<2>, ComplexTensor, Grid>;

  // The diagonal of an antisymmetric tensor is identically zero, and asking
  // for it is a compile error rather than a silently zero field. Exposed as a
  // trait so a compile-time traversal can skip them.
  static_assert(T::Vanishes<0, 0>);
  static_assert(T::Vanishes<1, 1>);
  static_assert(!T::Vanishes<0, 1>);

  using S = TensorField<2, Symmetric<2>, ComplexTensor, Grid>;
  static_assert(!S::Vanishes<0, 0>);
  SUCCEED();
}

TEST(TensorField, SignReversedComponentsAreNotWritable) {
  using T = TensorField<2, Antisymmetric<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  const auto& constT = t;

  // The representative is writable; the component derived from it by an
  // antisymmetric permutation is not, because a view cannot negate on the way
  // in. Reading it is fine either way.
  static_assert(T::Writable<-1, 0>);
  static_assert(!T::Writable<0, -1>);
  static_assert(T::Represents<0, -1>);

  // Asking a non-const tensor for a sign-reversed component is not an error:
  // overload resolution falls through to the const accessor, which reads it.
  // What is unavailable is writing, and the returned type says so -- an
  // expression rather than a view.
  static_assert((Readable<T, -1, 0>));
  static_assert((Readable<T, 0, -1>));
  static_assert(std::same_as<decltype(t.Component<0, -1>()),
                             decltype(constT.Component<0, -1>())>);
  static_assert(!std::same_as<decltype(t.Component<-1, 0>()),
                              decltype(constT.Component<-1, 0>())>);

  // A vanishing component is not readable at all, and neither is a
  // multi-index of the wrong length -- the latter without a hard error inside
  // std::array, which is what the pack-size guard on the constraints is for.
  static_assert(!(Readable<T, 0, 0>));
  static_assert(!(Readable<T, 0>));
  static_assert(!(Readable<T, 0, 0, 0>));
  static_assert(!(Assignable<T, 0, 0>));

  // Every component of a symmetric tensor is writable, which is the common
  // case and the reason this restriction costs little.
  using S = TensorField<2, Symmetric<2>, ComplexTensor, Grid>;
  static_assert(S::Writable<0, 1>);
  static_assert(S::Writable<1, 0>);
  static_assert(S::Writable<0, 0>);
  SUCCEED();
}

//--------------------------------------------------------------------------//
//                     Grouping components for a transform                   //
//--------------------------------------------------------------------------//

// A batch shares grid, degree and upper index, so this is the grouping the
// transform can consume, and its members sit FieldSize apart in the buffer.
TEST(TensorField, StoredComponentsGroupByUpperIndex) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;

  constexpr auto counts = std::array<Int, 5>{
      T::StoredAtUpperIndex(-2).second, T::StoredAtUpperIndex(-1).second,
      T::StoredAtUpperIndex(0).second, T::StoredAtUpperIndex(1).second,
      T::StoredAtUpperIndex(2).second};

  // The trinomial coefficients again, now counted over stored components.
  static_assert(counts[0] == 1);
  static_assert(counts[1] == 2);
  static_assert(counts[2] == 3);
  static_assert(counts[3] == 2);
  static_assert(counts[4] == 1);

  // Symmetry removes some of them, and the total always matches.
  using S = TensorField<2, Symmetric<2>, ComplexTensor, Grid>;
  constexpr auto total = S::StoredAtUpperIndex(-2).second +
                         S::StoredAtUpperIndex(-1).second +
                         S::StoredAtUpperIndex(0).second +
                         S::StoredAtUpperIndex(1).second +
                         S::StoredAtUpperIndex(2).second;
  static_assert(total == S::StoredComponents);
  SUCCEED();
}

//--------------------------------------------------------------------------//
//                            The named ranks                                //
//--------------------------------------------------------------------------//

TEST(TensorField, TheNamedRanksAreWhatTheySay) {
  static_assert(std::same_as<VectorField<ComplexTensor, Grid>,
                             TensorField<1, NoSymmetry<1>, ComplexTensor,
                                         Grid>>);
  static_assert(ElasticTensorField<ComplexTensor, Grid>::StoredComponents ==
                21);

  // A vector is the one rank where the multi-index and the upper index
  // coincide, which is the coincidence that makes rank 2 surprising.
  auto grid = TestGrid();
  auto v = VectorField<ComplexTensor, Grid>(grid);
  static_assert(decltype(v.Component<-1>())::UpperIndex == -1);
  static_assert(decltype(v.Component<0>())::UpperIndex == 0);
  static_assert(decltype(v.Component<1>())::UpperIndex == 1);
  SUCCEED();
}

//--------------------------------------------------------------------------//
//                             The transform                                 //
//--------------------------------------------------------------------------//

// Components sharing an upper index have to be contiguous in the buffer, or
// no (count, stride, dist) descriptor covers them and the batching step F
// exists for is unreachable. This is the property the layout is chosen for.
TEST(TensorField, ComponentsSharingAnUpperIndexAreContiguous) {
  constexpr auto contiguousByUpperIndex = []<typename T>() {
    auto seen = std::ptrdiff_t{0};
    for (auto n = -T::Rank; n <= T::Rank; n++) {
      const auto [first, count] = T::StoredAtUpperIndex(n);
      if (count == 0) continue;
      if (first != seen) return false;
      for (auto slot = first; slot < first + count; slot++) {
        if (T::ComponentLayout.upperIndexOfSlot[slot] != n) return false;
      }
      seen += count;
    }
    return seen == T::StoredComponents;
  };

  static_assert(contiguousByUpperIndex
                    .template operator()<
                        TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>>());
  static_assert(contiguousByUpperIndex
                    .template operator()<
                        TensorField<2, Symmetric<2>, ComplexTensor, Grid>>());
  static_assert(contiguousByUpperIndex
                    .template operator()<
                        TensorField<2, Antisymmetric<2>, ComplexTensor,
                                    Grid>>());
  static_assert(contiguousByUpperIndex
                    .template operator()<
                        TensorField<4, ElasticSymmetry, ComplexTensor, Grid>>());
  SUCCEED();
}

// The round trip, which is the whole bridge in one assertion: fill every
// stored component, transform the tensor, transform it back, and get the same
// fields.
TEST(TensorField, RoundTripsThroughTheSpectralDomain) {
  using T = TensorField<2, Symmetric<2>, ComplexTensor, Grid>;
  constexpr auto lMax = Int{6};

  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto t = T(grid);

  // Band-limited data, so that the round trip is exact up to rounding: a
  // transform back and forth reproduces a field only if the field is in the
  // span of the harmonics the grid resolves.
  auto coefficients = std::vector<Complex>(t.CoefficientSize(lMax));
  for (auto i = std::size_t{0}; i < coefficients.size(); i++) {
    coefficients[i] = Complex{std::cos(0.37 * i), std::sin(0.21 * i)};
  }
  t.InverseTransformation(lMax, coefficients);

  auto back = std::vector<Complex>(coefficients.size());
  t.ForwardTransformation(lMax, back);

  for (auto i = std::size_t{0}; i < coefficients.size(); i++) {
    EXPECT_NEAR(back[i].real(), coefficients[i].real(), 1.0e-11) << "at " << i;
    EXPECT_NEAR(back[i].imag(), coefficients[i].imag(), 1.0e-11) << "at " << i;
  }
}

// The batched call has to give what the components would give one at a time,
// exactly: batching widens the inner loop without reordering any sum.
TEST(TensorField, BatchedComponentsMatchComponentByComponentTransforms) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  constexpr auto lMax = Int{5};

  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto t = T(grid);
  const auto fieldSize = t.FieldSize();

  for (auto i = Int{0}; i < t.Size(); i++) {
    t.Data()[i] = Complex{std::cos(0.11 * i), std::sin(0.29 * i)};
  }

  auto batched = std::vector<Complex>(t.CoefficientSize(lMax));
  t.ForwardTransformation(lMax, batched);

  // The same components, transformed one at a time through the unbatched
  // entry point, in the buffer's order.
  auto offset = std::size_t{0};
  for (auto slot = Int{0}; slot < T::StoredComponents; slot++) {
    const auto n = T::ComponentLayout.upperIndexOfSlot[slot];
    const auto coefficientSize =
        static_cast<Int>(grid.CoefficientSize(lMax, n));
    auto one = std::vector<Complex>(fieldSize);
    std::copy_n(t.Data().begin() + slot * fieldSize, fieldSize, one.begin());
    auto expected = std::vector<Complex>(coefficientSize);
    grid.ForwardTransformation(lMax, n, one, expected);
    for (auto j = Int{0}; j < coefficientSize; j++) {
      EXPECT_EQ(batched[offset + j], expected[j])
          << "slot " << slot << ", coefficient " << j;
    }
    offset += static_cast<std::size_t>(coefficientSize);
  }
  EXPECT_EQ(offset, batched.size());
}

TEST(TensorField, CoefficientSizeAccountsForEachUpperIndex) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  constexpr auto lMax = Int{4};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto t = T(grid);

  // One block per stored component, each sized by its own upper index -- the
  // blocks are not all the same length, which is why this is computed.
  auto expected = Int{0};
  for (auto n = Int{-2}; n <= 2; n++) {
    expected += ComponentsAtUpperIndex<2>(n) *
                static_cast<Int>(grid.CoefficientSize(lMax, n));
  }
  EXPECT_EQ(t.CoefficientSize(lMax), expected);

  auto tooSmall = std::vector<Complex>(t.CoefficientSize(lMax) - 1);
  EXPECT_THROW(t.ForwardTransformation(lMax, tooSmall), std::invalid_argument);
}

TEST(TensorField, RejectsAGridThatCannotCarryItsUpperIndices) {
  // A rank-2 tensor has components at N = +-2, so a grid built for one is not
  // enough. Caught at construction rather than at the first component that
  // asks.
  auto narrow = Grid(6, 1, FFTWpp::Estimate);
  EXPECT_THROW((TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>(narrow)),
               std::invalid_argument);

  auto wide = Grid(6, 2, FFTWpp::Estimate);
  EXPECT_NO_THROW((TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>(wide)));
}

//--------------------------------------------------------------------------//
//                             The layout policy                             //
//--------------------------------------------------------------------------//
//
// The two layouts must be interchangeable in every respect except where the
// numbers sit. Nothing above this line asked which one it was using, and the
// transform did not need a repack for either.

namespace {

using PointMajorTensor =
    TensorField<2, Symmetric<2>, ComplexTensor, Grid, PointMajor>;
using ComponentMajorTensor =
    TensorField<2, Symmetric<2>, ComplexTensor, Grid, ComponentMajor>;

}  // namespace

TEST(TensorField, PointMajorInterleavesWhatComponentMajorSeparates) {
  auto grid = TestGrid();
  auto point = PointMajorTensor(grid);
  auto component = ComponentMajorTensor(grid);

  EXPECT_EQ(point.Size(), component.Size());
  static_assert(PointMajorTensor::StoredComponents ==
                ComponentMajorTensor::StoredComponents);

  // Writing through a component view puts the samples in different places,
  // and the view is what knows where.
  point.Component<0, 1>()[0, 0] = Complex{7.0, -2.0};
  component.Component<0, 1>()[0, 0] = Complex{7.0, -2.0};

  const auto stored = PointMajorTensor::StoredComponents;
  const auto slot = PointMajorTensor::SlotOfFlat(
      MultiIndex<2>(std::array<Int, 2>{0, 1}).Flat());

  EXPECT_EQ(point.Data()[slot], (Complex{7.0, -2.0}));
  EXPECT_EQ(component.Data()[slot * component.FieldSize()],
            (Complex{7.0, -2.0}));

  // The second sample of the same component is one field away in one layout
  // and one component away in the other.
  point.Component<0, 1>()[0, 1] = Complex{1.0, 1.0};
  EXPECT_EQ(point.Data()[slot + stored], (Complex{1.0, 1.0}));
}

TEST(TensorField, ComponentViewsAreStridedInPointMajor) {
  auto grid = TestGrid();
  auto t = PointMajorTensor(grid);

  auto u = t.Component<0, 1>();
  EXPECT_EQ(u.Stride(), PointMajorTensor::StoredComponents);
  EXPECT_EQ(u.Size(), t.FieldSize());

  auto v = ComponentMajorTensor(grid).Component<0, 1>();
  EXPECT_EQ(v.Stride(), 1);

  // A strided view is a node like any other: it evaluates, and it composes.
  for (auto i = Int{0}; i < t.FieldSize(); i++) {
    u[i / grid.NumberOfLongitudes(), i % grid.NumberOfLongitudes()] =
        Complex{static_cast<Real>(i), 0.0};
  }
  auto target = std::vector<Complex>(t.FieldSize());
  const auto& constU = u;
  constU.EvaluateInto(std::span(target));
  for (auto i = Int{0}; i < t.FieldSize(); i++) {
    EXPECT_EQ(target[i], (Complex{static_cast<Real>(i), 0.0})) << "at " << i;
  }

  static_assert(SpinWeighted<decltype(u)>);
}

// The point of [C9]: a point-major tensor is transformable in place, with the
// batch descriptor doing the work a repack would otherwise have to.
TEST(TensorField, BothLayoutsTransformToTheSameCoefficients) {
  constexpr auto lMax = Int{5};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  auto point = TensorField<2, Symmetric<2>, ComplexTensor, Grid, PointMajor>(
      grid);
  auto component =
      TensorField<2, Symmetric<2>, ComplexTensor, Grid, ComponentMajor>(grid);

  // The same field in each, written through the component views so that
  // neither test knows the layout.
  const auto Fill = [&](auto& t) {
    const auto nPhi = grid.NumberOfLongitudes();
    const auto write = [&](auto&& u, Int tag) {
      for (auto iTheta : grid.CoLatitudeIndices()) {
        for (auto iPhi : grid.LongitudeIndices()) {
          u[iTheta, iPhi] = Complex{std::cos(0.3 * (iTheta * nPhi + iPhi) + tag),
                                    std::sin(0.7 * iPhi - tag)};
        }
      }
    };
    write(t.template Component<-1, -1>(), 0);
    write(t.template Component<-1, 0>(), 1);
    write(t.template Component<-1, 1>(), 2);
    write(t.template Component<0, 0>(), 3);
    write(t.template Component<0, 1>(), 4);
    write(t.template Component<1, 1>(), 5);
  };
  Fill(point);
  Fill(component);

  auto a = std::vector<Complex>(component.CoefficientSize(lMax));
  auto b = std::vector<Complex>(point.CoefficientSize(lMax));
  ASSERT_EQ(a.size(), b.size());

  component.ForwardTransformation(lMax, a);
  point.ForwardTransformation(lMax, b);

  // Exactly equal: the strided read happens at the pack seam, which copies
  // into the plan's own buffers either way, so nothing downstream of it can
  // tell the difference.
  for (auto i = std::size_t{0}; i < a.size(); i++) {
    EXPECT_EQ(a[i], b[i]) << "coefficient " << i;
  }

  // And back again, into the interleaved layout.
  point.InverseTransformation(lMax, b);
  component.InverseTransformation(lMax, a);
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      EXPECT_EQ((point.Component<0, 1>()[iTheta, iPhi]),
                (component.Component<0, 1>()[iTheta, iPhi]));
    }
  }
}

//--------------------------------------------------------------------------//
//                            Tangential tensors                             //
//--------------------------------------------------------------------------//

// A tangential tensor has no radial slot: its indices are drawn from {-1, +1}
// and it has 2^Rank components rather than 3^Rank. It is another object in
// another bundle, not a general tensor that happens to vanish in some
// directions (field-algebra-plan.md section 18.2 [D8]).
//
// Nothing below is a special case inside the library. Every count comes out of
// the same orbit walk over whichever multi-indices exist, which is the check
// that the alphabet is a generalisation rather than a second implementation.
namespace {

using Tangential2 = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid,
                                ComponentMajor, TangentialSlots>;
using RealTangential2 = TensorField<2, NoSymmetry<2>, RealTensor, Grid,
                                    ComponentMajor, TangentialSlots>;
using RealTangentialSym2 = TensorField<2, Symmetric<2>, RealTensor, Grid,
                                       ComponentMajor, TangentialSlots>;

}  // namespace

TEST(TensorField, ATangentialTensorHasTwoLettersPerSlot) {
  static_assert(std::same_as<Tangential2::SlotSet, TangentialSlots>);
  static_assert(Tangential2::Components == 4);
  static_assert(Tangential2::StoredComponents == 4);
  static_assert(Tangential2::RealComponents == 0);
  static_assert(Tangential2::RealsPerPoint == 8);

  // The upper index is still the signed sum, so it has the parity of the rank
  // and there is no component at an odd one.
  static_assert(Tangential2::UpperIndexOf<-1, -1> == -2);
  static_assert(Tangential2::UpperIndexOf<-1, 1> == 0);
  static_assert(Tangential2::UpperIndexOf<1, -1> == 0);
  static_assert(Tangential2::UpperIndexOf<1, 1> == 2);

  // And the default is unchanged, which is what makes the parameter additive.
  using General2 = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  static_assert(std::same_as<General2::SlotSet, AllSlots>);
  static_assert(General2::Components == 9);
  SUCCEED();
}

// The negative case the letter check exists for. Without it, asking a
// tangential tensor for a radial component is a hard error inside the
// multi-index constructor rather than an unsatisfied constraint -- and every
// assertion below would be vacuous, because a requires-expression cannot see
// a throw in a constant expression (see IsSlotLetter in MultiIndex.h).
TEST(TensorField, ARadialComponentOfATangentialTensorIsNotAComponent) {
  static_assert(Tangential2::Represents<-1, 1>);
  static_assert(!Tangential2::Represents<0, 1>);
  static_assert(!Tangential2::Represents<1, 0>);
  static_assert(!Tangential2::Represents<0, 0>);
  static_assert(!Tangential2::Writable<0, 1>);

  // Not vanishing: a vanishing component is one the tensor has and that is
  // identically zero, which is a different thing from one it does not have.
  static_assert(!Tangential2::Vanishes<0, 1>);

  static_assert(Readable<Tangential2, -1, 1>);
  static_assert(!Readable<Tangential2, 0, 1>);
  static_assert(!Assignable<Tangential2, 0, 1>);

  // The wrong number of indices is rejected as before, and for the same
  // reason it always was.
  static_assert(!Readable<Tangential2, -1>);
  static_assert(!Readable<Tangential2, -1, 1, 1>);
  SUCCEED();
}

// The payoff of the alphabet, and the thing phase 4 could not do: negation
// has no fixed point when there is no zero letter, so every orbit has size
// two, nothing is pinned, and the second buffer is empty
// (field-algebra-plan.md section 18.2 [D6]).
TEST(TensorField, ARealTangentialTensorHasNoPinnedComponents) {
  static_assert(RealTangential2::StoredComponents == 2);
  static_assert(RealTangential2::ComplexComponents == 2);
  static_assert(RealTangential2::RealComponents == 0);
  static_assert(RealTangential2::RealsPerPoint == 4);

  auto grid = TestGrid();
  auto t = RealTangential2(grid);
  EXPECT_EQ(t.RealSize(), 0);
  EXPECT_EQ(t.Size(), 2 * t.FieldSize());

  // The stored components are the two representatives; their partners are
  // derived by T^{-alpha} = (-1)^N conj(T^{alpha}), which at N = -2 and at
  // N = 0 is plain conjugation.
  static_assert(RealTangential2::Writable<-1, -1>);
  static_assert(RealTangential2::Writable<-1, 1>);
  static_assert(!RealTangential2::Writable<1, 1>);
  static_assert(!RealTangential2::Writable<1, -1>);

  t.Component<-1, -1>()[1, 2] = Complex{3.0, -4.0};
  t.Component<-1, 1>()[1, 2] = Complex{5.0, 6.0};

  const auto& tensor = t;
  EXPECT_EQ((tensor.Component<1, 1>()[1, 2]), (Complex{3.0, 4.0}));
  EXPECT_EQ((tensor.Component<1, -1>()[1, 2]), (Complex{5.0, -6.0}));
}

// Under a permutation symmetry a self-paired component can still appear, and
// the orbit walk has to find that unaided: negation maps (-+) to (+-) and the
// symmetry maps it back, so that component is pinned real. One complex plus
// one real is three reals a point, which is a real symmetric 2x2 matrix.
TEST(TensorField, ASymmetricRealTangentialTensorIsARealSymmetricTwoByTwo) {
  static_assert(RealTangentialSym2::StoredComponents == 2);
  static_assert(RealTangentialSym2::ComplexComponents == 1);
  static_assert(RealTangentialSym2::RealComponents == 1);
  static_assert(RealTangentialSym2::RealsPerPoint == 3);

  auto grid = TestGrid();
  auto t = RealTangentialSym2(grid);
  EXPECT_EQ(t.RealSize(), t.FieldSize());

  // The pinned component is a real-valued field, and it is at upper index
  // zero -- which it has to be, since permutation preserves the slot sum and
  // negation reverses it.
  auto pinned = t.Component<-1, 1>();
  static_assert(std::same_as<decltype(pinned)::Value, RealValued>);
  static_assert(decltype(pinned)::UpperIndex == 0);

  pinned[1, 2] = 7.0;
  const auto& tensor = t;
  EXPECT_EQ((tensor.Component<1, -1>()[1, 2]), 7.0);
}

// The transform with an empty real buffer, which nothing has exercised before:
// phase 4 introduced that buffer and no tensor until now has had none of it.
TEST(TensorField, ARealTangentialTensorRoundTripsWithNoRealBuffer) {
  constexpr auto lMax = Int{6};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto t = RealTangential2(grid);

  const auto size = t.CoefficientSize(lMax);
  // Two complex blocks and no real one: the stored components sit at upper
  // index -2 and 0, and a block's length depends on its upper index.
  EXPECT_EQ(size, static_cast<Int>(grid.CoefficientSize(lMax, -2)) +
                      static_cast<Int>(grid.CoefficientSize(lMax, 0)));

  auto coefficients = std::vector<Complex>(size);
  for (auto i = std::size_t{0}; i < coefficients.size(); i++) {
    coefficients[i] = Complex{std::cos(0.31 * i), std::sin(0.17 * i)};
  }
  t.InverseTransformation(lMax, coefficients);

  auto back = std::vector<Complex>(coefficients.size());
  t.ForwardTransformation(lMax, back);

  for (auto i = std::size_t{0}; i < coefficients.size(); i++) {
    EXPECT_NEAR(back[i].real(), coefficients[i].real(), 1.0e-11) << "at " << i;
    EXPECT_NEAR(back[i].imag(), coefficients[i].imag(), 1.0e-11) << "at " << i;
  }
}
