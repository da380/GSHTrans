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
