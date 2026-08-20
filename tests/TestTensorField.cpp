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

  // The component at multi-index (+1) is the third stored one, since the flat
  // order runs from all -1.
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
