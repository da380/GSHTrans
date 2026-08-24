#include <gtest/gtest.h>

#include <GSHTrans/All>
#include <array>
#include <complex>
#include <cstddef>
#include <type_traits>

namespace {

using namespace GSHTrans;
using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

auto TestGrid() { return Grid(4, 4, FFTWpp::Estimate); }

// Distinct data per component, so that a permuted read can be told from an
// unpermuted one.
template <typename T, Int... Alphas>
void Fill(T& t, Real tag) {
  auto u = t.template Component<Alphas...>();
  const auto& grid = t.Grid();
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      u[iTheta, iPhi] = Complex{tag + iTheta, static_cast<Real>(iPhi)};
    }
  }
}

}  // namespace

//--------------------------------------------------------------------------//
//                            Permuting the slots                            //
//--------------------------------------------------------------------------//

TEST(TensorAlgebra, TensorFieldIsATensorExpression) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  static_assert(TensorExpr<T>);
  static_assert(IsTerminal<T>);
  SUCCEED();
}

TEST(TensorAlgebra, TransposeRelabelsWithoutCopying) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);

  Fill<T, 0, 1>(t, 1.0);
  Fill<T, 1, 0>(t, 2.0);

  const auto& tensor = t;
  auto transposed = Transpose(tensor);
  static_assert(TensorExpr<decltype(transposed)>);
  static_assert(decltype(transposed)::Rank == 2);

  // Slot 0 of the transpose is slot 1 of the operand.
  EXPECT_EQ((transposed.Component<0, 1>()[2, 3]),
            (tensor.Component<1, 0>()[2, 3]));
  EXPECT_EQ((transposed.Component<1, 0>()[2, 3]),
            (tensor.Component<0, 1>()[2, 3]));
  EXPECT_NE((transposed.Component<0, 1>()[2, 3]),
            (tensor.Component<0, 1>()[2, 3]));

  // The upper index follows the multi-index, not the slot order, so a
  // transposed component carries what its own indices say.
  static_assert(decltype(transposed.Component<0, 1>())::UpperIndex == 1);
  static_assert(decltype(transposed.Component<-1, 1>())::UpperIndex == 0);

  // Twice is the identity.
  auto back = Transpose(Transpose(tensor));
  EXPECT_EQ((back.Component<0, 1>()[2, 3]), (tensor.Component<0, 1>()[2, 3]));
}

TEST(TensorAlgebra, TransposingASymmetricTensorChangesNothing) {
  using T = TensorField<2, Symmetric<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  Fill<T, 0, 1>(t, 3.0);

  const auto& tensor = t;
  auto transposed = Transpose(tensor);
  EXPECT_EQ((transposed.Component<0, 1>()[1, 1]),
            (tensor.Component<0, 1>()[1, 1]));
}

TEST(TensorAlgebra, TransposingAnAntisymmetricTensorNegatesIt) {
  using T = TensorField<2, Antisymmetric<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  Fill<T, -1, 0>(t, 4.0);

  const auto& tensor = t;
  auto transposed = Transpose(tensor);
  EXPECT_EQ((transposed.Component<-1, 0>()[1, 2]),
            -(tensor.Component<-1, 0>()[1, 2]));

  // A component the operand cannot represent, the transpose cannot either.
  static_assert(!(decltype(transposed)::Represents<0, 0>));
  static_assert((decltype(transposed)::Represents<-1, 0>));
}

TEST(TensorAlgebra, HigherRankSlotsPermuteAsAsked) {
  using T = TensorField<4, NoSymmetry<4>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  Fill<T, 1, 0, -1, 0>(t, 5.0);

  const auto& tensor = t;

  // Send slot 0 to where slot 2 was, and so on: a cyclic shift.
  auto rotated = Permute<std::array<Int, 4>{2, 3, 0, 1}>(tensor);
  EXPECT_EQ((rotated.Component<-1, 0, 1, 0>()[0, 1]),
            (tensor.Component<1, 0, -1, 0>()[0, 1]));

  // The elastic symmetry contains exactly that permutation, so applying it
  // there is the identity.
  using E = TensorField<4, ElasticSymmetry, ComplexTensor, Grid>;
  auto e = E(grid);
  Fill<E, 1, 0, -1, 0>(e, 6.0);
  const auto& elastic = e;
  auto same = Permute<std::array<Int, 4>{2, 3, 0, 1}>(elastic);
  EXPECT_EQ((same.Component<1, 0, -1, 0>()[0, 1]),
            (elastic.Component<1, 0, -1, 0>()[0, 1]));
}

TEST(TensorAlgebra, ExpressionsComposeWithThePhaseOneAlgebra) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  Fill<T, 0, 1>(t, 1.0);
  Fill<T, 1, 0>(t, 2.0);

  const auto& tensor = t;
  auto transposed = Transpose(tensor);

  // The antisymmetric part of one component pair, formed as a phase-1
  // expression over two tensor expressions. Both terms are at upper index 1,
  // so the subtraction is admissible; had they not been, this would not
  // compile.
  auto skew = 0.5 * (tensor.Component<0, 1>() - transposed.Component<0, 1>());
  static_assert(decltype(skew)::UpperIndex == 1);

  const auto expected =
      0.5 * (tensor.Component<0, 1>()[2, 2] - tensor.Component<1, 0>()[2, 2]);
  EXPECT_EQ((skew[2, 2]), expected);
}

//--------------------------------------------------------------------------//
//                             The tensor product                            //
//--------------------------------------------------------------------------//

TEST(TensorAlgebra, TheProductConcatenatesMultiIndices) {
  using V = TensorField<1, NoSymmetry<1>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto a = V(grid);
  auto b = V(grid);
  Fill<V, 1>(a, 1.0);
  Fill<V, -1>(b, 2.0);

  const auto& u = a;
  const auto& v = b;
  auto product = TensorProduct(u, v);
  static_assert(TensorExpr<decltype(product)>);
  static_assert(decltype(product)::Rank == 2);

  // Upper indices add, which for a product of components is eq:N applied to
  // the concatenated multi-index.
  static_assert(decltype(product.Component<1, -1>())::UpperIndex == 0);
  static_assert(decltype(product.Component<1, 1>())::UpperIndex == 2);

  EXPECT_EQ((product.Component<1, -1>()[2, 1]),
            (u.Component<1>()[2, 1] * v.Component<-1>()[2, 1]));
}

TEST(TensorAlgebra, ProductsOfDifferentRanksSplitTheIndexCorrectly) {
  using V = TensorField<1, NoSymmetry<1>, ComplexTensor, Grid>;
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto v = V(grid);
  auto t = T(grid);
  Fill<V, 0>(v, 3.0);
  Fill<T, 1, -1>(t, 4.0);

  const auto& vector = v;
  const auto& tensor = t;

  // Rank 1 then rank 2: the first slot is the vector's.
  auto left = TensorProduct(vector, tensor);
  static_assert(decltype(left)::Rank == 3);
  EXPECT_EQ((left.Component<0, 1, -1>()[1, 1]),
            (vector.Component<0>()[1, 1] * tensor.Component<1, -1>()[1, 1]));

  // Rank 2 then rank 1: the last slot is.
  auto right = TensorProduct(tensor, vector);
  EXPECT_EQ((right.Component<1, -1, 0>()[1, 1]),
            (tensor.Component<1, -1>()[1, 1] * vector.Component<0>()[1, 1]));
}

TEST(TensorAlgebra, ProductsRejectOperandsOnDifferentGrids) {
  using V = TensorField<1, NoSymmetry<1>, ComplexTensor, Grid>;
  auto one = TestGrid();
  auto other = TestGrid();
  auto a = V(one);
  auto b = V(other);
  const auto& u = a;
  const auto& v = b;

  // Same parameters, different implementations: identity, not structure.
  EXPECT_THROW(TensorProduct(u, v), std::invalid_argument);
  EXPECT_NO_THROW(TensorProduct(u, u));
}

//--------------------------------------------------------------------------//
//                                Contraction                                //
//--------------------------------------------------------------------------//

TEST(TensorAlgebra, TraceIsTheMetricContraction) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  Fill<T, -1, 1>(t, 1.0);
  Fill<T, 0, 0>(t, 2.0);
  Fill<T, 1, -1>(t, 3.0);

  const auto& tensor = t;
  auto trace = Trace(tensor);

  // g_{ab} = (-1)^a delta_{a+b,0}, so the trace is -T^{-+} + T^{00} - T^{+-}.
  static_assert(decltype(trace)::UpperIndex == 0);
  const auto expected = -tensor.Component<-1, 1>()[2, 2] +
                        tensor.Component<0, 0>()[2, 2] -
                        tensor.Component<1, -1>()[2, 2];
  EXPECT_EQ((trace[2, 2]), expected);

  // Only the N = 0 components contribute, which is what makes the trace a
  // scalar: the contracted pair adds a + (-a) = 0 whatever a is.
  Fill<T, 1, 1>(t, 9.0);
  EXPECT_EQ((trace[2, 2]), expected);
}

TEST(TensorAlgebra, ContractionOfHigherRankLeavesTheSurvivingSlots) {
  using T = TensorField<4, NoSymmetry<4>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  Fill<T, 1, -1, 0, 1>(t, 1.0);
  Fill<T, 1, 0, 0, 1>(t, 2.0);
  Fill<T, 1, 1, 0, 1>(t, 3.0);

  const auto& tensor = t;

  // Contract slots 1 and 2, leaving slots 0 and 3 in order.
  auto contracted = Contract<1, 2>(tensor);
  static_assert(decltype(contracted)::Rank == 2);
  static_assert(decltype(contracted.Component<1, 1>())::UpperIndex == 2);

  // The surviving indices are (1, 1); the contracted pair runs over a.
  const auto expected = -tensor.Component<1, -1, 1, 1>()[1, 2] +
                        tensor.Component<1, 0, 0, 1>()[1, 2] -
                        tensor.Component<1, 1, -1, 1>()[1, 2];
  EXPECT_EQ((contracted.Component<1, 1>()[1, 2]), expected);
}

TEST(TensorAlgebra, ContractingAProductIsADoubleSum) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto a = T(grid);
  auto b = T(grid);
  for (auto i = Int{0}; i < a.Size(); i++) {
    a.Data()[i] = Complex{std::cos(0.3 * i), std::sin(0.7 * i)};
    b.Data()[i] = Complex{std::sin(0.2 * i), std::cos(0.5 * i)};
  }

  const auto& s = a;
  const auto& t = b;

  // S^{ab} T^{cd} contracted on b and c: the matrix product in canonical
  // components, which carries the metric with it.
  auto product = Contract<1, 2>(TensorProduct(s, t));
  static_assert(decltype(product)::Rank == 2);

  const auto expected =
      -s.Component<0, -1>()[1, 1] * t.Component<1, 0>()[1, 1] +
      s.Component<0, 0>()[1, 1] * t.Component<0, 0>()[1, 1] -
      s.Component<0, 1>()[1, 1] * t.Component<-1, 0>()[1, 1];
  const auto got = product.Component<0, 0>()[1, 1];
  EXPECT_NEAR(got.real(), expected.real(), 1.0e-13);
  EXPECT_NEAR(got.imag(), expected.imag(), 1.0e-13);
}

//--------------------------------------------------------------------------//
//                      Symmetrisation and materialisation                   //
//--------------------------------------------------------------------------//

TEST(TensorAlgebra, SymmetrisationProjectsOntoTheSymmetry) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  Fill<T, 0, 1>(t, 1.0);
  Fill<T, 1, 0>(t, 5.0);

  const auto& tensor = t;
  auto sym = Symmetrise<Symmetric<2>>(tensor);
  auto skew = Symmetrise<Antisymmetric<2>>(tensor);

  static_assert(TensorExpr<decltype(sym)>);
  static_assert(decltype(sym)::Rank == 2);

  const auto a = tensor.Component<0, 1>()[2, 2];
  const auto b = tensor.Component<1, 0>()[2, 2];

  EXPECT_NEAR((sym.Component<0, 1>()[2, 2]).real(), (0.5 * (a + b)).real(),
              1.0e-13);
  EXPECT_NEAR((skew.Component<0, 1>()[2, 2]).real(), (0.5 * (a - b)).real(),
              1.0e-13);

  // The symmetric part is symmetric, and the antisymmetric part changes sign.
  EXPECT_NEAR((sym.Component<1, 0>()[2, 2]).real(),
              (sym.Component<0, 1>()[2, 2]).real(), 1.0e-13);
  EXPECT_NEAR((skew.Component<1, 0>()[2, 2]).real(),
              -(skew.Component<0, 1>()[2, 2]).real(), 1.0e-13);

  // And they add back up to the original.
  EXPECT_NEAR(
      (sym.Component<0, 1>()[2, 2] + skew.Component<0, 1>()[2, 2]).real(),
      a.real(), 1.0e-13);
}

TEST(TensorAlgebra, SymmetrisingAnAlreadySymmetricTensorIsTheIdentity) {
  using T = TensorField<2, Symmetric<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  Fill<T, 0, 1>(t, 2.0);

  const auto& tensor = t;
  auto sym = Symmetrise<Symmetric<2>>(tensor);
  EXPECT_NEAR((sym.Component<0, 1>()[1, 1]).real(),
              (tensor.Component<0, 1>()[1, 1]).real(), 1.0e-13);

  // Antisymmetrising it gives zero, which is the check that the signs are the
  // right way round.
  auto skew = Symmetrise<Antisymmetric<2>>(tensor);
  EXPECT_NEAR((skew.Component<0, 1>()[1, 1]).real(), 0.0, 1.0e-13);
  EXPECT_NEAR((skew.Component<0, 1>()[1, 1]).imag(), 0.0, 1.0e-13);
}

TEST(TensorAlgebra, TheGroupIsClosedUnderItsGenerators) {
  // Two generators give the whole symmetric group on three slots, and six
  // elements is what that is.
  static_assert(TensorDetails::GroupElements<3, Symmetric<3>>().second == 6);
  static_assert(TensorDetails::GroupElements<2, Symmetric<2>>().second == 2);
  static_assert(TensorDetails::GroupElements<4, Symmetric<4>>().second == 24);
  static_assert(TensorDetails::GroupElements<2, NoSymmetry<2>>().second == 1);

  // The elastic symmetry is a group of order 8: two independent pair swaps
  // and the exchange of the pairs.
  static_assert(TensorDetails::GroupElements<4, ElasticSymmetry>().second == 8);
  SUCCEED();
}

TEST(TensorAlgebra, MaterialiseEvaluatesIntoAField) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  for (auto i = Int{0}; i < t.Size(); i++) {
    t.Data()[i] = Complex{std::cos(0.4 * i), std::sin(0.9 * i)};
  }

  const auto& tensor = t;
  auto transposed = Materialise(Transpose(tensor));
  static_assert(
      std::same_as<decltype(transposed),
                   TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>>);

  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      EXPECT_EQ((transposed.Component<0, 1>()[iTheta, iPhi]),
                (tensor.Component<1, 0>()[iTheta, iPhi]));
      EXPECT_EQ((transposed.Component<-1, 1>()[iTheta, iPhi]),
                (tensor.Component<1, -1>()[iTheta, iPhi]));
    }
  }
}

TEST(TensorAlgebra, MaterialiseStoresOnlyWhatTheAskedSymmetryKeeps) {
  using T = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);
  Fill<T, 0, 1>(t, 1.0);
  Fill<T, 1, 0>(t, 5.0);

  const auto& tensor = t;

  // Asking for a symmetry is an assertion about the value, honoured by
  // storing only the components it keeps: six rather than nine.
  auto sym = Materialise<Symmetric<2>>(Symmetrise<Symmetric<2>>(tensor));
  static_assert(decltype(sym)::StoredComponents == 6);
  EXPECT_EQ(sym.Size(), 6 * sym.FieldSize());

  const auto expected =
      0.5 * (tensor.Component<0, 1>()[2, 2] + tensor.Component<1, 0>()[2, 2]);
  EXPECT_NEAR((sym.Component<0, 1>()[2, 2]).real(), expected.real(), 1.0e-13);
  EXPECT_NEAR((sym.Component<1, 0>()[2, 2]).real(), expected.real(), 1.0e-13);
}

// The composite the whole layer exists for: an elastic tensor applied to a
// strain, which is a double contraction of a tensor product.
TEST(TensorAlgebra, ElasticTensorAppliedToAStrain) {
  using C = TensorField<4, ElasticSymmetry, ComplexTensor, Grid>;
  using E = TensorField<2, Symmetric<2>, ComplexTensor, Grid>;
  auto grid = TestGrid();

  auto c = C(grid);
  auto e = E(grid);
  for (auto i = Int{0}; i < c.Size(); i++) {
    c.Data()[i] = Complex{std::cos(0.13 * i), 0.0};
  }
  for (auto i = Int{0}; i < e.Size(); i++) {
    e.Data()[i] = Complex{std::sin(0.27 * i), 0.0};
  }

  const auto& elastic = c;
  const auto& strain = e;

  // c^{ijkl} e^{mn} contracted on (k, m) and then on what was (l, n): the
  // stress, of rank 2.
  auto product = TensorProduct(elastic, strain);
  static_assert(decltype(product)::Rank == 6);
  auto once = Contract<2, 4>(product);
  static_assert(decltype(once)::Rank == 4);
  auto stress = Contract<2, 3>(once);
  static_assert(decltype(stress)::Rank == 2);

  // It is still a tensor, its components are still spin-weighted fields at
  // the upper index their multi-index implies, and it can be materialised.
  static_assert(decltype(stress.Component<1, 1>())::UpperIndex == 2);
  auto materialised = Materialise(stress);
  EXPECT_EQ((materialised.Component<1, 1>()[1, 1]),
            (stress.Component<1, 1>()[1, 1]));
}

//--------------------------------------------------------------------------//
//                            Tangential tensors                             //
//--------------------------------------------------------------------------//

namespace {

using Tangential2 = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid,
                                ComponentMajor, TangentialSlots>;
using Tangential1 = TensorField<1, NoSymmetry<1>, ComplexTensor, Grid,
                                ComponentMajor, TangentialSlots>;
using General2 = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;

// Asked as a concept rather than as a bare requires-expression, for the same
// reason the component accessors are: in a non-template context GCC reports
// "no matching function" instead of an unsatisfied requirement.
template <typename L, typename R>
concept Multipliable =
    requires(const L& left, const R& right) { TensorProduct(left, right); };

}  // namespace

TEST(TensorAlgebra, TangentialTensorsAreTensorExpressions) {
  static_assert(TensorExpr<Tangential2>);
  static_assert(IsTerminal<Tangential2>);
  static_assert(std::same_as<Tangential2::SlotSet, TangentialSlots>);

  auto grid = TestGrid();
  auto t = Tangential2(grid);
  Fill<Tangential2, -1, 1>(t, 1.0);

  // Permutation is a relabelling and knows nothing about the alphabet.
  const auto& tensor = t;
  auto transposed = Transpose(tensor);
  static_assert(std::same_as<decltype(transposed)::SlotSet, TangentialSlots>);
  EXPECT_EQ((transposed.Component<1, -1>()[2, 2]),
            (tensor.Component<-1, 1>()[2, 2]));

  // And a slot the tangential bundle does not have is not a component of the
  // permuted tensor either.
  static_assert(!decltype(transposed)::Represents<0, 1>);
}

// The metric contraction runs over the alphabet's letters, so on a tangential
// tensor it is the induced metric of the sphere: the radial term is absent
// because there is no radial slot to contribute one.
TEST(TensorAlgebra, TheTangentialTraceHasNoRadialTerm) {
  auto grid = TestGrid();
  auto t = Tangential2(grid);
  Fill<Tangential2, -1, 1>(t, 1.0);
  Fill<Tangential2, 1, -1>(t, 3.0);
  Fill<Tangential2, -1, -1>(t, 5.0);
  Fill<Tangential2, 1, 1>(t, 7.0);

  const auto& tensor = t;
  auto trace = Trace(tensor);
  static_assert(decltype(trace)::UpperIndex == 0);

  const auto expected =
      -tensor.Component<-1, 1>()[2, 2] - tensor.Component<1, -1>()[2, 2];
  EXPECT_EQ((trace[2, 2]), expected);
}

// Crossing bundles is done by embedding at the call site and never
// implicitly, so a product of operands from different alphabets does not
// compile (field-algebra-plan.md section 18.2 [D10]).
TEST(TensorAlgebra, ProductsRejectOperandsFromDifferentBundles) {
  static_assert(Multipliable<Tangential1, Tangential1>);
  static_assert(Multipliable<General2, General2>);
  static_assert(!Multipliable<Tangential1, General2>);
  static_assert(!Multipliable<General2, Tangential1>);

  auto grid = TestGrid();
  auto u = Tangential1(grid);
  Fill<Tangential1, -1>(u, 1.0);
  Fill<Tangential1, 1>(u, 2.0);

  const auto& vector = u;
  auto product = TensorProduct(vector, vector);
  static_assert(decltype(product)::Rank == 2);
  static_assert(std::same_as<decltype(product)::SlotSet, TangentialSlots>);
  static_assert(!decltype(product)::Represents<0, 1>);
  EXPECT_EQ((product.Component<-1, 1>()[2, 2]),
            (vector.Component<-1>()[2, 2] * vector.Component<1>()[2, 2]));
}

// Materialising cannot move a tensor between bundles: the field it produces
// is over the expression's own alphabet.
TEST(TensorAlgebra, MaterialiseKeepsTheOperandsAlphabet) {
  auto grid = TestGrid();
  auto t = Tangential2(grid);
  Fill<Tangential2, -1, 1>(t, 1.0);
  Fill<Tangential2, 1, -1>(t, 3.0);

  const auto& tensor = t;
  auto symmetrised = Symmetrise<Symmetric<2>>(tensor);
  auto field = Materialise<Symmetric<2>>(symmetrised);

  static_assert(std::same_as<decltype(field)::SlotSet, TangentialSlots>);
  static_assert(decltype(field)::Components == 4);
  static_assert(decltype(field)::StoredComponents == 3);

  const auto expected =
      0.5 * (tensor.Component<-1, 1>()[2, 2] + tensor.Component<1, -1>()[2, 2]);
  EXPECT_EQ((field.Component<-1, 1>()[2, 2]), expected);
  EXPECT_EQ((field.Component<1, -1>()[2, 2]), expected);
}

//--------------------------------------------------------------------------//
//                       The maps between the bundles                        //
//--------------------------------------------------------------------------//

namespace {

template <typename T, Int... Alphas>
concept ReadableComponent =
    requires(const T& t) { t.template Component<Alphas...>(); };

template <typename T>
concept Materialisable = requires(const T& t) { Materialise(t); };

template <typename T>
concept Embeddable = requires(const T& t) { Embed(t); };

template <typename T>
concept Projectable = requires(const T& t) { Tangential(t); };

}  // namespace

// Embed widens the alphabet and nothing else. The components with a radial
// slot are ones the embedded tensor does not have, rather than components it
// has and that are zero -- which is the same statement the orbit table makes
// about a component an antisymmetry annihilates, and it is why nothing is
// allocated.
TEST(TensorAlgebra, EmbedWidensTheAlphabetWithoutStoringAnything) {
  auto grid = TestGrid();
  auto t = Tangential2(grid);
  Fill<Tangential2, -1, 1>(t, 1.0);
  Fill<Tangential2, 1, -1>(t, 3.0);

  const auto& tensor = t;
  auto embedded = Embed(tensor);
  static_assert(TensorExpr<decltype(embedded)>);
  static_assert(std::same_as<decltype(embedded)::SlotSet, AllSlots>);
  static_assert(decltype(embedded)::Rank == 2);

  // The tangential components come through unchanged...
  static_assert(decltype(embedded)::Represents<-1, 1>);
  EXPECT_EQ((embedded.Component<-1, 1>()[2, 2]),
            (tensor.Component<-1, 1>()[2, 2]));

  // ...and the ones with a radial slot are not represented at all.
  static_assert(!decltype(embedded)::Represents<0, 1>);
  static_assert(!decltype(embedded)::Represents<1, 0>);
  static_assert(!decltype(embedded)::Represents<0, 0>);
  static_assert(ReadableComponent<decltype(embedded), -1, 1>);
  static_assert(!ReadableComponent<decltype(embedded), 0, 1>);

  // A general tensor is already where Embed would send it, so it is refused.
  static_assert(Embeddable<Tangential2>);
  static_assert(!Embeddable<General2>);
}

// The projection is the adjoint: it drops exactly those components.
TEST(TensorAlgebra, TangentialDropsTheComponentsWithARadialSlot) {
  auto grid = TestGrid();
  auto t = General2(grid);
  Fill<General2, -1, 1>(t, 1.0);
  Fill<General2, 0, 1>(t, 2.0);
  Fill<General2, 1, -1>(t, 3.0);

  const auto& tensor = t;
  auto projected = Tangential(tensor);
  static_assert(TensorExpr<decltype(projected)>);
  static_assert(std::same_as<decltype(projected)::SlotSet, TangentialSlots>);

  EXPECT_EQ((projected.Component<-1, 1>()[2, 2]),
            (tensor.Component<-1, 1>()[2, 2]));
  static_assert(!decltype(projected)::Represents<0, 1>);

  static_assert(Projectable<General2>);
  static_assert(!Projectable<Tangential2>);
}

// Projecting an embedded tensor is the identity, which is what says the two
// maps are the pair they claim to be.
TEST(TensorAlgebra, ProjectingAnEmbeddedTensorGivesItBack) {
  auto grid = TestGrid();
  auto t = Tangential2(grid);
  Fill<Tangential2, -1, -1>(t, 1.0);
  Fill<Tangential2, -1, 1>(t, 2.0);
  Fill<Tangential2, 1, -1>(t, 3.0);
  Fill<Tangential2, 1, 1>(t, 4.0);

  const auto& tensor = t;
  auto round = Tangential(Embed(tensor));
  static_assert(std::same_as<decltype(round)::SlotSet, TangentialSlots>);

  const auto same = [&]<Int A, Int B>() {
    EXPECT_EQ((round.Component<A, B>()[2, 2]),
              (tensor.Component<A, B>()[2, 2]));
  };
  same.template operator()<-1, -1>();
  same.template operator()<-1, 1>();
  same.template operator()<1, -1>();
  same.template operator()<1, 1>();
}

// What the pair is for: [D10] says a product across bundles is written by
// embedding at the call site, and this is that sentence as code.
TEST(TensorAlgebra, EmbeddingIsHowAProductCrossesBundles) {
  auto grid = TestGrid();
  auto u = Tangential1(grid);
  Fill<Tangential1, -1>(u, 1.0);
  Fill<Tangential1, 1>(u, 2.0);

  auto v = TensorField<1, NoSymmetry<1>, ComplexTensor, Grid>(grid);
  Fill<decltype(v), -1>(v, 3.0);
  Fill<decltype(v), 0>(v, 4.0);
  Fill<decltype(v), 1>(v, 5.0);

  const auto& tangential = u;
  const auto& general = v;
  static_assert(!Multipliable<Tangential1, decltype(v)>);

  auto product = TensorProduct(Embed(tangential), general);
  static_assert(decltype(product)::Rank == 2);
  static_assert(std::same_as<decltype(product)::SlotSet, AllSlots>);

  EXPECT_EQ((product.Component<-1, 0>()[2, 2]),
            (tangential.Component<-1>()[2, 2] * general.Component<0>()[2, 2]));

  // The half of the product that the embedded factor does not reach is not
  // represented, so a traversal skips it rather than evaluating a zero.
  static_assert(!decltype(product)::Represents<0, 0>);
  static_assert(decltype(product)::Represents<1, 0>);
}

// Materialising through the maps lands in the bundle the expression is in,
// and the storage is the smaller one on the way down.
TEST(TensorAlgebra, MaterialisingAProjectionStoresTheTangentialSet) {
  auto grid = TestGrid();
  auto t = General2(grid);
  Fill<General2, -1, 1>(t, 1.0);
  Fill<General2, 1, -1>(t, 3.0);

  const auto& tensor = t;
  auto field = Materialise(Tangential(tensor));
  static_assert(std::same_as<decltype(field)::SlotSet, TangentialSlots>);
  static_assert(decltype(field)::StoredComponents == 4);
  EXPECT_EQ((field.Component<-1, 1>()[2, 2]),
            (tensor.Component<-1, 1>()[2, 2]));
}

// A spectral tensor is not a tensor expression, and that has to be asserted
// rather than assumed. TensorExpansion answers every other question the
// concept asks -- it has a rank, a grid, a slot alphabet and a Represents --
// so without the truncation-degree discriminator it satisfies TensorExpr, and
// Permute, the tensor product, Materialise and the bundle maps all accept one
// while composing the wrong Component: a view over coefficients rather than a
// spin-weighted node.
//
// Found the way these things are found. Tangential(SurfaceGradient(Embed(t)))
// resolved to the *spatial* projection for a spectral operand, because a
// forwarding reference binds a prvalue better than a const reference does, and
// the error was that the resulting node had no Coefficient.
TEST(TensorAlgebra, ASpectralTensorIsNotATensorExpression) {
  using Field = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  using Expansion = TensorExpansion<2, NoSymmetry<2>, ComplexTensor, Grid>;

  static_assert(TensorExpr<Field>);
  static_assert(!TensorExpr<Expansion>);

  static_assert(Materialisable<Field>);
  static_assert(!Materialisable<Expansion>);
  SUCCEED();
}
