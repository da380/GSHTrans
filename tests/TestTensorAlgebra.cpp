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
  EXPECT_EQ((transposed.Component<0, 1>()[2, 3]), (tensor.Component<1, 0>()[2, 3]));
  EXPECT_EQ((transposed.Component<1, 0>()[2, 3]), (tensor.Component<0, 1>()[2, 3]));
  EXPECT_NE((transposed.Component<0, 1>()[2, 3]), (tensor.Component<0, 1>()[2, 3]));

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
