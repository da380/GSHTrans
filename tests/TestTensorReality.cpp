#include <gtest/gtest.h>

#include <GSHTrans/All>

#include <array>
#include <cmath>
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

auto TestGrid() { return Grid(4, 4, FFTWpp::Estimate); }

template <typename T, Int... Alphas>
void Fill(T& t, Real tag) {
  auto u = t.template Component<Alphas...>();
  const auto& grid = t.Grid();
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      if constexpr (std::same_as<typename std::remove_cvref_t<decltype(u)>::Value,
                                 RealValued>) {
        u[iTheta, iPhi] = tag + iTheta + 0.5 * iPhi;
      } else {
        u[iTheta, iPhi] = Complex{tag + iTheta, 0.5 * iPhi - tag};
      }
    }
  }
}

}  // namespace

//--------------------------------------------------------------------------//
//                          What a real tensor stores                        //
//--------------------------------------------------------------------------//

TEST(TensorReality, StorageIsTheRealDegreesOfFreedom) {
  // Two reals for each complex component and one for each the reality
  // condition pins, which is 3^p -- the real degrees of freedom of a real
  // rank-p tensor, and the theory note's claim that "the reduction loses
  // nothing".
  static_assert((TensorField<1, NoSymmetry<1>, RealTensor, Grid>::RealsPerPoint)
                == 3);
  static_assert((TensorField<2, NoSymmetry<2>, RealTensor, Grid>::RealsPerPoint)
                == 9);
  static_assert((TensorField<4, NoSymmetry<4>, RealTensor, Grid>::RealsPerPoint)
                == 81);

  // Symmetry composes with it without being special-cased.
  static_assert((TensorField<2, Symmetric<2>, RealTensor, Grid>::RealsPerPoint)
                == 6);
  static_assert(
      (TensorField<2, Antisymmetric<2>, RealTensor, Grid>::RealsPerPoint) == 3);
  static_assert((TensorField<3, Symmetric<3>, RealTensor, Grid>::RealsPerPoint)
                == 10);
  static_assert((TensorField<4, ElasticSymmetry, RealTensor, Grid>::RealsPerPoint)
                == 21);

  // Against a complex tensor of the same rank, which stores 2 * 3^p.
  static_assert(
      (TensorField<4, NoSymmetry<4>, ComplexTensor, Grid>::RealsPerPoint) ==
      162);
  SUCCEED();
}

TEST(TensorReality, PinnedComponentsAreRealFieldsAtUpperIndexZero) {
  using T = TensorField<2, NoSymmetry<2>, RealTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);

  // Four complex components and one pinned, and the buffers say so.
  static_assert(T::ComplexComponents == 4);
  static_assert(T::RealComponents == 1);
  EXPECT_EQ(t.Size(), 4 * t.FieldSize());
  EXPECT_EQ(t.RealSize(), 1 * t.FieldSize());

  // The all-zero component is the pinned one, and it is a real-valued field.
  // Phase 1 admits RealValued only at upper index zero, which is exactly
  // where a self-paired component must sit.
  auto zero = t.Component<0, 0>();
  static_assert(std::same_as<typename decltype(zero)::Value, RealValued>);
  static_assert(decltype(zero)::UpperIndex == 0);
  static_assert(std::same_as<typename decltype(zero)::Scalar, Real>);
}

TEST(TensorReality, ComplexTensorsAreUnchanged) {
  using C = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  static_assert(C::ComplexComponents == 9);
  static_assert(C::RealComponents == 0);

  auto grid = TestGrid();
  auto t = C(grid);
  EXPECT_EQ(t.RealSize(), 0);
  auto zero = t.Component<0, 0>();
  static_assert(
      std::same_as<typename decltype(zero)::Value, ComplexValued>);
}

//--------------------------------------------------------------------------//
//                            The reality condition                          //
//--------------------------------------------------------------------------//

// The relation the whole reduction rests on, checked at every component of
// every tensor below: T^{-alpha} = (-1)^N conj(T^{alpha}).
template <typename T, Int... Alphas>
void ExpectRealityCondition(const T& t, Int iTheta, Int iPhi) {
  constexpr auto negated =
      MultiIndex<T::Rank>(std::array<Int, T::Rank>{Alphas...}).Negated();
  constexpr auto N = MultiIndex<T::Rank>(std::array<Int, T::Rank>{Alphas...})
                         .UpperIndex();

  const auto value = t.template Component<Alphas...>()[iTheta, iPhi];
  const auto derived = [&]<std::size_t... I>(std::index_sequence<I...>) {
    return t.template Component<negated[I]...>()[iTheta, iPhi];
  }(std::make_index_sequence<T::Rank>{});

  const auto expected =
      static_cast<Real>(MinusOneToPower(N)) * std::conj(Complex{value});
  EXPECT_NEAR(std::real(Complex{derived}), std::real(expected), 1.0e-13);
  EXPECT_NEAR(std::imag(Complex{derived}), std::imag(expected), 1.0e-13);
}

TEST(TensorReality, EveryComponentSatisfiesTheRealityCondition) {
  using T = TensorField<2, NoSymmetry<2>, RealTensor, Grid>;
  auto grid = TestGrid();
  auto t = T(grid);

  Fill<T, -1, -1>(t, 1.0);
  Fill<T, -1, 0>(t, 2.0);
  Fill<T, -1, 1>(t, 3.0);
  Fill<T, 0, -1>(t, 4.0);
  Fill<T, 0, 0>(t, 5.0);

  const auto& tensor = t;
  ExpectRealityCondition<T, -1, -1>(tensor, 2, 3);
  ExpectRealityCondition<T, -1, 0>(tensor, 2, 3);
  ExpectRealityCondition<T, -1, 1>(tensor, 2, 3);
  ExpectRealityCondition<T, 0, -1>(tensor, 2, 3);
  ExpectRealityCondition<T, 0, 0>(tensor, 2, 3);
  ExpectRealityCondition<T, 0, 1>(tensor, 2, 3);
  ExpectRealityCondition<T, 1, -1>(tensor, 2, 3);
  ExpectRealityCondition<T, 1, 0>(tensor, 2, 3);
  ExpectRealityCondition<T, 1, 1>(tensor, 2, 3);
}

TEST(TensorReality, TheVectorCaseIsTheOneTheTheoryNoteSpellsOut) {
  // "For a vector this reads u^0 = conj(u^0) and u^- = -conj(u^+): the radial
  // component is a real field, and the two transverse components are not
  // independent."
  using V = TensorField<1, NoSymmetry<1>, RealTensor, Grid>;
  auto grid = TestGrid();
  auto v = V(grid);

  Fill<V, -1>(v, 1.0);
  Fill<V, 0>(v, 2.0);

  const auto& u = v;

  // u^0 is a real field, not merely a complex one that happens to be real.
  static_assert(
      std::same_as<typename decltype(u.Component<0>())::Value, RealValued>);

  // u^- = -conj(u^+), equivalently u^+ = -conj(u^-).
  const auto minus = u.Component<-1>()[1, 2];
  const auto plus = u.Component<1>()[1, 2];
  EXPECT_NEAR(std::real(plus), std::real(-std::conj(minus)), 1.0e-14);
  EXPECT_NEAR(std::imag(plus), std::imag(-std::conj(minus)), 1.0e-14);

  // The derived component carries the reversed upper index, which is what
  // conj does and what the layer it replaced got wrong.
  static_assert(decltype(u.Component<-1>())::UpperIndex == -1);
  static_assert(decltype(u.Component<1>())::UpperIndex == 1);
}

TEST(TensorReality, DerivedComponentsAreNotWritable) {
  using T = TensorField<2, NoSymmetry<2>, RealTensor, Grid>;

  // The stored half is writable; the half the reality condition derives is
  // not, since writing it would mean conjugating on the way in.
  static_assert(T::Writable<-1, -1>);
  static_assert(!T::Writable<1, 1>);
  static_assert(T::Represents<1, 1>);

  // The pinned component is writable as the one real number it is.
  static_assert(T::Writable<0, 0>);
}

//--------------------------------------------------------------------------//
//                      Symmetry and reality together                       //
//--------------------------------------------------------------------------//

TEST(TensorReality, SymmetricRealTensorMatchesTheWorkedExample) {
  // The theory note's worked check: four orbits, six reals, and two of the
  // four stored components real-valued -- {(00)} and the self-paired
  // {(-+), (+-)}.
  using T = TensorField<2, Symmetric<2>, RealTensor, Grid>;
  static_assert(T::StoredComponents == 4);
  static_assert(T::ComplexComponents == 2);
  static_assert(T::RealComponents == 2);
  static_assert(T::RealsPerPoint == 6);

  auto grid = TestGrid();
  auto t = T(grid);
  Fill<T, -1, -1>(t, 1.0);
  Fill<T, -1, 0>(t, 2.0);
  Fill<T, -1, 1>(t, 3.0);
  Fill<T, 0, 0>(t, 4.0);

  const auto& tensor = t;

  // (-+) is self-paired, so it is real, and it equals (+-) by symmetry.
  static_assert(std::same_as<
                typename decltype(tensor.Component<-1, 1>())::Value, RealValued>);
  EXPECT_EQ((tensor.Component<-1, 1>()[1, 1]),
            (tensor.Component<1, -1>()[1, 1]));

  // And the reality condition still holds everywhere.
  ExpectRealityCondition<T, -1, -1>(tensor, 1, 2);
  ExpectRealityCondition<T, -1, 0>(tensor, 1, 2);
  ExpectRealityCondition<T, -1, 1>(tensor, 1, 2);
  ExpectRealityCondition<T, 0, 0>(tensor, 1, 2);
  ExpectRealityCondition<T, 0, 1>(tensor, 1, 2);
}

TEST(TensorReality, AntisymmetricRealTensorHasAnImaginaryComponent) {
  // The one place an Imaginary constraint arises: a real antisymmetric rank-2
  // tensor, whose self-paired component satisfies T = -conj(T).
  using T = TensorField<2, Antisymmetric<2>, RealTensor, Grid>;
  static_assert(T::StoredComponents == 2);
  static_assert(T::ComplexComponents == 1);
  static_assert(T::RealComponents == 1);
  static_assert(T::RealsPerPoint == 3);

  auto grid = TestGrid();
  auto t = T(grid);
  Fill<T, -1, 0>(t, 1.0);
  Fill<T, -1, 1>(t, 2.0);

  const auto& tensor = t;

  // The pinned component is stored as the real coefficient of i, so its value
  // is purely imaginary.
  const auto value = tensor.Component<-1, 1>()[1, 1];
  EXPECT_NEAR(std::real(value), 0.0, 1.0e-14);
  EXPECT_NE(std::imag(value), 0.0);

  ExpectRealityCondition<T, -1, 0>(tensor, 1, 1);
  ExpectRealityCondition<T, -1, 1>(tensor, 1, 1);
}

//--------------------------------------------------------------------------//
//                              The transform                                //
//--------------------------------------------------------------------------//

TEST(TensorReality, RoundTripsWithTheRealComponentsOnTheRealPath) {
  using T = TensorField<2, Symmetric<2>, RealTensor, Grid>;
  constexpr auto lMax = Int{6};

  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto t = T(grid);

  // A pinned component uses the reduced m >= 0 storage, so the coefficient
  // count is not the complex count times the number of components.
  const auto expected =
      static_cast<Int>(grid.CoefficientSize(lMax, -2)) +
      static_cast<Int>(grid.CoefficientSize(lMax, -1)) +
      2 * static_cast<Int>(grid.RealCoefficientSize(lMax));
  EXPECT_EQ(t.CoefficientSize(lMax), expected);

  // The round trip starts from a *field*, not from arbitrary coefficients.
  // A real field's m = 0 coefficients are themselves real -- f_{l,-m} =
  // (-1)^m conj(f_{lm}) at m = 0 -- so an arbitrary complex block is not the
  // transform of anything the reduced storage can hold, and filling one would
  // be testing the library against invalid data.
  Fill<T, -1, -1>(t, 1.0);
  Fill<T, -1, 0>(t, 2.0);
  Fill<T, -1, 1>(t, 3.0);
  Fill<T, 0, 0>(t, 4.0);

  auto coefficients = std::vector<Complex>(t.CoefficientSize(lMax));
  t.ForwardTransformation(lMax, coefficients);
  t.InverseTransformation(lMax, coefficients);

  // The reality condition survives a transform, which is the check that the
  // pinned components really did go through the real path: a complex
  // transform of the same data would not keep them real.
  const auto& tensor = t;
  ExpectRealityCondition<T, -1, -1>(tensor, 2, 2);
  ExpectRealityCondition<T, -1, 0>(tensor, 2, 2);
  ExpectRealityCondition<T, -1, 1>(tensor, 2, 2);
  ExpectRealityCondition<T, 0, 0>(tensor, 2, 2);

  auto back = std::vector<Complex>(coefficients.size());
  t.ForwardTransformation(lMax, back);
  for (auto i = std::size_t{0}; i < coefficients.size(); i++) {
    EXPECT_NEAR(back[i].real(), coefficients[i].real(), 1.0e-11) << "at " << i;
    EXPECT_NEAR(back[i].imag(), coefficients[i].imag(), 1.0e-11) << "at " << i;
  }
}

//--------------------------------------------------------------------------//
//                        Producing a real tensor                            //
//--------------------------------------------------------------------------//

// The payoff of the reduction for a computed result: a real tensor built from
// an expression stores half as much and derives the rest.
TEST(TensorReality, MaterialiseCanProduceARealTensor) {
  using Source = TensorField<2, NoSymmetry<2>, RealTensor, Grid>;
  auto grid = TestGrid();
  auto s = Source(grid);
  Fill<Source, -1, -1>(s, 1.0);
  Fill<Source, -1, 0>(s, 2.0);
  Fill<Source, -1, 1>(s, 3.0);
  Fill<Source, 0, -1>(s, 4.0);
  Fill<Source, 0, 0>(s, 5.0);

  const auto& tensor = s;

  // The transpose of a real tensor is a real tensor, which is a fact about
  // the value that only the caller can assert.
  auto transposed = Materialise<NoSymmetry<2>, RealTensor>(Transpose(tensor));
  static_assert(std::same_as<decltype(transposed),
                             TensorField<2, NoSymmetry<2>, RealTensor, Grid>>);
  static_assert(decltype(transposed)::RealsPerPoint == 9);

  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      EXPECT_NEAR(std::real(Complex{transposed.Component<0, 1>()[iTheta, iPhi]}),
                  std::real(Complex{tensor.Component<1, 0>()[iTheta, iPhi]}),
                  1.0e-13);
      EXPECT_NEAR(std::imag(Complex{transposed.Component<0, 1>()[iTheta, iPhi]}),
                  std::imag(Complex{tensor.Component<1, 0>()[iTheta, iPhi]}),
                  1.0e-13);
    }
  }

  // And the result is a real tensor in its own right, so its derived half
  // still satisfies the reality condition.
  const auto& result = transposed;
  ExpectRealityCondition<decltype(transposed), -1, 0>(result, 1, 1);
  ExpectRealityCondition<decltype(transposed), 0, 0>(result, 1, 1);
}

TEST(TensorReality, TheSymmetricPartOfARealTensorIsRealAndSymmetric) {
  using Source = TensorField<2, NoSymmetry<2>, RealTensor, Grid>;
  auto grid = TestGrid();
  auto s = Source(grid);
  Fill<Source, -1, -1>(s, 1.0);
  Fill<Source, -1, 0>(s, 2.0);
  Fill<Source, -1, 1>(s, 3.0);
  Fill<Source, 0, -1>(s, 4.0);
  Fill<Source, 0, 0>(s, 5.0);

  const auto& tensor = s;
  auto sym = Materialise<Symmetric<2>, RealTensor>(
      Symmetrise<Symmetric<2>>(tensor));

  // Six reals a point rather than eighteen, which is what a rank-2 tensor
  // costs stored naively as nine complex fields.
  static_assert(decltype(sym)::RealsPerPoint == 6);
  static_assert(decltype(sym)::StoredComponents == 4);

  const auto& result = sym;
  const auto expected = 0.5 * (Complex{tensor.Component<0, 1>()[1, 2]} +
                               Complex{tensor.Component<1, 0>()[1, 2]});
  EXPECT_NEAR(std::real(Complex{result.Component<0, 1>()[1, 2]}),
              std::real(expected), 1.0e-13);
  EXPECT_NEAR(std::imag(Complex{result.Component<0, 1>()[1, 2]}),
              std::imag(expected), 1.0e-13);

  // Symmetric and real at once: transposing does nothing and the reality
  // condition holds.
  EXPECT_EQ((result.Component<0, 1>()[1, 2]), (result.Component<1, 0>()[1, 2]));
  ExpectRealityCondition<decltype(sym), -1, 0>(result, 1, 2);
  ExpectRealityCondition<decltype(sym), -1, 1>(result, 1, 2);
}
