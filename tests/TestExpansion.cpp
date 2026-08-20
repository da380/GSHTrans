#include <gtest/gtest.h>

#include <GSHTrans/All>

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

}  // namespace

//--------------------------------------------------------------------------//
//                             The expansion                                 //
//--------------------------------------------------------------------------//

TEST(SpinExpansion, HoldsTheCoefficientsItsUpperIndexAllows) {
  constexpr auto lMax = Int{5};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  auto e = SpinExpansion<2, Grid>(grid, lMax);
  static_assert(decltype(e)::UpperIndex == 2);

  // The degrees start at |N|: below that the harmonics do not exist, so
  // there is no coefficient to hold.
  EXPECT_EQ(e.MinDegree(), 2);
  EXPECT_EQ(e.MaxDegree(), lMax);
  EXPECT_EQ(e.Size(), static_cast<Int>(grid.CoefficientSize(lMax, 2)));

  // Asking for an expansion below its own upper index is an error rather
  // than an empty one.
  EXPECT_THROW((SpinExpansion<2, Grid>(grid, 1)), std::invalid_argument);
}

TEST(SpinExpansion, RealFieldsUseTheReducedStorage) {
  constexpr auto lMax = Int{5};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  auto complexOne = SpinExpansion<0, Grid, ComplexValued>(grid, lMax);
  auto realOne = SpinExpansion<0, Grid, RealValued>(grid, lMax);

  // The reduced m >= 0 storage is a type distinction rather than a
  // convention, and it is about half the size.
  static_assert(std::same_as<decltype(complexOne)::MRange, All>);
  static_assert(std::same_as<decltype(realOne)::MRange, NonNegative>);
  EXPECT_EQ(realOne.Size(), static_cast<Int>(grid.RealCoefficientSize(lMax)));
  EXPECT_LT(realOne.Size(), complexOne.Size());

  // And it is available only at upper index zero, which is the same
  // constraint phase 1 puts on a real-valued field. That one is a
  // static_assert in the class body rather than a requires-clause, since
  // there is no overload to fall through to -- so it is a hard error and not
  // something a negative test can probe.
}

TEST(SpinExpansion, IndexingReachesTheCoefficientTheTransformWrote) {
  constexpr auto lMax = Int{4};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  auto e = SpinExpansion<1, Grid>(grid, lMax);
  for (auto l : e.Degrees()) {
    for (auto m : e.Orders(l)) {
      e[l, m] = Complex{static_cast<Real>(l), static_cast<Real>(m)};
    }
  }

  // The flat layout is the transform's, so a raw read agrees with an
  // indexed one.
  auto indices = GSHIndices<All>(lMax, lMax, 1);
  for (auto l : e.Degrees()) {
    for (auto m : e.Orders(l)) {
      EXPECT_EQ((e[l, m]), e.Data()[indices.Index(l, m)]);
      EXPECT_EQ((e[l, m]), (Complex{static_cast<Real>(l),
                                    static_cast<Real>(m)}));
    }
  }
}

TEST(SpinExpansion, RoundTripsAgainstTheField) {
  constexpr auto lMax = Int{6};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  auto e = SpinExpansion<2, Grid>(grid, lMax);
  for (auto l : e.Degrees()) {
    for (auto m : e.Orders(l)) {
      e[l, m] = Complex{std::cos(0.3 * (l + m)), std::sin(0.7 * l - m)};
    }
  }

  auto field = Evaluate(e);
  static_assert(decltype(field)::UpperIndex == 2);
  static_assert(SpinWeighted<decltype(field)>);

  auto back = Expand(field, lMax);
  for (auto l : e.Degrees()) {
    for (auto m : e.Orders(l)) {
      EXPECT_NEAR((back[l, m]).real(), (e[l, m]).real(), 1.0e-11);
      EXPECT_NEAR((back[l, m]).imag(), (e[l, m]).imag(), 1.0e-11);
    }
  }
}

TEST(SpinExpansion, ExpandsAnExpressionAsReadilyAsAField) {
  constexpr auto lMax = Int{5};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  auto u = SpinField<1, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::cos(theta), std::sin(phi)};
  });
  auto v = SpinField<1, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::sin(theta) * std::cos(phi), 0.5};
  });

  // Expand takes any spin-weighted node, so a lazy expression needs no
  // materialising first -- which is what EvaluateInto was the seam for.
  const auto& a = u;
  const auto& b = v;
  auto sum = Expand(a + b, lMax);
  static_assert(decltype(sum)::UpperIndex == 1);

  auto separately = Expand(a, lMax);
  auto other = Expand(b, lMax);
  for (auto l : sum.Degrees()) {
    for (auto m : sum.Orders(l)) {
      const auto expected = separately[l, m] + other[l, m];
      EXPECT_NEAR((sum[l, m]).real(), expected.real(), 1.0e-12);
      EXPECT_NEAR((sum[l, m]).imag(), expected.imag(), 1.0e-12);
    }
  }
}

//--------------------------------------------------------------------------//
//                       Raising and lowering the index                      //
//--------------------------------------------------------------------------//

TEST(Eth, ChangesTheUpperIndexByOne) {
  constexpr auto lMax = Int{6};
  auto grid = Grid(lMax, 3, FFTWpp::Estimate);

  auto e = SpinExpansion<1, Grid>(grid, lMax);
  auto raised = Raise(e);
  auto lowered = Lower(e);

  static_assert(decltype(raised)::UpperIndex == 2);
  static_assert(decltype(lowered)::UpperIndex == 0);

  // Raising narrows the degree range and lowering widens it, both because a
  // field of upper index N has no content below degree |N|.
  EXPECT_EQ(raised.MinDegree(), 2);
  EXPECT_EQ(lowered.MinDegree(), 0);
}

// The identity that pins the relative sign of the two operators and both of
// their magnitudes: eth-bar eth is the surface Laplacian on a scalar, whose
// eigenvalue on Y_{lm} is -l(l+1).
TEST(Eth, LoweringARaisedScalarIsTheSurfaceLaplacian) {
  constexpr auto lMax = Int{7};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  auto f = SpinExpansion<0, Grid>(grid, lMax);
  for (auto l : f.Degrees()) {
    for (auto m : f.Orders(l)) {
      f[l, m] = Complex{std::cos(0.4 * l + m), std::sin(0.2 * l - m)};
    }
  }

  auto laplacian = Lower(Raise(f));
  static_assert(decltype(laplacian)::UpperIndex == 0);

  for (auto l : f.Degrees()) {
    for (auto m : f.Orders(l)) {
      const auto expected =
          -static_cast<Real>(l * (l + 1)) * Complex{f[l, m]};
      EXPECT_NEAR((laplacian[l, m]).real(), expected.real(), 1.0e-12)
          << "l = " << l << ", m = " << m;
      EXPECT_NEAR((laplacian[l, m]).imag(), expected.imag(), 1.0e-12)
          << "l = " << l << ", m = " << m;
    }
  }
}

// The other identity that survives the convention: [eth, eth-bar] = -2N.
TEST(Eth, TheCommutatorIsMinusTwiceTheUpperIndex) {
  constexpr auto lMax = Int{6};
  auto grid = Grid(lMax, 3, FFTWpp::Estimate);

  constexpr auto N = Int{1};
  auto f = SpinExpansion<N, Grid>(grid, lMax);
  for (auto l : f.Degrees()) {
    for (auto m : f.Orders(l)) {
      f[l, m] = Complex{std::sin(0.5 * l + m), std::cos(0.3 * l - m)};
    }
  }

  auto raiseThenLower = Lower(Raise(f));
  auto lowerThenRaise = Raise(Lower(f));

  for (auto l : f.Degrees()) {
    for (auto m : f.Orders(l)) {
      const auto commutator =
          Complex{lowerThenRaise[l, m]} - Complex{raiseThenLower[l, m]};
      const auto expected = -2.0 * static_cast<Real>(N) * Complex{f[l, m]};
      EXPECT_NEAR(commutator.real(), expected.real(), 1.0e-11)
          << "l = " << l << ", m = " << m;
      EXPECT_NEAR(commutator.imag(), expected.imag(), 1.0e-11)
          << "l = " << l << ", m = " << m;
    }
  }
}

// The coefficient the raised field cannot carry was zero anyway, because the
// factor vanishes exactly where the target expansion has no room for it.
TEST(Eth, TheFactorVanishesWhereTheTargetHasNoRoom) {
  constexpr auto lMax = Int{5};
  auto grid = Grid(lMax, 3, FFTWpp::Estimate);

  auto f = SpinExpansion<2, Grid>(grid, lMax);
  for (auto l : f.Degrees()) {
    for (auto m : f.Orders(l)) f[l, m] = Complex{1.0, 1.0};
  }

  // Raising from N = 2 loses degree 2, and the factor there is
  // -sqrt((2 - 2)(2 + 3)) = 0, so nothing was lost.
  EXPECT_EQ(EthDetails::RaisingFactor<Real>(2, 2), 0.0);
  auto raised = Raise(f);
  EXPECT_EQ(raised.MinDegree(), 3);

  // Lowering from N = 0 likewise: the l = 0 coefficient of the source cannot
  // reach the target at N = -1, and its factor is zero.
  EXPECT_EQ(EthDetails::LoweringFactor<Real>(0, 0), 0.0);
}

TEST(Eth, LoweringAScalarReadsTheOrdersARealFieldDoesNotStore) {
  constexpr auto lMax = Int{5};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  // A real scalar stores only m >= 0; the negative orders follow from
  // f_{l,-m} = (-1)^m conj(f_{lm}), and lowering has to use them rather than
  // read off the end.
  auto real = SpinExpansion<0, Grid, RealValued>(grid, lMax);
  auto complexOne = SpinExpansion<0, Grid, ComplexValued>(grid, lMax);
  for (auto l : real.Degrees()) {
    for (auto m : real.Orders(l)) {
      const auto value = Complex{std::cos(0.4 * l + m), std::sin(0.6 * m)};
      real[l, m] = m == 0 ? Complex{value.real(), 0.0} : value;
    }
  }
  for (auto l : complexOne.Degrees()) {
    for (auto m : complexOne.Orders(l)) {
      complexOne[l, m] =
          m >= 0 ? Complex{real[l, m]}
                 : static_cast<Real>(MinusOneToPower(m)) *
                       std::conj(Complex{real[l, -m]});
    }
  }

  auto fromReal = Lower(real);
  auto fromComplex = Lower(complexOne);
  for (auto l : fromReal.Degrees()) {
    for (auto m : fromReal.Orders(l)) {
      EXPECT_NEAR((fromReal[l, m]).real(), (fromComplex[l, m]).real(), 1.0e-13)
          << "l = " << l << ", m = " << m;
      EXPECT_NEAR((fromReal[l, m]).imag(), (fromComplex[l, m]).imag(), 1.0e-13)
          << "l = " << l << ", m = " << m;
    }
  }
}

//--------------------------------------------------------------------------//
//                       A tensor in the spectral domain                     //
//--------------------------------------------------------------------------//

TEST(TensorExpansion, MirrorsTheFieldWithoutTheSecondBuffer) {
  constexpr auto lMax = Int{5};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  using Field = TensorField<2, Symmetric<2>, RealTensor, Grid>;
  using Expansion = TensorExpansion<2, Symmetric<2>, RealTensor, Grid>;

  auto field = Field(grid);
  auto e = Expansion(grid, lMax);

  // The same stored set as the field -- reality reduces both sides.
  static_assert(Expansion::StoredComponents == Field::StoredComponents);
  static_assert(Expansion::RealComponents == 2);

  // But one buffer, not two: a pinned component is a real *field*, and its
  // coefficients are complex numbers in the reduced m >= 0 storage.
  EXPECT_EQ(e.Size(), field.CoefficientSize(lMax));

  // Its block is about half the length of a complex one at the same degree,
  // which is the same saving in the spectral domain as in the spatial one.
  auto pinned = e.Component<0, 0>();
  static_assert(std::same_as<decltype(pinned)::Value, RealValued>);
  EXPECT_EQ(pinned.Size(), static_cast<Int>(grid.RealCoefficientSize(lMax)));
  EXPECT_LT(pinned.Size(), static_cast<Int>(grid.CoefficientSize(lMax, 0)));
}

TEST(TensorExpansion, ComponentsCarryTheirOwnUpperIndex) {
  constexpr auto lMax = Int{4};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto e = TensorExpansion<2, NoSymmetry<2>, ComplexTensor, Grid>(grid, lMax);

  static_assert(decltype(e.Component<1, 1>())::UpperIndex == 2);
  static_assert(decltype(e.Component<0, 1>())::UpperIndex == 1);
  static_assert(decltype(e.Component<-1, 1>())::UpperIndex == 0);

  // The degrees of a block start at the component's own upper index.
  EXPECT_EQ((e.Component<1, 1>().MinDegree()), 2);
  EXPECT_EQ((e.Component<-1, 1>().MinDegree()), 0);

  // Writing through a component reaches the buffer the transform uses.
  e.Component<1, 1>()[3, -2] = Complex{1.0, -2.0};
  const auto& expansion = e;
  EXPECT_EQ((expansion.Component<1, 1>()[3, -2]), (Complex{1.0, -2.0}));
}

TEST(TensorExpansion, RoundTripsAgainstTheTensorField) {
  constexpr auto lMax = Int{6};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  using Field = TensorField<2, Symmetric<2>, RealTensor, Grid>;
  auto field = Field(grid);

  // Fill through the components, so that whatever the storage is, the values
  // are ones the tensor can hold.
  const auto write = [&](auto&& u, Real tag) {
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        using Node = std::remove_cvref_t<decltype(u)>;
        if constexpr (std::same_as<typename Node::Value, RealValued>) {
          u[iTheta, iPhi] = tag + iTheta + 0.25 * iPhi;
        } else {
          u[iTheta, iPhi] = Complex{tag + iTheta, 0.25 * iPhi - tag};
        }
      }
    }
  };
  write(field.Component<-1, -1>(), 1.0);
  write(field.Component<-1, 0>(), 2.0);
  write(field.Component<-1, 1>(), 3.0);
  write(field.Component<0, 0>(), 4.0);

  auto e = Expand(field, lMax);
  static_assert(std::same_as<decltype(e),
                             TensorExpansion<2, Symmetric<2>, RealTensor,
                                             Grid>>);

  auto back = Evaluate(e);
  static_assert(std::same_as<decltype(back), Field>);

  // Once band-limited, the two representations agree.
  auto again = Expand(back, lMax);
  for (auto i = Int{0}; i < e.Size(); i++) {
    EXPECT_NEAR(again.Data()[i].real(), e.Data()[i].real(), 1.0e-11)
        << "at " << i;
    EXPECT_NEAR(again.Data()[i].imag(), e.Data()[i].imag(), 1.0e-11)
        << "at " << i;
  }

  // And the reality condition survives the round trip.
  const auto& result = back;
  const auto value = Complex{result.Component<-1, -1>()[2, 2]};
  const auto derived = Complex{result.Component<1, 1>()[2, 2]};
  EXPECT_NEAR(derived.real(), std::conj(value).real(), 1.0e-11);
  EXPECT_NEAR(derived.imag(), std::conj(value).imag(), 1.0e-11);
}

TEST(TensorExpansion, DerivedComponentsAreNotOfferedInTheSpectralDomain) {
  using E = TensorExpansion<2, NoSymmetry<2>, RealTensor, Grid>;

  // A derived component has no block: deriving it here means applying
  // eq:complevel, T^{-N}_{l,-m} = (-1)^m conj(T^N_{lm}), which reverses the
  // order index rather than acting pointwise. That is not a view over
  // anything, so the accessor is not offered.
  static_assert(E::Writable<-1, -1>);
  static_assert(!E::Writable<1, 1>);
  SUCCEED();
}

// Every component is readable in the spectral domain, stored or not. The
// oracle is the same tensor widened to a complex one, where every component
// *is* stored -- so the reduced expansion's derived components are checked
// against directly stored ones rather than against the relation they were
// computed from.
TEST(TensorExpansion, EveryComponentIsReadable) {
  constexpr auto lMax = Int{6};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  using RealField = TensorField<2, NoSymmetry<2>, RealTensor, Grid>;
  auto t = RealField(grid);

  const auto write = [&](auto&& u, Real tag) {
    using Node = std::remove_cvref_t<decltype(u)>;
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        if constexpr (std::same_as<typename Node::Value, RealValued>) {
          u[iTheta, iPhi] = tag + std::cos(0.3 * iTheta) * std::sin(0.2 * iPhi);
        } else {
          u[iTheta, iPhi] = Complex{tag + std::cos(0.3 * iTheta),
                                    std::sin(0.2 * iPhi) - tag};
        }
      }
    }
  };
  write(t.Component<-1, -1>(), 1.0);
  write(t.Component<-1, 0>(), 2.0);
  write(t.Component<-1, 1>(), 3.0);
  write(t.Component<0, -1>(), 4.0);
  write(t.Component<0, 0>(), 5.0);

  const auto& tensor = t;
  auto reduced = Expand(tensor, lMax);

  // The same field as a complex tensor: nine stored components, nothing
  // derived.
  auto widened = Materialise<NoSymmetry<2>, ComplexTensor>(tensor);
  auto full = Expand(widened, lMax);

  const auto compare = [&]<Int A, Int B>() {
    for (auto l = Int{0}; l <= lMax; l++) {
      for (auto m = -l; m <= l; m++) {
        const auto got = reduced.Coefficient<A, B>(l, m);
        const auto expected = full.Coefficient<A, B>(l, m);
        EXPECT_NEAR(got.real(), expected.real(), 1.0e-11)
            << "component (" << A << "," << B << ") at l = " << l
            << ", m = " << m;
        EXPECT_NEAR(got.imag(), expected.imag(), 1.0e-11)
            << "component (" << A << "," << B << ") at l = " << l
            << ", m = " << m;
      }
    }
  };

  compare.template operator()<-1, -1>();
  compare.template operator()<-1, 0>();
  compare.template operator()<-1, 1>();
  compare.template operator()<0, -1>();
  compare.template operator()<0, 0>();
  compare.template operator()<0, 1>();   // derived by reality
  compare.template operator()<1, -1>();  // derived
  compare.template operator()<1, 0>();   // derived
  compare.template operator()<1, 1>();   // derived
}

TEST(TensorExpansion, CoefficientsVanishBelowTheirOwnDegree) {
  constexpr auto lMax = Int{5};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto e = TensorExpansion<2, NoSymmetry<2>, ComplexTensor, Grid>(grid, lMax);

  // A component at upper index 2 has no content below degree 2: not a
  // missing value but an absent one.
  EXPECT_EQ((e.Coefficient<1, 1>(0, 0)), Complex{});
  EXPECT_EQ((e.Coefficient<1, 1>(1, 0)), Complex{});
  EXPECT_EQ((e.Coefficient<1, 1>(6, 0)), Complex{});  // above lMax
  EXPECT_EQ((e.Coefficient<0, 0>(2, 3)), Complex{});  // |m| > l
}

TEST(TensorExpansion, SymmetryRelativesAgreeWithTheirRepresentative) {
  constexpr auto lMax = Int{5};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  auto skew = TensorExpansion<2, Antisymmetric<2>, ComplexTensor, Grid>(grid,
                                                                       lMax);
  skew.Component<-1, 0>()[3, 1] = Complex{2.0, -1.0};

  const auto& e = skew;
  EXPECT_EQ((e.Coefficient<-1, 0>(3, 1)), (Complex{2.0, -1.0}));
  EXPECT_EQ((e.Coefficient<0, -1>(3, 1)), (Complex{-2.0, 1.0}));
  EXPECT_EQ((e.Coefficient<0, 0>(3, 1)), Complex{});  // vanishing orbit
}

//--------------------------------------------------------------------------//
//                       The contravariant derivative                        //
//--------------------------------------------------------------------------//

// A scalar has no slots and therefore no connection terms, so there the
// operator *is* eth, up to the sqrt(2) that is the normalisation of e_{+-}.
// This checks the Omega factors and the plumbing; it says nothing about the
// connection terms, which is what the tests after it are for.
TEST(ContravariantDerivative, OnAScalarItIsEthUpToRootTwo) {
  constexpr auto lMax = Int{8};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  auto scalar = TensorExpansion<0, NoSymmetry<0>, ComplexTensor, Grid>(grid,
                                                                       lMax);
  auto asSpin = SpinExpansion<0, Grid>(grid, lMax);
  for (auto l : asSpin.Degrees()) {
    for (auto m : asSpin.Orders(l)) {
      const auto value = Complex{std::cos(0.3 * l + m), std::sin(0.2 * l - m)};
      asSpin[l, m] = value;
      scalar.Component<>()[l, m] = value;
    }
  }

  auto gradient = SurfaceGradient(scalar);
  static_assert(decltype(gradient)::Rank == 1);

  auto raised = Raise(asSpin);
  auto lowered = Lower(asSpin);
  const auto rootTwo = std::numbers::sqrt2_v<Real>;

  for (auto l = Int{1}; l <= lMax; l++) {
    for (auto m = -l; m <= l; m++) {
      // d^+ = -eth / sqrt(2), d^- = +eth-bar / sqrt(2).
      const auto plus = gradient.Coefficient<1>(l, m);
      const auto minus = gradient.Coefficient<-1>(l, m);
      EXPECT_NEAR(plus.real(), (-Complex{raised[l, m]} / rootTwo).real(),
                  1.0e-12) << "l = " << l << ", m = " << m;
      EXPECT_NEAR(minus.real(), (Complex{lowered[l, m]} / rootTwo).real(),
                  1.0e-12) << "l = " << l << ", m = " << m;
      EXPECT_NEAR(plus.imag(), (-Complex{raised[l, m]} / rootTwo).imag(),
                  1.0e-12);
      EXPECT_NEAR(minus.imag(), (Complex{lowered[l, m]} / rootTwo).imag(),
                  1.0e-12);

      // The surface gradient has no radial component at all.
      EXPECT_EQ((gradient.Coefficient<0>(l, m)), Complex{});
    }
  }
}

// The check that the connection terms are right, and it is not circular: the
// metric trace of the second surface gradient of a scalar is the surface
// Laplacian, whose eigenvalue on Y_{lm} is -l(l+1).
//
// The connection term enters through the -v^0 subtraction in d^- v^+ and
// d^+ v^-, so this would fail if the operator were eth applied component by
// component.
TEST(ContravariantDerivative, TheTraceOfTheSecondGradientIsTheLaplacian) {
  constexpr auto lMax = Int{8};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  auto scalar = TensorExpansion<0, NoSymmetry<0>, ComplexTensor, Grid>(grid,
                                                                       lMax);
  for (auto l = Int{0}; l <= lMax; l++) {
    for (auto m = -l; m <= l; m++) {
      scalar.Component<>()[l, m] =
          Complex{std::sin(0.4 * l + m), std::cos(0.15 * l - m)};
    }
  }

  auto second = SurfaceGradient(SurfaceGradient(scalar));
  static_assert(decltype(second)::Rank == 2);

  // The metric contraction g_{ab} = (-1)^a delta_{a+b,0}, so the trace is
  // -T^{-+} + T^{00} - T^{+-}; the middle term vanishes because the surface
  // gradient has no radial component.
  for (auto l = Int{0}; l <= lMax; l++) {
    for (auto m = -l; m <= l; m++) {
      const auto trace = -second.Coefficient<-1, 1>(l, m) +
                         second.Coefficient<0, 0>(l, m) -
                         second.Coefficient<1, -1>(l, m);
      const auto expected = -static_cast<Real>(l * (l + 1)) *
                            scalar.Coefficient<>(l, m);
      EXPECT_NEAR(trace.real(), expected.real(), 1.0e-11)
          << "l = " << l << ", m = " << m;
      EXPECT_NEAR(trace.imag(), expected.imag(), 1.0e-11)
          << "l = " << l << ", m = " << m;
    }
  }
}

// The chain rule, D&T (C.155), which is the check of the connection terms at
// rank one and above.
//
// It has to be done in the *spatial* domain. The relation holds for the
// operator on fields and does not hold coefficient by coefficient, because a
// product of fields is not a product of coefficients -- which is the same
// fact that makes "the gradient of a product" an operation in two
// representations.
TEST(ContravariantDerivative, ObeysTheChainRuleOnFields) {
  constexpr auto lMax = Int{8};
  constexpr auto band = Int{3};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  using Scalar = TensorExpansion<0, NoSymmetry<0>, ComplexTensor, Grid>;
  using Vector = TensorExpansion<1, NoSymmetry<1>, ComplexTensor, Grid>;

  // Band-limited inputs, so that the products below stay inside the grid and
  // every transform in this test is exact.
  auto f = Scalar(grid, lMax);
  auto u = Vector(grid, lMax);
  for (auto l = Int{0}; l <= band; l++) {
    for (auto m = -l; m <= l; m++) {
      f.Component<>()[l, m] = Complex{std::cos(0.7 * l + m), 0.3 * m};
    }
  }
  for (auto alpha : {Int{-1}, Int{0}, Int{1}}) {
    for (auto l = std::abs(alpha); l <= band; l++) {
      for (auto m = -l; m <= l; m++) {
        const auto value = Complex{0.5 * l - m, std::sin(0.4 * l + alpha)};
        if (alpha == -1) u.Component<-1>()[l, m] = value;
        if (alpha == 0) u.Component<0>()[l, m] = value;
        if (alpha == 1) u.Component<1>()[l, m] = value;
      }
    }
  }

  auto fField = Evaluate(f);
  auto uField = Evaluate(u);

  // S = f u, formed pointwise.
  auto s = TensorField<1, NoSymmetry<1>, ComplexTensor, Grid>(grid);
  const auto& scalarField = fField;
  const auto& vectorField = uField;
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      const auto scale = scalarField.Component<>()[iTheta, iPhi];
      s.Component<-1>()[iTheta, iPhi] =
          scale * vectorField.Component<-1>()[iTheta, iPhi];
      s.Component<0>()[iTheta, iPhi] =
          scale * vectorField.Component<0>()[iTheta, iPhi];
      s.Component<1>()[iTheta, iPhi] =
          scale * vectorField.Component<1>()[iTheta, iPhi];
    }
  }

  // The gradient of the product, and the two gradients it should decompose
  // into, all brought back to the sphere.
  auto gradS = Evaluate(SurfaceGradient(Expand(s, lMax)));
  auto gradF = Evaluate(SurfaceGradient(f));
  auto gradU = Evaluate(SurfaceGradient(u));

  const auto& product = gradS;
  const auto& scalarGradient = gradF;
  const auto& vectorGradient = gradU;

  const auto check = [&]<Int Sigma, Int Alpha>() {
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        const auto got = product.template Component<Sigma, Alpha>()[iTheta, iPhi];
        const auto expected =
            scalarGradient.template Component<Sigma>()[iTheta, iPhi] *
                vectorField.template Component<Alpha>()[iTheta, iPhi] +
            scalarField.template Component<>()[iTheta, iPhi] *
                vectorGradient.template Component<Sigma, Alpha>()[iTheta, iPhi];
        ASSERT_NEAR(got.real(), expected.real(), 1.0e-10)
            << "sigma = " << Sigma << ", alpha = " << Alpha;
        ASSERT_NEAR(got.imag(), expected.imag(), 1.0e-10)
            << "sigma = " << Sigma << ", alpha = " << Alpha;
      }
    }
  };

  check.template operator()<-1, -1>();
  check.template operator()<-1, 0>();
  check.template operator()<-1, 1>();
  check.template operator()<1, -1>();
  check.template operator()<1, 0>();
  check.template operator()<1, 1>();
}

// The gradient of a real tensor is real, so only the reduced set is computed
// and the rest follows. The oracle is the same field widened to a complex
// tensor, where nothing is derived on either side.
TEST(ContravariantDerivative, AgreesOnARealTensorAndItsWidening) {
  constexpr auto lMax = Int{7};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  using RealVector = TensorField<1, NoSymmetry<1>, RealTensor, Grid>;
  auto v = RealVector(grid);

  const auto write = [&](auto&& u, Real tag) {
    using Node = std::remove_cvref_t<decltype(u)>;
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        if constexpr (std::same_as<typename Node::Value, RealValued>) {
          u[iTheta, iPhi] = tag + std::cos(0.3 * iTheta) * std::sin(0.2 * iPhi);
        } else {
          u[iTheta, iPhi] = Complex{tag + std::cos(0.3 * iTheta),
                                    std::sin(0.2 * iPhi) - tag};
        }
      }
    }
  };
  write(v.Component<-1>(), 1.0);
  write(v.Component<0>(), 2.0);

  const auto& vector = v;
  auto reduced = SurfaceGradient(Expand(vector, lMax));
  static_assert(std::same_as<decltype(reduced)::Reality, RealTensor>);

  auto widened = Materialise<NoSymmetry<1>, ComplexTensor>(vector);
  auto full = SurfaceGradient(Expand(widened, lMax));

  const auto compare = [&]<Int A, Int B>() {
    for (auto l = Int{0}; l <= lMax; l++) {
      for (auto m = -l; m <= l; m++) {
        const auto got = reduced.Coefficient<A, B>(l, m);
        const auto expected = full.Coefficient<A, B>(l, m);
        ASSERT_NEAR(got.real(), expected.real(), 1.0e-11)
            << "(" << A << "," << B << ") at l = " << l << ", m = " << m;
        ASSERT_NEAR(got.imag(), expected.imag(), 1.0e-11)
            << "(" << A << "," << B << ") at l = " << l << ", m = " << m;
      }
    }
  };

  compare.template operator()<-1, -1>();
  compare.template operator()<-1, 0>();
  compare.template operator()<-1, 1>();
  compare.template operator()<0, -1>();
  compare.template operator()<0, 0>();
  compare.template operator()<0, 1>();
  compare.template operator()<1, -1>();
  compare.template operator()<1, 0>();
  compare.template operator()<1, 1>();

  // Half the storage, as the reduction promises.
  EXPECT_LT(reduced.Size(), full.Size());
}
