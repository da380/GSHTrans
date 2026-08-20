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
