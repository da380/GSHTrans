#include <gtest/gtest.h>

#include <GSHTrans/GSHTrans.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

// Wigner tables at large degree, and in single precision.
//
// The recursion in degree is seeded at l0 = max(|n|, |m|) with a value that is
// exponentially small near the poles, and the column it starts grows back to
// order one only at its turning point. A seed that has underflowed *and* still
// matters is a column silently wrong by order one, so a table whose shape
// allows one is refused: SafeDegree finds the first degree at which a seed
// below min / eps meets its turning point, over the table's columns and all
// colatitudes. With every order and upper index zero that is 1822 in double;
// with the orders restricted it rises steeply; with large upper indices it
// falls, to 913 at |n| = |m| = 336, the worst shape of all.
//
// Single precision is recursed in double and narrowed, because what single
// precision is good for is what is stored and moved, and the recursion runs
// once. So its limit is the double one and not 194.
//
// GSHTRANS_TEST_THOROUGH=1 widens the sets, as in the 3-j tests.

using namespace GSHTrans;

namespace {

using Int = std::ptrdiff_t;

bool Thorough() {
  const auto* value = std::getenv("GSHTRANS_TEST_THOROUGH");
  return value != nullptr && std::string(value) != "0";
}

// The defect of sum_m |d^l_{nm}|^2 = 1 at one degree, from the stored values
// sqrt((2l+1)/(4 pi)) d^l_{nm}.
template <typename Table>
double UnitarityDefect(const Table& table, Int l) {
  auto sum = static_cast<long double>(0);
  for (auto m = -l; m <= l; m++) {
    const auto value = static_cast<long double>(table[l][m]);
    sum += value * value;
  }
  return std::abs(static_cast<double>(
      sum * 4 * std::numbers::pi_v<long double> / (2 * l + 1) - 1));
}

// The largest difference between a double table and the long double one of
// the same shape, over every stored entry, in d^l_{nm} itself rather than in
// the orthonormalised value. Long double's own ceiling is far above every
// degree asked of it here, so it is the truth to double's rounding.
double DifferenceFromLongDouble(Int lMax, Int mMax, Int n, double theta) {
  const auto inDouble =
      Wigner<double, All, Single, Single>(lMax, mMax, n, theta);
  const auto inLong = Wigner<long double, All, Single, Single>(
      lMax, mMax, n, static_cast<long double>(theta));
  auto worst = 0.0;
  for (auto l = std::abs(n); l <= lMax; l++) {
    const auto norm = std::sqrt((2.0 * l + 1) / (4 * std::numbers::pi));
    const auto mBound = std::min(l, mMax);
    for (auto m = -mBound; m <= mBound; m++) {
      const auto difference =
          std::abs(static_cast<long double>(inDouble[l][m]) - inLong[l][m]);
      worst = std::max(worst, static_cast<double>(difference) / norm);
    }
  }
  return worst;
}

constexpr auto unbounded = std::numeric_limits<Int>::max();

}  // namespace

TEST(LargeDegree, TheSafeDegreeDependsOnTheShape) {
  // Every order, at upper index zero.
  EXPECT_EQ(MaxSafeDegree<double>(), 1822);
  EXPECT_EQ(MaxSafeDegree<long double>(), 30740);
  // Recursed in double, so limited as double is.
  EXPECT_EQ(MaxSafeDegree<float>(), 1822);
  EXPECT_EQ(SafeDegree<double>(unbounded, 0, 0), MaxSafeDegree<double>());

  // Small upper indices lower it a little, for the grids' usual ranges.
  EXPECT_EQ(SafeDegree<double>(unbounded, 0, 2), 1818);
  EXPECT_EQ(SafeDegree<double>(unbounded, 0, 4), 1814);
  EXPECT_EQ(SafeDegree<double>(unbounded, 2, 2), 1818);

  // Restricted orders raise it steeply, and without bound once the seeds
  // cannot underflow where they matter at all.
  EXPECT_EQ(SafeDegree<double>(500, 0, 2), 1903);
  EXPECT_EQ(SafeDegree<double>(200, 0, 0), 5721);
  EXPECT_EQ(SafeDegree<double>(100, 0, 2), 72314);
  EXPECT_EQ(SafeDegree<double>(2, 0, 2), unbounded);
  EXPECT_EQ(SafeDegree<double>(0, 0, 0), unbounded);
  for (auto mMax = Int{100}; mMax < 2000; mMax += 100) {
    EXPECT_GE(SafeDegree<double>(mMax, 0, 2),
              SafeDegree<double>(mMax + 100, 0, 2));
  }

  // Large upper indices lower it, to the worst shape of all at
  // |n| = |m| = 336, which is the closed form every table is safe below.
  EXPECT_EQ(SafeDegreeForAnyShape<double>(), 913);
  EXPECT_EQ(SafeDegreeForAnyShape<long double>(), 15373);
  EXPECT_EQ(SafeDegree<double>(unbounded, 336, 336), 913);
  EXPECT_EQ(SafeDegree<double>(unbounded, 0, unbounded), 913);
  EXPECT_EQ(SafeDegree<double>(unbounded, 1300, 1300), 1683);

  EXPECT_THROW(SafeDegree<double>(-1, 0, 0), std::invalid_argument);
  EXPECT_THROW(SafeDegree<double>(0, 2, 1), std::invalid_argument);
}

TEST(LargeDegree, ATableAboveItsSafeDegreeIsRefused) {
  // Every order at one colatitude, so that the tables which are built cost
  // little: a table of order zero alone is safe at any degree.
  using Table = Wigner<double, All, Single, Single>;
  const auto limit = MaxSafeDegree<double>();
  EXPECT_NO_THROW(Table(limit, limit, 0, 1.0));
  EXPECT_THROW(Table(limit + 1, limit + 1, 0, 1.0), std::invalid_argument);
  EXPECT_NO_THROW(Table(4 * limit, 0, 0, 1.0));

  using Matrices = WignerMatrices<double, All, Single>;
  const auto theta = std::vector<double>{1.0};
  EXPECT_NO_THROW(Matrices(limit, limit, 0, theta));
  EXPECT_THROW(Matrices(limit + 1, limit + 1, 0, theta), std::invalid_argument);

  // A grid refuses before it builds anything, which matters: a generating
  // grid has no table to do the refusing.
  using Grid = GaussLegendreGrid<double, All, All>;
  EXPECT_THROW(Grid(limit + 1, 0, FFTWpp::Estimate), std::invalid_argument);
  EXPECT_THROW(Grid(limit + 1, 0, FFTWpp::Estimate, Chunking::Automatic(),
                    WignerValues::Generated()),
               std::invalid_argument);
}

TEST(LargeDegree, RestrictedOrdersAreRightFarAboveTheFullCeiling) {
  // Orders up to 200 are safe to 5721 in double, three times the ceiling of
  // a table of every order, and right there: the colatitudes are around the
  // ones at which the columns near |m| = 200 first underflow.
  constexpr auto mMax = Int{200};
  const auto lMax = SafeDegree<double>(mMax, 0, 0);
  ASSERT_EQ(lMax, 5721);
  const auto angles = Thorough()
                          ? std::vector<double>{0.020, 0.024, 0.028, 0.035}
                          : std::vector<double>{0.024};
  for (auto theta : angles) {
    EXPECT_LT(DifferenceFromLongDouble(lMax, mMax, 0, theta), 1e-10)
        << "theta = " << theta;
  }
  using Table = Wigner<double, All, Single, Single>;
  EXPECT_THROW(Table(lMax + 1, mMax, 0, 0.024), std::invalid_argument);
}

TEST(LargeDegree, ALargeUpperIndexLowersTheCeiling) {
  // At |n| = 336 the column of order -336 is seeded by sin^672(theta/2) and
  // is wrong from degree 929 in the tenth place, and by 8e-2 at 1003 -- both
  // under the 1827 that once let it through. Its own ceiling is 913, and the
  // table is right there.
  constexpr auto n = Int{336};
  const auto lMax = SafeDegree<double>(unbounded, n, n);
  ASSERT_EQ(lMax, 913);
  const auto angles = Thorough() ? std::vector<double>{0.55, 0.62, 0.67, 0.75}
                                 : std::vector<double>{0.67};
  for (auto theta : angles) {
    EXPECT_LT(DifferenceFromLongDouble(lMax, lMax, n, theta), 1e-11)
        << "theta = " << theta;
  }
  using Table = Wigner<double, All, Single, Single>;
  EXPECT_THROW(Table(lMax + 1, lMax + 1, n, 0.67), std::invalid_argument);
}

TEST(LargeDegree, AnUnderflowingSeedFallsToZeroAndDoesNotStick) {
  // At |n| = 1300 the seed row holds cos^2600(theta/2), e^-1600 at
  // theta = 2. Formed by repeated multiplication in double it once stuck at
  // the smallest denormal -- each factor, above one half, rounded it back --
  // and the columns seeded from it came back as numbers near 1e114 at degree
  // 1650, inside this shape's ceiling of 1683. Carried with an exponent of its
  // own, the seed rounds to zero once, and the table is right.
  constexpr auto n = Int{1300};
  constexpr auto lMax = Int{1650};
  ASSERT_LE(lMax, SafeDegree<double>(unbounded, n, n));
  const auto angles = Thorough() ? std::vector<double>{1.2, 1.5, 1.7, 2.0, 2.4}
                                 : std::vector<double>{1.7, 2.0};
  for (auto theta : angles) {
    EXPECT_LT(DifferenceFromLongDouble(lMax, lMax, n, theta), 1e-11)
        << "theta = " << theta;
  }
}

TEST(LargeDegree, TheTableIsStillUnitaryAtTheSafeDegree) {
  // At about the worst colatitude, asin(1/e), and at each shape's own limit.
  // A hundred degrees above it this is wrong in the ninth place, and by
  // lMax = 2000 in the third.
  using Table = Wigner<double, All, Single, Single>;
  const auto worstAngle = std::asin(std::exp(-1.0));
  const auto angles =
      Thorough() ? std::vector<double>{worstAngle, 0.30, 0.45, 0.9, 1.5707}
                 : std::vector<double>{worstAngle};
  const auto indices =
      Thorough() ? std::vector<Int>{0, 1, 2, 3, 4} : std::vector<Int>{0, 2};
  for (auto theta : angles) {
    for (auto n : indices) {
      const auto lMax = SafeDegree<double>(unbounded, n, n);
      const auto table = Table(lMax, lMax, n, theta);
      for (auto l : {lMax, lMax - lMax / 7}) {
        EXPECT_LT(UnitarityDefect(table, l), 1e-11)
            << "n = " << n << ", theta = " << theta << ", l = " << l;
      }
    }
  }
}

TEST(LargeDegree, ALargeUpperIndexHasAFiniteSeedRow) {
  // The seed row at l = |n| carries sqrt(C(2|n|, .)). Formed before its root
  // is taken, the binomial overflows beyond |n| = 510 or so, and two thirds
  // of this table would come back non-finite.
  constexpr auto lMax = Int{700};
  constexpr auto n = Int{600};
  const auto table = Wigner<double, All, Single, Single>(lMax, lMax, n, 1.0);
  for (auto l = n; l <= lMax; l++) {
    for (auto m = -l; m <= l; m++) {
      ASSERT_TRUE(std::isfinite(table[l][m])) << "l = " << l << ", m = " << m;
    }
  }
  EXPECT_LT(UnitarityDefect(table, n), 1e-11);
  EXPECT_LT(UnitarityDefect(table, lMax), 1e-11);

  // And at the pole, where d^l_{nm}(0) is the identity: the one entry at
  // m = n, and zero everywhere else. The logarithms cannot say this -- the
  // logarithm of the vanishing half-angle sine is not stored as minus
  // infinity -- so the seed row has to ask whether it is at a pole.
  const auto atPole = Wigner<double, All, Single, Single>(lMax, lMax, n, 0.0);
  for (auto l : {n, lMax}) {
    for (auto m = -l; m <= l; m++) {
      const auto wanted =
          m == n ? std::sqrt((2 * l + 1) / (4 * std::numbers::pi)) : 0.0;
      // The seed here comes from logarithms and a hundred degrees of
      // recursion follow it, so this is rounding and not exactness.
      ASSERT_NEAR(atPole[l][m], wanted, 1e-10) << "l = " << l << ", m = " << m;
    }
  }
}

TEST(LargeDegree, ASinglePrecisionTableIsTheDoubleOneRounded) {
  constexpr auto lMax = Int{48};
  constexpr auto nMax = Int{2};
  auto single = std::vector<float>{};
  auto widened = std::vector<double>{};
  for (auto i = 0; i < 9; i++) {
    single.push_back(0.2f + 0.31f * static_cast<float>(i));
    widened.push_back(static_cast<double>(single.back()));
  }

  const auto inSingle =
      Wigner<float, All, All, Multiple>(lMax, lMax, nMax, single);
  const auto inDouble =
      Wigner<double, All, All, Multiple>(lMax, lMax, nMax, widened);
  const auto matrices =
      WignerMatrices<float, All, All>(lMax, lMax, nMax, single);

  for (auto n = -nMax; n <= nMax; n++) {
    for (auto iTheta = Int{0}; iTheta < 9; iTheta++) {
      const auto a = inSingle[n, iTheta];
      const auto b = inDouble[n, iTheta];
      for (auto l = std::abs(n); l <= lMax; l++) {
        for (auto m = -l; m <= l; m++) {
          ASSERT_EQ((a[l, m]), static_cast<float>((b[l, m])))
              << "n = " << n << ", l = " << l << ", m = " << m;
          // And the transform-major table holds the same numbers.
          ASSERT_EQ(
              (matrices[n, m])[static_cast<std::size_t>(
                  (l - std::max<Int>(std::abs(n), std::abs(m))) * 9 + iTheta)],
              (a[l, m]))
              << "n = " << n << ", l = " << l << ", m = " << m;
        }
      }
    }
  }
}

TEST(LargeDegree, ASinglePrecisionGridWorksBeyondItsOwnUnderflow) {
  // A recursion *in* single precision fails from degree 237, and not subtly:
  // whole columns of the table come back zero. Measured on a round trip, the
  // error is 9e-4 at lMax = 256, where few columns have yet gone, 13 at 288,
  // and 3e11 at 320. Recursed in double it is 8e-4 at 288, which is single
  // precision's own rounding. Every kernel is run, since each reaches the
  // recursion by its own route, and the stored and generated routes must
  // agree to the bit.
  using Real = float;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr auto lMax = Int{288};
  constexpr auto n = Int{0};

  auto stored = Grid(lMax, n, FFTWpp::Estimate);
  auto generated = Grid(lMax, n, FFTWpp::Estimate, Chunking::Automatic(),
                        WignerValues::Generated());

  const auto size = static_cast<std::size_t>(stored.CoefficientSize(lMax, n));
  auto given = std::vector<Complex>(size);
  for (std::size_t j = 0; j < size; j++) {
    given[j] =
        Complex{static_cast<Real>(std::sin(0.37 * static_cast<double>(j))),
                static_cast<Real>(std::cos(0.11 * static_cast<double>(j)))};
  }

  const auto roundTrip = [&](const Grid& grid) {
    auto field = std::vector<Complex>(grid.FieldSize());
    auto back = std::vector<Complex>(size);
    grid.InverseTransformation(lMax, n, given, field);
    grid.ForwardTransformation(lMax, n, field, back);
    return std::pair(std::move(field), std::move(back));
  };

  const auto [storedField, storedBack] = roundTrip(stored);
  auto worst = 0.0;
  for (std::size_t j = 0; j < size; j++) {
    worst = std::max<double>(worst, std::abs(storedBack[j] - given[j]));
  }
  EXPECT_LT(worst, 5e-3) << "stored, loop kernel";

  const auto [generatedField, generatedBack] = roundTrip(generated);
  for (std::size_t j = 0; j < storedField.size(); j++) {
    ASSERT_EQ(storedField[j], generatedField[j]) << "sample " << j;
  }
  for (std::size_t j = 0; j < size; j++) {
    ASSERT_EQ(storedBack[j], generatedBack[j]) << "coefficient " << j;
  }

#ifdef GSHTRANS_HAVE_BLAS
  auto matrix = Grid(lMax, n, FFTWpp::Estimate, Chunking::Automatic(),
                     WignerValues::Stored(), TransformKernel::Matrix());
  const auto [matrixField, matrixBack] = roundTrip(matrix);
  auto worstMatrix = 0.0;
  for (std::size_t j = 0; j < size; j++) {
    worstMatrix =
        std::max<double>(worstMatrix, std::abs(matrixBack[j] - given[j]));
  }
  EXPECT_LT(worstMatrix, 5e-3) << "matrix kernel";
#endif
}
