#include <gtest/gtest.h>

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

// Wigner tables at large degree, and in single precision.
//
// The recursion in degree is seeded at l = |m| with a value of about
// (sin theta)^m, and the column it starts grows back to order one only once
// l sin theta reaches m. A seed that has underflowed *and* still matters needs
// lMax * s |ln s| > |ln min| for some s = sin theta, and s |ln s| is largest
// at s = 1/e -- so the plain recursion is safe if and only if
//
//     lMax  <  e * |ln min|,
//
// and a little below that, e * (|ln min| - |ln eps|), if the seed is not to go
// denormal and lose bits on the way. Measured, the formula puts the onset
// within a per cent in both precisions. Above it the table was silently wrong
// by order one; now it is refused.
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

}  // namespace

TEST(LargeDegree, TheSafeDegreeIsWhatTheFormulaGives) {
  EXPECT_EQ(MaxSafeDegree<double>(), 1827);
  EXPECT_EQ(MaxSafeDegree<long double>(), 30747);
  // Recursed in double, so limited as double is.
  EXPECT_EQ(MaxSafeDegree<float>(), 1827);
}

TEST(LargeDegree, ATableAboveTheSafeDegreeIsRefused) {
  // One order and one colatitude, so that the tables which are built cost
  // nothing: the limit is on the degree whatever else is asked for.
  using Table = Wigner<double, All, Single, Single>;
  const auto limit = MaxSafeDegree<double>();
  EXPECT_NO_THROW(Table(limit, 0, 0, 1.0));
  EXPECT_THROW(Table(limit + 1, 0, 0, 1.0), std::invalid_argument);

  using Matrices = WignerMatrices<double, All, Single>;
  const auto theta = std::vector<double>{1.0};
  EXPECT_NO_THROW(Matrices(limit, 0, 0, theta));
  EXPECT_THROW(Matrices(limit + 1, 0, 0, theta), std::invalid_argument);

  // A grid refuses before it builds anything, which matters: a generating
  // grid has no table to do the refusing.
  using Grid = GaussLegendreGrid<double, All, All>;
  EXPECT_THROW(Grid(limit + 1, 0, FFTWpp::Estimate), std::invalid_argument);
  EXPECT_THROW(Grid(limit + 1, 0, FFTWpp::Estimate, Chunking::Automatic(),
                    WignerValues::Generated()),
               std::invalid_argument);
}

TEST(LargeDegree, TheTableIsStillUnitaryAtTheSafeDegree) {
  // At the worst colatitude, asin(1/e), and at the limit itself. One step
  // above the bare e |ln min| this is wrong in the ninth place, and by
  // lMax = 2000 in the third.
  using Table = Wigner<double, All, Single, Single>;
  const auto lMax = MaxSafeDegree<double>();
  const auto worstAngle = std::asin(std::exp(-1.0));
  const auto angles =
      Thorough() ? std::vector<double>{worstAngle, 0.30, 0.45, 0.9, 1.5707}
                 : std::vector<double>{worstAngle};
  const auto indices =
      Thorough() ? std::vector<Int>{0, 1, 2, 3, 4} : std::vector<Int>{0, 2};
  for (auto theta : angles) {
    for (auto n : indices) {
      const auto table = Table(lMax, lMax, n, theta);
      for (auto l : {lMax, lMax - lMax / 7}) {
        EXPECT_LT(UnitarityDefect(table, l), 1e-11)
            << "n = " << n << ", theta = " << theta << ", l = " << l;
      }
    }
  }
}

TEST(LargeDegree, ALargeUpperIndexHasAFiniteSeedRow) {
  // The seed row at l = |n| carries sqrt(C(2|n|, .)), and the binomial was
  // formed before its root was taken: beyond |n| = 510 or so it overflowed,
  // and two thirds of this table came back non-finite.
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
  // error was 9e-4 at lMax = 256, where few columns have yet gone, 13 at 288,
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

  const auto size = stored.CoefficientSize(lMax, n);
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
