#include <gtest/gtest.h>

#include <GSHTrans/All>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>

namespace {

using namespace GSHTrans;

using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

constexpr std::ptrdiff_t lMax = 5;
constexpr std::ptrdiff_t n = 2;
constexpr Real tolerance = 2.0e-11;

template <std::ptrdiff_t N>
using ReducedExpansion = RealCanonicalComponentExpansion<N, Grid>;

template <std::ptrdiff_t N>
using FullExpansion = ComplexCanonicalComponentExpansion<N, Grid>;

int Phase(std::ptrdiff_t exponent) {
  return exponent % 2 == 0 ? 1 : -1;
}

void ExpectNear(Complex actual, Complex expected) {
  EXPECT_NEAR(actual.real(), expected.real(), tolerance);
  EXPECT_NEAR(actual.imag(), expected.imag(), tolerance);
}

void MakeRealSamples(const Grid& grid, FFTWpp::vector<Real>& realSamples,
                     FFTWpp::vector<Complex>& complexSamples) {
  auto i = std::ptrdiff_t{0};
  for (auto [theta, phi] : grid.Points()) {
    const auto value =
        1.25 + 0.4 * std::cos(theta) +
        0.3 * std::sin(theta) * std::cos(phi) -
        0.2 * std::sin(2.0 * theta) * std::sin(2.0 * phi) +
        0.15 * std::cos(3.0 * theta) * std::cos(3.0 * phi) +
        0.45 * (1.0 + 0.3 * std::cos(theta)) * std::cos(lMax * phi);
    realSamples[i] = value;
    complexSamples[i] = Complex{value, 0.0};
    ++i;
  }
}

template <std::ptrdiff_t N>
void Forward(const Grid& grid, const FFTWpp::vector<Real>& realSamples,
             const FFTWpp::vector<Complex>& complexSamples,
             ReducedExpansion<N>& reduced, FullExpansion<N>& full) {
  auto reducedData = reduced.Data();
  auto fullData = full.Data();
  std::ranges::fill(reducedData, Complex{});
  std::ranges::fill(fullData, Complex{});
  grid.ForwardTransformation(lMax, N, realSamples, reducedData);
  grid.ForwardTransformation(lMax, N, complexSamples, fullData);
}

}  // namespace

TEST(RealFieldSymmetry, ReducedSpectrumMatchesFullPositiveOrders) {
  auto grid = Grid(lMax, n, FFTWpp::Estimate);
  auto realSamples = FFTWpp::vector<Real>(grid.FieldSize());
  auto complexSamples = FFTWpp::vector<Complex>(grid.FieldSize());
  MakeRealSamples(grid, realSamples, complexSamples);

  auto reducedPlus = ReducedExpansion<n>(grid);
  auto fullPlus = FullExpansion<n>(grid);
  auto reducedMinus = ReducedExpansion<-n>(grid);
  auto fullMinus = FullExpansion<-n>(grid);
  Forward(grid, realSamples, complexSamples, reducedPlus, fullPlus);
  Forward(grid, realSamples, complexSamples, reducedMinus, fullMinus);

  for (auto [l, m] : reducedPlus.Indices()) {
    if (l != lMax || m != lMax) {
      ExpectNear(reducedPlus[l, m], fullPlus[l, m]);
      ExpectNear(reducedMinus[l, m], fullMinus[l, m]);
    }
  }

  for (auto l : reducedPlus.Degrees()) {
    const auto plusZonal = reducedPlus[l, 0];
    const auto minusZonal = reducedMinus[l, 0];
    EXPECT_NEAR(plusZonal.imag(), 0.0, tolerance);
    EXPECT_NEAR(minusZonal.imag(), 0.0, tolerance);
  }

  const auto plusNyquist = reducedPlus[lMax, lMax];
  const auto minusNyquist = reducedMinus[lMax, lMax];
  EXPECT_NEAR(plusNyquist.imag(), 0.0, tolerance);
  EXPECT_NEAR(minusNyquist.imag(), 0.0, tolerance);
  EXPECT_GT(std::abs(plusNyquist), tolerance);
  EXPECT_GT(std::abs(minusNyquist), tolerance);

  ExpectNear(fullMinus[lMax, -lMax],
             static_cast<Real>(Phase(lMax - n)) *
                 std::conj(plusNyquist));
  ExpectNear(fullPlus[lMax, -lMax],
             static_cast<Real>(Phase(lMax + n)) *
                 std::conj(minusNyquist));
  ExpectNear(fullPlus[lMax, lMax], Complex{});
  ExpectNear(fullMinus[lMax, lMax], Complex{});
}

TEST(RealFieldSymmetry, FullTransformsSatisfyCrossUpperIndexIdentity) {
  auto grid = Grid(lMax, n, FFTWpp::Estimate);
  auto realSamples = FFTWpp::vector<Real>(grid.FieldSize());
  auto complexSamples = FFTWpp::vector<Complex>(grid.FieldSize());
  MakeRealSamples(grid, realSamples, complexSamples);

  auto reducedPlus = ReducedExpansion<n>(grid);
  auto reducedMinus = ReducedExpansion<-n>(grid);
  auto fullPlus = FullExpansion<n>(grid);
  auto fullMinus = FullExpansion<-n>(grid);
  Forward(grid, realSamples, complexSamples, reducedPlus, fullPlus);
  Forward(grid, realSamples, complexSamples, reducedMinus, fullMinus);

  for (auto [l, m] : fullPlus.Indices()) {
    if (l == lMax && std::abs(m) == lMax) {
      continue;
    }
    const auto expected =
        static_cast<Real>(Phase(m - n)) * std::conj(fullPlus[l, m]);
    ExpectNear(fullMinus[l, -m], expected);
  }
}

TEST(RealFieldSymmetry, ReducedInverseUsesTheImplicitHermitianPair) {
  auto grid = Grid(lMax, n, FFTWpp::Estimate);
  auto realSamples = FFTWpp::vector<Real>(grid.FieldSize());
  auto complexSamples = FFTWpp::vector<Complex>(grid.FieldSize());
  MakeRealSamples(grid, realSamples, complexSamples);

  auto reduced = ReducedExpansion<n>(grid);
  auto unusedFull = FullExpansion<n>(grid);
  Forward(grid, realSamples, complexSamples, reduced, unusedFull);

  auto positive = FullExpansion<n>(grid);
  auto negative = FullExpansion<-n>(grid);
  auto positiveData = positive.Data();
  auto negativeData = negative.Data();
  std::ranges::fill(positiveData, Complex{});
  std::ranges::fill(negativeData, Complex{});

  for (auto [l, m] : reduced.Indices()) {
    positive[l, m] = reduced[l, m];
    if (m > 0 && m < lMax) {
      negative[l, -m] =
          static_cast<Real>(Phase(m - n)) * std::conj(reduced[l, m]);
    }
  }

  auto reducedField = FFTWpp::vector<Real>(grid.FieldSize());
  auto positiveField = FFTWpp::vector<Complex>(grid.FieldSize());
  auto negativeField = FFTWpp::vector<Complex>(grid.FieldSize());
  auto reducedData = reduced.Data();
  grid.InverseTransformation(lMax, n, reducedData, reducedField);
  grid.InverseTransformation(lMax, n, positiveData, positiveField);
  grid.InverseTransformation(lMax, -n, negativeData, negativeField);

  for (auto i = std::ptrdiff_t{0}; i < grid.FieldSize(); ++i) {
    const auto hermitianValue = positiveField[i] + negativeField[i];
    EXPECT_NEAR(hermitianValue.imag(), 0.0, tolerance);
    EXPECT_NEAR(reducedField[i], hermitianValue.real(), tolerance);
  }
}
