#include <gtest/gtest.h>

#include <limits>
#include <numbers>
#include <stdexcept>

#include "CheckCoeff2Coeff.h"

TEST(GaussLegendreGrid, DegreeZeroGeometryAndWeights) {
  using Grid = GaussLegendreGrid<double, All, All>;
  auto grid = Grid(0, 0, FFTWpp::Estimate);

  ASSERT_EQ(grid.NumberOfCoLatitudes(), 1);
  ASSERT_EQ(grid.NumberOfLongitudes(), 1);
  ASSERT_EQ(grid.FieldSize(), 1);

  EXPECT_DOUBLE_EQ(*grid.Longitudes().begin(), 0.0);
  EXPECT_NEAR(*grid.LongitudeWeights().begin(),
              2.0 * std::numbers::pi_v<double>,
              16.0 * std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(*grid.Weights().begin(), 4.0 * std::numbers::pi_v<double>,
              32.0 * std::numeric_limits<double>::epsilon());
}

TEST(GaussLegendreGrid, DegreeZeroRealKnownAnswer) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  auto grid = Grid(0, 0, FFTWpp::Estimate);
  auto field = FFTWpp::vector<Real>(1);
  auto coefficients = FFTWpp::vector<Complex>(1);
  field[0] = 3.25;

  grid.ForwardTransformation(0, 0, field, coefficients);
  const auto expectedCoefficient =
      field[0] * 2.0 / std::numbers::inv_sqrtpi_v<Real>;
  EXPECT_NEAR(coefficients[0].real(), expectedCoefficient,
              32.0 * std::numeric_limits<Real>::epsilon());
  EXPECT_DOUBLE_EQ(coefficients[0].imag(), 0.0);

  auto reconstructed = FFTWpp::vector<Real>(1);
  grid.InverseTransformation(0, 0, coefficients, reconstructed);
  EXPECT_NEAR(reconstructed[0], field[0],
              32.0 * std::numeric_limits<Real>::epsilon());
}

TEST(GaussLegendreGrid, DegreeZeroComplexKnownAnswer) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  auto grid = Grid(0, 0, FFTWpp::Estimate);
  auto field = FFTWpp::vector<Complex>(1);
  auto coefficients = FFTWpp::vector<Complex>(1);
  field[0] = Complex(2.5, -0.75);

  grid.ForwardTransformation(0, 0, field, coefficients);
  const auto expectedCoefficient =
      field[0] * 2.0 / std::numbers::inv_sqrtpi_v<Real>;
  EXPECT_NEAR(coefficients[0].real(), expectedCoefficient.real(),
              32.0 * std::numeric_limits<Real>::epsilon());
  EXPECT_NEAR(coefficients[0].imag(), expectedCoefficient.imag(),
              32.0 * std::numeric_limits<Real>::epsilon());

  auto reconstructed = FFTWpp::vector<Complex>(1);
  grid.InverseTransformation(0, 0, coefficients, reconstructed);
  EXPECT_NEAR(reconstructed[0].real(), field[0].real(),
              32.0 * std::numeric_limits<Real>::epsilon());
  EXPECT_NEAR(reconstructed[0].imag(), field[0].imag(),
              32.0 * std::numeric_limits<Real>::epsilon());
}

TEST(GaussLegendreGrid, DegreeZeroTruncationUsesTheWholeGrid) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr auto tolerance = 64.0 * std::numeric_limits<Real>::epsilon();

  auto grid = Grid(2, 0, FFTWpp::Estimate);
  auto realField = FFTWpp::vector<Real>(grid.FieldSize());
  auto complexField = FFTWpp::vector<Complex>(grid.FieldSize());
  for (auto i = std::ptrdiff_t{0}; i < grid.FieldSize(); ++i) {
    const auto value = static_cast<Real>(i + 1);
    realField[i] = value;
    complexField[i] = Complex{value, -0.5 * value};
  }

  auto expectedReal = Real{};
  auto expectedComplex = Complex{};
  auto i = std::ptrdiff_t{0};
  for (const auto weight : grid.Weights()) {
    expectedReal += realField[i] * weight;
    expectedComplex += complexField[i] * weight;
    ++i;
  }
  const auto y00 = std::numbers::inv_sqrtpi_v<Real> / 2.0;
  expectedReal *= y00;
  expectedComplex *= y00;

  auto realCoefficient = FFTWpp::vector<Complex>(1);
  auto complexCoefficient = FFTWpp::vector<Complex>(1);
  grid.ForwardTransformation(0, 0, realField, realCoefficient);
  grid.ForwardTransformation(0, 0, complexField, complexCoefficient);
  EXPECT_NEAR(realCoefficient[0].real(), expectedReal, tolerance);
  EXPECT_NEAR(realCoefficient[0].imag(), 0.0, tolerance);
  EXPECT_NEAR(complexCoefficient[0].real(), expectedComplex.real(), tolerance);
  EXPECT_NEAR(complexCoefficient[0].imag(), expectedComplex.imag(), tolerance);

  const auto constant = Complex{2.0, -0.75};
  complexCoefficient[0] = constant / y00;
  realCoefficient[0] = Complex{constant.real() / y00, 0.0};
  auto reconstructedReal = FFTWpp::vector<Real>(grid.FieldSize());
  auto reconstructedComplex = FFTWpp::vector<Complex>(grid.FieldSize());
  grid.InverseTransformation(0, 0, realCoefficient, reconstructedReal);
  grid.InverseTransformation(0, 0, complexCoefficient, reconstructedComplex);
  for (auto j = std::ptrdiff_t{0}; j < grid.FieldSize(); ++j) {
    EXPECT_NEAR(reconstructedReal[j], constant.real(), tolerance);
    EXPECT_NEAR(reconstructedComplex[j].real(), constant.real(), tolerance);
    EXPECT_NEAR(reconstructedComplex[j].imag(), constant.imag(), tolerance);
  }
}

TEST(GaussLegendreGrid, RejectsUnsupportedTransformRequests) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  auto onePointGrid = Grid(0, 0, FFTWpp::Estimate);
  auto onePointField = FFTWpp::vector<Real>(onePointGrid.FieldSize());
  auto degreeOneCoefficients = FFTWpp::vector<Complex>(3);

  EXPECT_THROW(onePointGrid.ForwardTransformation(
                   1, 0, onePointField, degreeOneCoefficients),
               std::invalid_argument);
  EXPECT_THROW(onePointGrid.InverseTransformation(
                   1, 0, degreeOneCoefficients, onePointField),
               std::invalid_argument);

  auto largerGrid = Grid(2, 2, FFTWpp::Estimate);
  auto largerField = FFTWpp::vector<Real>(largerGrid.FieldSize());
  auto degreeThreeCoefficients = FFTWpp::vector<Complex>(10);
  EXPECT_THROW(largerGrid.ForwardTransformation(
                   3, 0, largerField, degreeThreeCoefficients),
               std::invalid_argument);
  EXPECT_THROW(largerGrid.InverseTransformation(
                   3, 0, degreeThreeCoefficients, largerField),
               std::invalid_argument);
  EXPECT_THROW(largerGrid.ForwardTransformation(
                   1, 2, largerField, degreeOneCoefficients),
               std::invalid_argument);

  auto scalarGrid = Grid(2, 0, FFTWpp::Estimate);
  EXPECT_THROW(scalarGrid.ForwardTransformation(
                   1, 1, largerField, degreeOneCoefficients),
               std::invalid_argument);
}

// nPhi must exceed 2 * lMax so that the orders m = +-lMax are separate
// discrete modes, and should be a length FFTW transforms quickly
// (core-plan.md F2, step D).
static_assert(FastFFTSize(1) == 1);
static_assert(FastFFTSize(2) == 2);
static_assert(FastFFTSize(11) == 11);
static_assert(!IsFastFFTSize(121));  // 11^2: eleven may appear only once
static_assert(!IsFastFFTSize(143));  // 11 * 13: and not alongside thirteen
static_assert(FastFFTSize(143) == 144);
// The case that motivates the rule: 2 * 256 + 2 = 514 = 2 * 257, and 257 is
// prime. The next fast length is 520 = 2^3 * 5 * 13.
static_assert(!IsFastFFTSize(514));
static_assert(FastFFTSize(513) == 520);

TEST(GaussLegendreGrid, LongitudeCountResolvesTheHighestOrders) {
  using Grid = GaussLegendreGrid<double, All, All>;

  for (auto lMax : {0, 1, 2, 3, 5, 8, 16, 33}) {
    auto grid = Grid(lMax, 0, FFTWpp::Estimate);
    const auto nPhi = static_cast<std::ptrdiff_t>(grid.NumberOfLongitudes());
    EXPECT_GE(nPhi, 2 * lMax + 1) << "lMax = " << lMax;
    EXPECT_TRUE(IsFastFFTSize(nPhi)) << "lMax = " << lMax << ", nPhi = " << nPhi;
    EXPECT_EQ(grid.FieldSize(), (lMax + 1) * nPhi) << "lMax = " << lMax;
  }
}

// The round trip that the old sizing made impossible: a field carrying only
// the order m = +lMax must come back with that coefficient intact and its
// m = -lMax partner still zero. At nPhi = 2 * lMax the two were one mode, and
// the forward transform zeroed the (lMax, lMax) coefficient outright.
TEST(GaussLegendreGrid, HighestOrdersSurviveARoundTrip) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr std::ptrdiff_t lMax = 6;
  constexpr auto tolerance = 1.0e-12;

  auto grid = Grid(lMax, 0, FFTWpp::Estimate);
  const auto indices = GSHIndices<All>(lMax, lMax, 0);

  for (auto m : {lMax, -lMax}) {
    auto coefficients = FFTWpp::vector<Complex>(indices.Size());
    std::ranges::fill(coefficients, Complex{});
    const auto value = Complex{0.75, -0.25};
    coefficients[indices.Index(lMax, m)] = value;

    auto field = FFTWpp::vector<Complex>(grid.FieldSize());
    auto recovered = FFTWpp::vector<Complex>(indices.Size());
    std::ranges::fill(recovered, Complex{});
    grid.InverseTransformation(lMax, 0, coefficients, field);
    grid.ForwardTransformation(lMax, 0, field, recovered);

    EXPECT_NEAR(recovered[indices.Index(lMax, m)].real(), value.real(),
                tolerance)
        << "m = " << m;
    EXPECT_NEAR(recovered[indices.Index(lMax, m)].imag(), value.imag(),
                tolerance)
        << "m = " << m;
    EXPECT_NEAR(std::abs(recovered[indices.Index(lMax, -m)]), 0.0, tolerance)
        << "m = " << m << " leaked into its partner";
  }
}

TEST(GaussLegendreGrid, ForBandGivesRequestedHeadroom) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr std::ptrdiff_t band = 8;

  auto exact = Grid::ForBand(band, 0, 1.0, FFTWpp::Estimate);
  auto threeHalves = Grid::ForBand(band, 0, 1.5, FFTWpp::Estimate);
  auto doubled = Grid::ForBand(band, 0, 2.0, FFTWpp::Estimate);

  EXPECT_EQ(exact.MaxDegree(), band);
  EXPECT_EQ(threeHalves.MaxDegree(), 12);
  EXPECT_EQ(doubled.MaxDegree(), 2 * band);

  // A non-integer product rounds up, never down: the point of the parameter is
  // headroom, and rounding down would silently remove it.
  EXPECT_EQ(Grid::ForBand(band, 0, 1.1, FFTWpp::Estimate).MaxDegree(), 9);

  EXPECT_THROW(Grid::ForBand(band, 0, 0.5, FFTWpp::Estimate),
               std::invalid_argument);
  EXPECT_THROW(Grid::ForBand(-1, 0, 1.0, FFTWpp::Estimate),
               std::invalid_argument);

  // An oversampled grid still transforms at the band, which is the point:
  // headroom in the quadrature, truncation in the transform.
  auto field = FFTWpp::vector<Complex>(doubled.FieldSize());
  auto coefficients =
      FFTWpp::vector<Complex>(doubled.CoefficientSize(band, 0));
  EXPECT_NO_THROW(doubled.ForwardTransformation(band, 0, field, coefficients));
}

// The round-trip tests draw their grid, degree, upper index and coefficients
// from one seeded generator, and report the seed so that a failure can be
// reproduced with GSHTRANS_TEST_SEED (core-plan.md F12).
TEST(GaussLegendreGrid, Coeff2CoeffDoubleR2C) {
  const auto seed = GSHTransTest::TestSeed();
  auto gen = GSHTransTest::MakeGenerator(seed);
  EXPECT_FALSE((Coeff2Coeff<double, All, All>(gen))) << "seed = " << seed;
}

TEST(GaussLegendreGrid, Coeff2CoeffLongDoubleR2C) {
  const auto seed = GSHTransTest::TestSeed();
  auto gen = GSHTransTest::MakeGenerator(seed);
  EXPECT_FALSE((Coeff2Coeff<long double, All, All>(gen))) << "seed = " << seed;
}

TEST(GaussLegendreGrid, Coeff2CoeffDoubleC2C) {
  const auto seed = GSHTransTest::TestSeed();
  auto gen = GSHTransTest::MakeGenerator(seed);
  EXPECT_FALSE((Coeff2Coeff<std::complex<double>, All, All>(gen)))
      << "seed = " << seed;
}

TEST(GaussLegendreGrid, Coeff2CoeffLongDoubleC2C) {
  const auto seed = GSHTransTest::TestSeed();
  auto gen = GSHTransTest::MakeGenerator(seed);
  EXPECT_FALSE((Coeff2Coeff<std::complex<long double>, All, All>(gen)))
      << "seed = " << seed;
}
