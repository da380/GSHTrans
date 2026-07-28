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

TEST(GaussLegendreGrid, Coeff2CoeffDoubleR2C) {
  using Scalar = double;
  bool result = Coeff2Coeff<Scalar, All, All>();
  EXPECT_FALSE(result);
}

TEST(GaussLegendreGrid, Coeff2CoeffLongDoubleR2C) {
  using Scalar = long double;
  bool result = Coeff2Coeff<Scalar, All, All>();
  EXPECT_FALSE(result);
}

TEST(GaussLegendreGrid, Coeff2CoeffDoubleC2C) {
  using Scalar = std::complex<double>;
  bool result = Coeff2Coeff<Scalar, All, All>();
  EXPECT_FALSE(result);
}

TEST(GaussLegendreGrid, Coeff2CoeffLongDoubleC2C) {
  using Scalar = std::complex<long double>;
  bool result = Coeff2Coeff<Scalar, All, All>();
  EXPECT_FALSE(result);
}
