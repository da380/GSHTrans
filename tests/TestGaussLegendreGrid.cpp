#include <gtest/gtest.h>

#include <cstdint>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

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

// The colatitude loop accumulates, so out had to arrive zeroed -- an unstated,
// unchecked precondition that every caller met by accident. Transforming twice
// into one buffer doubled the answer (core-plan.md F1).
TEST(GaussLegendreGrid, ForwardTransformOwnsItsOutputBuffer) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr std::ptrdiff_t lMax = 6;
  constexpr auto tolerance = 1.0e-13;

  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto field = FFTWpp::vector<Complex>(grid.FieldSize());
  auto i = std::ptrdiff_t{0};
  for (auto [theta, phi] : grid.Points()) {
    field[i++] = Complex{std::cos(theta) + 0.5 * std::sin(3.0 * phi),
                         0.25 * std::sin(theta) * std::cos(phi)};
  }

  auto once = FFTWpp::vector<Complex>(grid.CoefficientSize(lMax, 1));
  auto twice = FFTWpp::vector<Complex>(grid.CoefficientSize(lMax, 1));
  std::ranges::fill(once, Complex{});
  std::ranges::fill(twice, Complex{});

  grid.ForwardTransformation(lMax, 1, field, once);
  grid.ForwardTransformation(lMax, 1, field, twice);
  grid.ForwardTransformation(lMax, 1, field, twice);

  auto largest = 0.0;
  for (auto j = std::size_t{0}; j < once.size(); ++j) {
    EXPECT_NEAR(std::abs(twice[j] - once[j]), 0.0, tolerance) << "j = " << j;
    largest = std::max(largest, std::abs(once[j]));
  }
  EXPECT_GT(largest, tolerance);  // and the spectrum is not simply empty

  // Also true of a buffer arriving with rubbish in it rather than a previous
  // answer, which is the case an accumulating routine cannot serve at all.
  auto dirty = FFTWpp::vector<Complex>(once.size());
  std::ranges::fill(dirty, Complex{1.0e3, -1.0e3});
  grid.ForwardTransformation(lMax, 1, field, dirty);
  for (auto j = std::size_t{0}; j < once.size(); ++j) {
    EXPECT_NEAR(std::abs(dirty[j] - once[j]), 0.0, tolerance) << "j = " << j;
  }
}

// Size mismatches were assert-only, so a short output range was a silent heap
// overflow under NDEBUG (core-plan.md F5). This test is meaningful only
// because the suite is run in Release as well as Debug.
TEST(GaussLegendreGrid, RejectsMismatchedRangeSizes) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr std::ptrdiff_t lMax = 4;

  auto grid = Grid(lMax, 1, FFTWpp::Estimate);
  const auto fieldSize = grid.FieldSize();
  const auto coefficientSize = grid.CoefficientSize(lMax, 0);

  auto field = FFTWpp::vector<Complex>(fieldSize);
  auto shortField = FFTWpp::vector<Complex>(fieldSize - 1);
  auto coefficients = FFTWpp::vector<Complex>(coefficientSize);
  auto shortCoefficients = FFTWpp::vector<Complex>(coefficientSize - 1);
  auto longCoefficients = FFTWpp::vector<Complex>(coefficientSize + 1);

  EXPECT_THROW(grid.ForwardTransformation(lMax, 0, shortField, coefficients),
               std::invalid_argument);
  EXPECT_THROW(grid.ForwardTransformation(lMax, 0, field, shortCoefficients),
               std::invalid_argument);
  EXPECT_THROW(grid.ForwardTransformation(lMax, 0, field, longCoefficients),
               std::invalid_argument);
  EXPECT_THROW(grid.InverseTransformation(lMax, 0, shortCoefficients, field),
               std::invalid_argument);
  EXPECT_THROW(grid.InverseTransformation(lMax, 0, coefficients, shortField),
               std::invalid_argument);

  // The message says which range and what was expected.
  try {
    grid.ForwardTransformation(lMax, 0, field, shortCoefficients);
    FAIL() << "expected a throw";
  } catch (const std::invalid_argument& error) {
    const auto message = std::string(error.what());
    EXPECT_NE(message.find("coefficient"), std::string::npos) << message;
    EXPECT_NE(message.find(std::to_string(coefficientSize)),
              std::string::npos)
        << message;
  }
}

// The FFT plans are executed on the grid's own aligned buffers, never on
// caller storage (core-plan.md F3), so a caller may hand over any storage
// aligned for its scalar type -- which is what the field layer promises about
// slice targets. Here the ranges are offset views into plain std::vectors,
// chosen so that they are not on a 64-byte boundary: fftw_malloc'd storage is,
// and FFTW's new-array execute is documented as valid only for buffers sharing
// the planning buffers' alignment.
namespace {

// The first offset, in elements, at which the pointer is not aligned to
// `boundary` bytes.
template <typename T>
std::size_t MisalignedOffset(const T* data, std::size_t boundary) {
  for (auto offset = std::size_t{0}; offset <= boundary / sizeof(T); ++offset) {
    if (reinterpret_cast<std::uintptr_t>(data + offset) % boundary != 0) {
      return offset;
    }
  }
  return 0;
}

// The first offset landing on an 8-byte but not 16-byte boundary. FFTW sorts
// pointers into alignment classes and refuses to say anything about executing
// a plan on storage in a different class from the buffer it was planned on;
// on this build fftw_alignment_of gives 0 for 16-byte-aligned storage and 8
// for this one, so this is the offset that makes new-array execute genuinely
// invalid rather than merely unusual.
std::size_t OddlyAlignedOffset(const double* data) {
  for (auto offset = std::size_t{0}; offset < 4; ++offset) {
    if (reinterpret_cast<std::uintptr_t>(data + offset) % 16 == 8) {
      return offset;
    }
  }
  return 0;
}

}  // namespace

TEST(GaussLegendreGrid, AcceptsUnalignedCallerStorage) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr std::ptrdiff_t lMax = 5;
  constexpr std::size_t boundary = 64;
  constexpr auto tolerance = 1.0e-12;

  auto grid = Grid(lMax, 1, FFTWpp::Estimate);
  const auto fieldSize = static_cast<std::size_t>(grid.FieldSize());

  // Complex field, upper index one: the complex-to-complex path.
  {
    const auto size = static_cast<std::size_t>(grid.CoefficientSize(lMax, 1));
    auto fieldStorage = std::vector<Complex>(fieldSize + 16);
    auto givenStorage = std::vector<Complex>(size + 16);
    auto backStorage = std::vector<Complex>(size + 16);
    auto field = std::span(fieldStorage)
                     .subspan(MisalignedOffset(fieldStorage.data(), boundary),
                              fieldSize);
    auto given = std::span(givenStorage)
                     .subspan(MisalignedOffset(givenStorage.data(), boundary),
                              size);
    auto back = std::span(backStorage)
                    .subspan(MisalignedOffset(backStorage.data(), boundary),
                             size);

    ASSERT_NE(reinterpret_cast<std::uintptr_t>(field.data()) % boundary, 0u);
    ASSERT_NE(reinterpret_cast<std::uintptr_t>(given.data()) % boundary, 0u);

    std::ranges::fill(given, Complex{});
    const auto indices = GSHIndices<All>(lMax, lMax, 1);
    given[indices.Index(3, -2)] = Complex{0.5, 0.25};
    given[indices.Index(4, 1)] = Complex{-0.75, 0.125};

    grid.InverseTransformation(lMax, 1, given, field);
    grid.ForwardTransformation(lMax, 1, field, back);

    for (auto j = std::size_t{0}; j < size; ++j) {
      EXPECT_NEAR(std::abs(back[j] - given[j]), 0.0, tolerance) << "j = " << j;
    }
    // Storage past the end of each view must be untouched. The offset may
    // legitimately be zero when the allocation is already misaligned, so the
    // tail, not the head, is the padding that is always outside.
    EXPECT_EQ(fieldStorage.back(), Complex{});
    EXPECT_EQ(backStorage.back(), Complex{});
  }

  // Real field at upper index zero: the real-to-complex path, and the one
  // whose elements are small enough to land off a 16-byte boundary too.
  {
    const auto size = static_cast<std::size_t>(grid.RealCoefficientSize(lMax));
    auto fieldStorage = std::vector<Real>(fieldSize + 16);
    auto givenStorage = std::vector<Complex>(size + 16);
    auto backStorage = std::vector<Complex>(size + 16);
    auto field = std::span(fieldStorage)
                     .subspan(OddlyAlignedOffset(fieldStorage.data()),
                              fieldSize);
    auto given = std::span(givenStorage)
                     .subspan(MisalignedOffset(givenStorage.data(), boundary),
                              size);
    auto back = std::span(backStorage)
                    .subspan(MisalignedOffset(backStorage.data(), boundary),
                             size);

    // Not merely off a cache line: in a different FFTW alignment class from
    // the fftw_malloc'd buffers the plans are made on.
    auto* planningStorage = static_cast<Real*>(fftw_malloc(sizeof(Real)));
    ASSERT_NE(fftw_alignment_of(field.data()),
              fftw_alignment_of(planningStorage));
    fftw_free(planningStorage);

    std::ranges::fill(given, Complex{});
    const auto indices = GSHIndices<NonNegative>(lMax, lMax, 0);
    given[indices.Index(2, 0)] = Complex{1.5, 0.0};
    given[indices.Index(4, 3)] = Complex{-0.25, 0.75};

    grid.InverseTransformation(lMax, 0, given, field);
    grid.ForwardTransformation(lMax, 0, field, back);

    for (auto j = std::size_t{0}; j < size; ++j) {
      EXPECT_NEAR(std::abs(back[j] - given[j]), 0.0, tolerance) << "j = " << j;
    }
    EXPECT_EQ(fieldStorage.back(), Real{});
  }
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
