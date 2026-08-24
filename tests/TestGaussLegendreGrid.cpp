#include <gtest/gtest.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
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

// The grid is a value-semantic handle over shared immutable state
// (core-plan.md step B). Copying one must not copy the Wigner table, which at
// production sizes is hundreds of megabytes against a couple for a field.
TEST(GaussLegendreGrid, IsAValueSemanticHandle) {
  using Grid = GaussLegendreGrid<double, All, All>;

  static_assert(!std::is_default_constructible_v<Grid>,
                "a grid without a quadrature is not a grid (F4)");
  static_assert(std::copy_constructible<Grid>);
  static_assert(std::is_copy_assignable_v<Grid>);
  static_assert(std::is_nothrow_move_constructible_v<Grid>);
  // Small and trivially cheap to copy, which is the property F9 is about:
  // copying a grid used to copy hundreds of megabytes of Wigner table.
  //
  // It used to be exactly a shared_ptr. [C18] added the chunking policy and
  // the planner flag beside it -- read per call, deciding nothing about the
  // table -- so that changing either is a pointer copy rather than a table
  // rebuild. That grows the handle by a few words and changes nothing about
  // what the assertion is for, so the bound is stated as a bound.
  static_assert(sizeof(Grid) <= 4 * sizeof(std::shared_ptr<void>),
                "a grid should be a handle and a few scalars, never a table");

  auto grid = Grid(8, 2, FFTWpp::Estimate);
  auto copy = grid;
  auto assigned = Grid(4, 0, FFTWpp::Estimate);
  assigned = grid;

  EXPECT_EQ(grid.Identity(), copy.Identity());
  EXPECT_EQ(grid.Identity(), assigned.Identity());

  // Equal parameters are not the same grid: identity is the handle, not the
  // configuration.
  auto other = Grid(8, 2, FFTWpp::Estimate);
  EXPECT_NE(grid.Identity(), other.Identity());
  EXPECT_EQ(grid.MaxDegree(), other.MaxDegree());
  EXPECT_EQ(grid.FieldSize(), other.FieldSize());
}

// A copy keeps the implementation alive after the grid it was copied from has
// gone. Without shared ownership this is a use-after-free, and it is exactly
// what the field layer does when a terminal holds a grid by value.
TEST(GaussLegendreGrid, CopiesOutliveTheGridTheyCameFrom) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr std::ptrdiff_t lMax = 6;
  constexpr auto tolerance = 1.0e-12;

  auto escaped = [] {
    auto local = Grid(lMax, 1, FFTWpp::Estimate);
    auto copy = local;
    return copy;
  }();

  auto given = FFTWpp::vector<Complex>(escaped.CoefficientSize(lMax, 1));
  auto field = FFTWpp::vector<Complex>(escaped.FieldSize());
  auto back = FFTWpp::vector<Complex>(given.size());
  std::ranges::fill(given, Complex{});
  const auto indices = GSHIndices<All>(lMax, lMax, 1);
  given[indices.Index(4, -3)] = Complex{0.5, -0.25};

  escaped.InverseTransformation(lMax, 1, given, field);
  escaped.ForwardTransformation(lMax, 1, field, back);
  for (auto j = std::size_t{0}; j < given.size(); ++j) {
    EXPECT_NEAR(std::abs(back[j] - given[j]), 0.0, tolerance) << "j = " << j;
  }
}

// Everything a transform needs is const, so a shared grid can be used through
// a const handle from several places at once.
TEST(GaussLegendreGrid, IsUsableThroughAConstHandle) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr std::ptrdiff_t lMax = 4;

  const auto grid = Grid(lMax, 1, FFTWpp::Estimate);
  auto field = FFTWpp::vector<Complex>(grid.FieldSize());
  auto coefficients = FFTWpp::vector<Complex>(grid.CoefficientSize(lMax, 0));

  EXPECT_NO_THROW(grid.ForwardTransformation(lMax, 0, field, coefficients));
  EXPECT_NO_THROW(grid.InverseTransformation(lMax, 0, coefficients, field));
  EXPECT_EQ(std::ranges::distance(grid.Points()), grid.FieldSize());
  EXPECT_EQ(std::ranges::distance(grid.ProjectFunction(
                [](auto theta, auto phi) { return theta + phi; })),
            grid.FieldSize());
}

// Plans and work buffers are made once per thread per shape and kept
// (core-plan.md step E). Two things that must stay true: a grid shared between
// threads still gives every thread the right answer, and a thread that uses
// several grids gets the right buffers for each.
TEST(GaussLegendreGrid, OneGridServesManyThreads) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr std::ptrdiff_t lMax = 8;
  constexpr std::ptrdiff_t n = 1;
  constexpr auto tolerance = 1.0e-12;

  auto grid = Grid(lMax, n, FFTWpp::Estimate);
  const auto indices = GSHIndices<All>(lMax, lMax, n);

  auto Given = [&](int seed) {
    auto given = FFTWpp::vector<Complex>(indices.Size());
    for (auto j = std::size_t{0}; j < given.size(); ++j) {
      given[j] = Complex{0.1 * seed + 0.01 * j, -0.05 * seed + 0.02 * j};
    }
    return given;
  };

  constexpr auto threadCount = 8;

  // What each thread should produce, computed serially first.
  auto reference = std::vector<FFTWpp::vector<Complex>>{};
  for (auto t = 0; t < threadCount; ++t) {
    auto given = Given(t);
    auto field = FFTWpp::vector<Complex>(grid.FieldSize());
    auto back = FFTWpp::vector<Complex>(indices.Size());
    grid.InverseTransformation(lMax, n, given, field);
    grid.ForwardTransformation(lMax, n, field, back);
    reference.push_back(back);
  }

  auto results = std::vector<FFTWpp::vector<Complex>>(
      threadCount, FFTWpp::vector<Complex>(indices.Size()));
  auto threads = std::vector<std::thread>{};
  for (auto t = 0; t < threadCount; ++t) {
    threads.emplace_back([&, t] {
      auto given = Given(t);
      auto field = FFTWpp::vector<Complex>(grid.FieldSize());
      // Several round trips, so the workspace is reused and not merely made.
      for (auto repeat = 0; repeat < 4; ++repeat) {
        grid.InverseTransformation(lMax, n, given, field);
        grid.ForwardTransformation(lMax, n, field, results[t]);
      }
    });
  }
  for (auto& thread : threads) thread.join();

  for (auto t = 0; t < threadCount; ++t) {
    for (auto j = std::size_t{0}; j < reference[t].size(); ++j) {
      EXPECT_NEAR(std::abs(results[t][j] - reference[t][j]), 0.0, tolerance)
          << "thread " << t << ", coefficient " << j;
    }
  }
}

TEST(GaussLegendreGrid, OneThreadServesManyGrids) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr auto tolerance = 1.0e-12;

  // Different degrees mean different nPhi, hence different cached shapes.
  for (auto repeat = 0; repeat < 3; ++repeat) {
    for (auto lMax :
         {std::ptrdiff_t{4}, std::ptrdiff_t{7}, std::ptrdiff_t{12}}) {
      auto grid = Grid(lMax, 0, FFTWpp::Estimate);
      const auto indices = GSHIndices<All>(lMax, lMax, 0);
      auto given = FFTWpp::vector<Complex>(indices.Size());
      for (auto j = std::size_t{0}; j < given.size(); ++j) {
        given[j] = Complex{0.25 + 0.01 * j, -0.125 + 0.02 * j};
      }
      auto field = FFTWpp::vector<Complex>(grid.FieldSize());
      auto back = FFTWpp::vector<Complex>(indices.Size());
      grid.InverseTransformation(lMax, 0, given, field);
      grid.ForwardTransformation(lMax, 0, field, back);
      for (auto j = std::size_t{0}; j < given.size(); ++j) {
        EXPECT_NEAR(std::abs(back[j] - given[j]), 0.0, tolerance)
            << "lMax " << lMax << ", coefficient " << j;
      }
    }
  }
}

// Threading (core-plan.md step H). The forward transform divides the
// colatitudes between threads and reduces private partial sums; the inverse
// divides them and writes disjoint rows.
TEST(GaussLegendreGrid, ParallelAgreesWithSequential) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr std::ptrdiff_t lMax = 24;
  constexpr std::ptrdiff_t n = 2;
  constexpr auto tolerance = 1.0e-12;

  auto grid = Grid(lMax, n, FFTWpp::Estimate);
  const auto indices = GSHIndices<All>(lMax, lMax, n);

  auto given = FFTWpp::vector<Complex>(indices.Size());
  for (auto j = std::size_t{0}; j < given.size(); ++j) {
    given[j] = Complex{0.3 + 0.01 * j, -0.2 + 0.017 * j};
  }

  auto fieldSequential = FFTWpp::vector<Complex>(grid.FieldSize());
  auto coefficientsSequential = FFTWpp::vector<Complex>(indices.Size());
  grid.InverseTransformation(lMax, n, given, fieldSequential);
  grid.ForwardTransformation(lMax, n, fieldSequential,
                             coefficientsSequential);

  for (auto threads : {1, 2, 3, 4, 8}) {
    const auto policy = Execution::Parallel(threads);

    auto field = FFTWpp::vector<Complex>(grid.FieldSize());
    grid.InverseTransformation(lMax, n, given, field, policy);
    // The inverse writes disjoint rows, each computed exactly as it would be
    // sequentially, so this is exact rather than close.
    for (auto i = std::size_t{0}; i < field.size(); ++i) {
      EXPECT_EQ(field[i], fieldSequential[i])
          << threads << " threads, point " << i;
    }

    auto coefficients = FFTWpp::vector<Complex>(indices.Size());
    grid.ForwardTransformation(lMax, n, field, coefficients, policy);
    // The forward transform sums the colatitudes in a different order, so
    // this is close and not exact.
    for (auto j = std::size_t{0}; j < coefficients.size(); ++j) {
      EXPECT_NEAR(std::abs(coefficients[j] - coefficientsSequential[j]), 0.0,
                  tolerance)
          << threads << " threads, coefficient " << j;
    }
  }
}

// Exactly one level threads. A transform asked to run in parallel from inside
// a parallel region must run sequentially rather than nest.
TEST(GaussLegendreGrid, NestedParallelismIsSuppressed) {
  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  constexpr std::ptrdiff_t lMax = 12;
  constexpr auto tolerance = 1.0e-12;

  auto grid = Grid(lMax, 0, FFTWpp::Estimate);
  const auto indices = GSHIndices<All>(lMax, lMax, 0);
  constexpr auto count = 8;

  auto Given = [&](int seed) {
    auto given = FFTWpp::vector<Complex>(indices.Size());
    for (auto j = std::size_t{0}; j < given.size(); ++j) {
      given[j] = Complex{0.1 * seed + 0.01 * j, -0.05 * seed + 0.013 * j};
    }
    return given;
  };

  auto reference = std::vector<FFTWpp::vector<Complex>>{};
  for (auto s = 0; s < count; ++s) {
    auto given = Given(s);
    auto field = FFTWpp::vector<Complex>(grid.FieldSize());
    auto back = FFTWpp::vector<Complex>(indices.Size());
    grid.InverseTransformation(lMax, 0, given, field);
    grid.ForwardTransformation(lMax, 0, field, back);
    reference.push_back(back);
  }

  // The outer loop owns the parallelism; the inner calls ask for it too and
  // must not get it.
  auto results = std::vector<FFTWpp::vector<Complex>>(
      count, FFTWpp::vector<Complex>(indices.Size()));
  auto sawNesting = 0;
#pragma omp parallel for schedule(static) reduction(+ : sawNesting)
  for (auto s = 0; s < count; ++s) {
    auto given = Given(s);
    auto field = FFTWpp::vector<Complex>(grid.FieldSize());
    grid.InverseTransformation(lMax, 0, given, field,
                               Execution::Parallel(4));
    grid.ForwardTransformation(lMax, 0, field, results[s],
                               Execution::Parallel(4));
    if (omp_get_level() > 1) sawNesting += 1;
  }

  EXPECT_EQ(sawNesting, 0);
  for (auto s = 0; s < count; ++s) {
    for (auto j = std::size_t{0}; j < reference[s].size(); ++j) {
      EXPECT_NEAR(std::abs(results[s][j] - reference[s][j]), 0.0, tolerance)
          << "problem " << s << ", coefficient " << j;
    }
  }
}

TEST(GaussLegendreGrid, ExecutionPolicyDefaultsToSequential) {
  EXPECT_FALSE(Execution::Sequential().IsParallel());
  EXPECT_TRUE(Execution::Parallel().IsParallel());
  EXPECT_TRUE(Execution::Parallel(4).IsParallel());
  EXPECT_EQ(Execution::Parallel(4).Threads(), 4);
  // Zero means "whatever OpenMP would choose", not "no threads".
  EXPECT_EQ(Execution::Parallel().Threads(), 0);
  EXPECT_EQ(Execution::Parallel(0).Threads(), 0);
  EXPECT_EQ(Execution::Sequential().Threads(), 1);
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

//--------------------------------------------------------------------------//
//                        The batched transform                              //
//--------------------------------------------------------------------------//
//
// The oracle throughout is the unbatched call: a batch of k fields must give
// exactly what k separate transforms give, in every layout the descriptor can
// express. Exactly, not approximately -- batching changes the order in which
// the Wigner values are fetched but not the order in which anything is summed,
// so any difference at all would mean the loop had been restructured rather
// than widened (core-plan.md step F, tier 1).

namespace {

using BatchReal = double;
using BatchComplex = std::complex<BatchReal>;
using BatchGrid = GaussLegendreGrid<BatchReal, All, All>;

// Distinct, reproducible data for field k.
auto BatchField(std::ptrdiff_t size, std::ptrdiff_t k) {
  auto field = std::vector<BatchComplex>(size);
  for (auto i = std::ptrdiff_t{0}; i < size; i++) {
    field[i] = BatchComplex{std::cos(0.3 * i + k), std::sin(0.7 * i - 2 * k)};
  }
  return field;
}

}  // namespace

// [C18]: the chunking policy and the planner flag live on the handle, not
// beside the table, so a grid with a different chunk is a pointer copy and
// shares one table. Two properties make that safe, and both are asserted
// rather than argued.
//
// The first is that a chunk is not observable in a result. It is how the
// inner loop schedules itself; changing it must change nothing, exactly, in
// either direction and batched or not.
TEST(BatchedTransform, ChunkingIsNotObservableInAnyResult) {
  constexpr auto lMax = std::ptrdiff_t{6};
  constexpr auto n = std::ptrdiff_t{2};
  constexpr auto count = std::ptrdiff_t{5};

  auto grid = BatchGrid(lMax, n, FFTWpp::Estimate);
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(grid.CoefficientSize(lMax, n));

  auto fields = FFTWpp::vector<BatchComplex>(count * fieldSize);
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    const auto one = BatchField(fieldSize, k);
    std::copy(one.begin(), one.end(), fields.begin() + k * fieldSize);
  }

  const auto Forward = [&](const auto& g) {
    auto out = FFTWpp::vector<BatchComplex>(count * coefficientSize);
    g.ForwardTransformation(lMax, n, fields,
                            Batch::Contiguous(count, fieldSize), out,
                            Batch::Contiguous(count, coefficientSize));
    return out;
  };
  const auto Inverse = [&](const auto& g, const auto& coefficients) {
    auto out = FFTWpp::vector<BatchComplex>(count * fieldSize);
    g.InverseTransformation(lMax, n, coefficients,
                            Batch::Contiguous(count, coefficientSize), out,
                            Batch::Contiguous(count, fieldSize));
    return out;
  };

  const auto reference = Forward(grid);
  const auto back = Inverse(grid, reference);

  for (auto chunk : {std::ptrdiff_t{1}, std::ptrdiff_t{2}, std::ptrdiff_t{5},
                     std::ptrdiff_t{16}}) {
    const auto tuned = grid.With(Chunking::Fixed(chunk));

    // The second property: a differently chunked grid is the *same* grid, so
    // fields and coefficients built on one are usable on the other. That is
    // what lets a tuned grid be substituted for an untuned one.
    EXPECT_EQ(tuned.Identity(), grid.Identity()) << "chunk " << chunk;

    const auto got = Forward(tuned);
    for (auto j = std::size_t{0}; j < got.size(); ++j) {
      EXPECT_EQ(got[j], reference[j]) << "chunk " << chunk << ", forward " << j;
    }
    const auto gotBack = Inverse(tuned, reference);
    for (auto j = std::size_t{0}; j < gotBack.size(); ++j) {
      EXPECT_EQ(gotBack[j], back[j]) << "chunk " << chunk << ", inverse " << j;
    }
  }

  // And it really is the policy that moved, not a copy of the default.
  EXPECT_EQ(grid.With(Chunking::Fixed(3)).ChunkingPolicy(),
            Chunking::Fixed(3));
  EXPECT_EQ(grid.ChunkingPolicy(), Chunking::Automatic());
}

TEST(BatchedTransform, ContiguousBatchMatchesSeparateCalls) {
  constexpr auto lMax = std::ptrdiff_t{6};
  constexpr auto n = std::ptrdiff_t{2};
  constexpr auto count = std::ptrdiff_t{4};

  auto grid = BatchGrid(lMax, n, FFTWpp::Estimate);
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(grid.CoefficientSize(lMax, n));

  auto fields = FFTWpp::vector<BatchComplex>(count * fieldSize);
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    const auto one = BatchField(fieldSize, k);
    std::copy(one.begin(), one.end(), fields.begin() + k * fieldSize);
  }

  auto batched = FFTWpp::vector<BatchComplex>(count * coefficientSize);
  grid.ForwardTransformation(lMax, n, fields,
                             Batch::Contiguous(count, fieldSize), batched,
                             Batch::Contiguous(count, coefficientSize));

  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    auto one = BatchField(fieldSize, k);
    auto expected = FFTWpp::vector<BatchComplex>(coefficientSize);
    grid.ForwardTransformation(lMax, n, one, expected);
    for (auto j = std::ptrdiff_t{0}; j < coefficientSize; j++) {
      EXPECT_EQ(batched[k * coefficientSize + j], expected[j])
          << "field " << k << ", coefficient " << j;
    }
  }
}

TEST(BatchedTransform, InterleavedBatchMatchesSeparateCalls) {
  constexpr auto lMax = std::ptrdiff_t{5};
  constexpr auto n = std::ptrdiff_t{1};
  constexpr auto count = std::ptrdiff_t{3};

  auto grid = BatchGrid(lMax, n, FFTWpp::Estimate);
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(grid.CoefficientSize(lMax, n));

  // Point-major storage five components wide, of which the call touches
  // three: this is the layout the field plan would otherwise have had to
  // repack before transforming ([C9]).
  constexpr auto width = std::ptrdiff_t{5};
  auto fields = FFTWpp::vector<BatchComplex>(width * fieldSize);
  const auto sentinel = BatchComplex{-7.0, 11.0};
  std::ranges::fill(fields, sentinel);
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    const auto one = BatchField(fieldSize, k);
    for (auto i = std::ptrdiff_t{0}; i < fieldSize; i++) {
      fields[i * width + k] = one[i];
    }
  }

  auto batched = FFTWpp::vector<BatchComplex>(width * coefficientSize);
  std::ranges::fill(batched, sentinel);
  grid.ForwardTransformation(lMax, n, fields, Batch::Interleaved(count, width),
                             batched, Batch::Interleaved(count, width));

  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    auto one = BatchField(fieldSize, k);
    auto expected = FFTWpp::vector<BatchComplex>(coefficientSize);
    grid.ForwardTransformation(lMax, n, one, expected);
    for (auto j = std::ptrdiff_t{0}; j < coefficientSize; j++) {
      EXPECT_EQ(batched[j * width + k], expected[j])
          << "field " << k << ", coefficient " << j;
    }
  }

  // The columns outside the batch belong to components this call knows
  // nothing about, and must be exactly as they were left.
  for (auto j = std::ptrdiff_t{0}; j < coefficientSize; j++) {
    for (auto k = count; k < width; k++) {
      EXPECT_EQ(batched[j * width + k], sentinel)
          << "column " << k << " is not part of the call";
    }
  }
}

TEST(BatchedTransform, InverseBatchMatchesSeparateCalls) {
  constexpr auto lMax = std::ptrdiff_t{6};
  constexpr auto n = std::ptrdiff_t{2};
  constexpr auto count = std::ptrdiff_t{3};

  auto grid = BatchGrid(lMax, n, FFTWpp::Estimate);
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(grid.CoefficientSize(lMax, n));

  // Coefficients that came from real fields, so the inverse is asked for
  // something a transform could actually have produced.
  auto coefficients = FFTWpp::vector<BatchComplex>(count * coefficientSize);
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    auto one = BatchField(fieldSize, k);
    auto block = FFTWpp::vector<BatchComplex>(coefficientSize);
    grid.ForwardTransformation(lMax, n, one, block);
    std::copy(block.begin(), block.end(),
              coefficients.begin() + k * coefficientSize);
  }

  auto batched = FFTWpp::vector<BatchComplex>(count * fieldSize);
  grid.InverseTransformation(lMax, n, coefficients,
                             Batch::Contiguous(count, coefficientSize), batched,
                             Batch::Contiguous(count, fieldSize));

  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    auto block = FFTWpp::vector<BatchComplex>(coefficientSize);
    std::copy_n(coefficients.begin() + k * coefficientSize, coefficientSize,
                block.begin());
    auto expected = FFTWpp::vector<BatchComplex>(fieldSize);
    grid.InverseTransformation(lMax, n, block, expected);
    for (auto i = std::ptrdiff_t{0}; i < fieldSize; i++) {
      EXPECT_EQ(batched[k * fieldSize + i], expected[i])
          << "field " << k << ", point " << i;
    }
  }
}

TEST(BatchedTransform, MixedLayoutsAndParallelAgreeToo) {
  constexpr auto lMax = std::ptrdiff_t{8};
  constexpr auto n = std::ptrdiff_t{0};
  constexpr auto count = std::ptrdiff_t{4};
  constexpr auto width = std::ptrdiff_t{6};

  auto grid = BatchGrid(lMax, n, FFTWpp::Estimate);
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(grid.CoefficientSize(lMax, n));

  // Point-major in, component-major out: the two sides are described
  // independently, which is why the call takes two descriptors.
  auto fields = FFTWpp::vector<BatchComplex>(width * fieldSize);
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    const auto one = BatchField(fieldSize, k);
    for (auto i = std::ptrdiff_t{0}; i < fieldSize; i++) {
      fields[i * width + k] = one[i];
    }
  }

  auto sequential = FFTWpp::vector<BatchComplex>(count * coefficientSize);
  auto parallel = FFTWpp::vector<BatchComplex>(count * coefficientSize);
  const auto inBatch = Batch::Interleaved(count, width);
  const auto outBatch = Batch::Contiguous(count, coefficientSize);

  grid.ForwardTransformation(lMax, n, fields, inBatch, sequential, outBatch);
  grid.ForwardTransformation(lMax, n, fields, inBatch, parallel, outBatch,
                             Execution::Parallel(4));

  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    auto one = BatchField(fieldSize, k);
    auto expected = FFTWpp::vector<BatchComplex>(coefficientSize);
    grid.ForwardTransformation(lMax, n, one, expected);
    for (auto j = std::ptrdiff_t{0}; j < coefficientSize; j++) {
      EXPECT_EQ(sequential[k * coefficientSize + j], expected[j]);
    }
  }

  // The partitioned reduction reassociates the colatitude sum, so the
  // parallel path is close rather than equal -- the same relation the
  // unbatched ParallelAgreesWithSequential test pins.
  for (auto i = std::ptrdiff_t{0}; i < count * coefficientSize; i++) {
    EXPECT_NEAR(parallel[i].real(), sequential[i].real(), 1.0e-12);
    EXPECT_NEAR(parallel[i].imag(), sequential[i].imag(), 1.0e-12);
  }
}

TEST(BatchedTransform, RejectsDescriptorsItCannotHonour) {
  constexpr auto lMax = std::ptrdiff_t{4};
  constexpr auto n = std::ptrdiff_t{0};

  auto grid = BatchGrid(lMax, n, FFTWpp::Estimate);
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(grid.CoefficientSize(lMax, n));

  auto fields = FFTWpp::vector<BatchComplex>(4 * fieldSize);
  auto coefficients = FFTWpp::vector<BatchComplex>(4 * coefficientSize);

  // Counts must agree between the two sides.
  EXPECT_THROW(grid.ForwardTransformation(
                   lMax, n, fields, Batch::Contiguous(3, fieldSize),
                   coefficients, Batch::Contiguous(2, coefficientSize)),
               std::invalid_argument);

  // A range too short for the span the batch describes.
  auto tooSmall = FFTWpp::vector<BatchComplex>(2 * coefficientSize);
  EXPECT_THROW(grid.ForwardTransformation(
                   lMax, n, fields, Batch::Contiguous(4, fieldSize), tooSmall,
                   Batch::Contiguous(4, coefficientSize)),
               std::invalid_argument);

  // Members that overlap: the transform writes every element of every field,
  // so this would silently give a wrong answer.
  EXPECT_THROW(grid.ForwardTransformation(
                   lMax, n, fields, Batch::Strided(4, 1, 2), coefficients,
                   Batch::Contiguous(4, coefficientSize)),
               std::invalid_argument);
}

TEST(BatchedTransform, SingleFieldIsTheBatchAtCountOne) {
  constexpr auto lMax = std::ptrdiff_t{6};
  constexpr auto n = std::ptrdiff_t{1};

  auto grid = BatchGrid(lMax, n, FFTWpp::Estimate);
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(grid.CoefficientSize(lMax, n));

  auto field = BatchField(fieldSize, 0);
  auto viaWrapper = FFTWpp::vector<BatchComplex>(coefficientSize);
  auto viaBatch = FFTWpp::vector<BatchComplex>(coefficientSize);

  grid.ForwardTransformation(lMax, n, field, viaWrapper);
  grid.ForwardTransformation(lMax, n, field, Batch::One(fieldSize), viaBatch,
                             Batch::One(coefficientSize));

  for (auto j = std::ptrdiff_t{0}; j < coefficientSize; j++) {
    EXPECT_EQ(viaWrapper[j], viaBatch[j]);
  }

  // The wrapper keeps the equality check its contract promises, where the
  // batched entry would accept a longer range.
  auto tooLong = FFTWpp::vector<BatchComplex>(coefficientSize + 1);
  EXPECT_THROW(grid.ForwardTransformation(lMax, n, field, tooLong),
               std::invalid_argument);
  EXPECT_NO_THROW(grid.ForwardTransformation(
      lMax, n, field, Batch::One(fieldSize), tooLong,
      Batch::One(coefficientSize)));
}

TEST(BatchedTransform, ChunkingDoesNotChangeTheAnswer) {
  constexpr auto lMax = std::ptrdiff_t{6};
  constexpr auto n = std::ptrdiff_t{2};
  constexpr auto count = std::ptrdiff_t{5};

  // Five fields in chunks of two: two full chunks and a short one, so the
  // boundary the loop has to get right is exercised rather than assumed. At
  // these degrees the automatic policy would take the whole batch at once,
  // which is why the chunk is pinned instead.
  auto chunked = BatchGrid(lMax, n, FFTWpp::Estimate, Chunking::Fixed(2));
  auto whole = BatchGrid(lMax, n, FFTWpp::Estimate, Chunking::Fixed(count));

  const auto fieldSize = static_cast<std::ptrdiff_t>(chunked.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(chunked.CoefficientSize(lMax, n));

  auto fields = FFTWpp::vector<BatchComplex>(count * fieldSize);
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    const auto one = BatchField(fieldSize, k);
    std::copy(one.begin(), one.end(), fields.begin() + k * fieldSize);
  }

  const auto inBatch = Batch::Contiguous(count, fieldSize);
  const auto outBatch = Batch::Contiguous(count, coefficientSize);
  auto inChunks = FFTWpp::vector<BatchComplex>(count * coefficientSize);
  auto inOne = FFTWpp::vector<BatchComplex>(count * coefficientSize);

  chunked.ForwardTransformation(lMax, n, fields, inBatch, inChunks, outBatch);
  whole.ForwardTransformation(lMax, n, fields, inBatch, inOne, outBatch);

  // Chunking partitions the batch; it does not touch the order of any sum.
  for (auto i = std::ptrdiff_t{0}; i < count * coefficientSize; i++) {
    EXPECT_EQ(inChunks[i], inOne[i]) << "element " << i;
  }

  // And the same for the inverse, back to the fields it came from.
  auto backChunked = FFTWpp::vector<BatchComplex>(count * fieldSize);
  auto backWhole = FFTWpp::vector<BatchComplex>(count * fieldSize);
  chunked.InverseTransformation(lMax, n, inChunks, outBatch, backChunked,
                                inBatch);
  whole.InverseTransformation(lMax, n, inOne, outBatch, backWhole, inBatch);
  for (auto i = std::ptrdiff_t{0}; i < count * fieldSize; i++) {
    EXPECT_EQ(backChunked[i], backWhole[i]) << "element " << i;
  }
}

TEST(BatchedTransform, ChunkingPolicyIsCarriedByTheGrid) {
  constexpr auto lMax = std::ptrdiff_t{6};
  constexpr auto n = std::ptrdiff_t{0};

  // A copied grid shares the implementation, and with it the policy.
  auto grid = BatchGrid(lMax, n, FFTWpp::Estimate, Chunking::Fixed(2));
  auto copy = grid;
  EXPECT_EQ(copy.Identity(), grid.Identity());

  // ForBand forwards it too, rather than silently resetting to Automatic.
  auto banded = BatchGrid::ForBand(4, n, 1.5, FFTWpp::Estimate,
                                   Chunking::Fixed(3));
  EXPECT_EQ(banded.MaxDegree(), 6);
}

//--------------------------------------------------------------------------//
//                          The generating grid                              //
//--------------------------------------------------------------------------//
//
// A grid asked for WignerValues::Generated() builds no table and runs the
// recursion inside each transform instead (core-plan.md step F', T11).
//
// The oracle throughout is the stored grid, and the comparison is **exact**.
// This is the same recursion, seeded the same way, evaluated in the same order
// and with the same rounding; only where the values are put differs. So a
// difference of one bit would mean the two paths had diverged arithmetically,
// which is precisely what must not happen -- and it is a far sharper
// instrument than agreeing to a round-trip tolerance, which would pass on
// almost any plausible mistake in the seeding or the truncation.

namespace {

using GenReal = double;
using GenComplex = std::complex<GenReal>;
using GenGrid = GaussLegendreGrid<GenReal, All, All>;

auto GenField(std::ptrdiff_t size, std::ptrdiff_t k = 0) {
  auto field = std::vector<GenComplex>(size);
  for (auto i = std::ptrdiff_t{0}; i < size; i++) {
    field[i] = GenComplex{std::cos(0.31 * i + k), std::sin(0.53 * i - k)};
  }
  return field;
}

auto GenCoefficients(std::ptrdiff_t size, std::ptrdiff_t k = 0) {
  auto c = std::vector<GenComplex>(size);
  for (auto j = std::ptrdiff_t{0}; j < size; j++) {
    c[j] = GenComplex{std::sin(0.17 * j - k), std::cos(0.41 * j + 2 * k)};
  }
  return c;
}

}  // namespace

TEST(GeneratingGrid, ForwardAgreesWithTheStoredTableExactly) {
  for (auto [lMax, n] : std::vector<std::pair<std::ptrdiff_t, std::ptrdiff_t>>{
           {7, 2}, {12, 0}, {9, -3}, {5, 5}}) {
    auto stored = GenGrid(lMax, std::abs(n), FFTWpp::Estimate);
    auto generated =
        GenGrid(lMax, std::abs(n), FFTWpp::Estimate, Chunking::Automatic(),
                WignerValues::Generated());

    const auto fieldSize = static_cast<std::ptrdiff_t>(stored.FieldSize());
    const auto coefficientSize =
        static_cast<std::ptrdiff_t>(stored.CoefficientSize(lMax, n));
    ASSERT_EQ(fieldSize, static_cast<std::ptrdiff_t>(generated.FieldSize()));

    const auto field = GenField(fieldSize);
    auto a = std::vector<GenComplex>(coefficientSize);
    auto b = std::vector<GenComplex>(coefficientSize);
    stored.ForwardTransformation(lMax, n, field, a);
    generated.ForwardTransformation(lMax, n, field, b);

    for (auto j = std::ptrdiff_t{0}; j < coefficientSize; j++) {
      EXPECT_EQ(a[j], b[j])
          << "lMax = " << lMax << ", n = " << n << ", coefficient " << j;
    }
  }
}

TEST(GeneratingGrid, InverseAgreesWithTheStoredTableExactly) {
  constexpr auto lMax = std::ptrdiff_t{8};
  constexpr auto n = std::ptrdiff_t{2};

  auto stored = GenGrid(lMax, n, FFTWpp::Estimate);
  auto generated = GenGrid(lMax, n, FFTWpp::Estimate, Chunking::Automatic(),
                           WignerValues::Generated());

  const auto fieldSize = static_cast<std::ptrdiff_t>(stored.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(stored.CoefficientSize(lMax, n));

  const auto coefficients = GenCoefficients(coefficientSize);
  auto a = std::vector<GenComplex>(fieldSize);
  auto b = std::vector<GenComplex>(fieldSize);
  stored.InverseTransformation(lMax, n, coefficients, a);
  generated.InverseTransformation(lMax, n, coefficients, b);

  for (auto i = std::ptrdiff_t{0}; i < fieldSize; i++) {
    EXPECT_EQ(a[i], b[i]) << "sample " << i;
  }
}

// A truncated call is where the two paths could most easily diverge: the
// stored grid reads a prefix of rows laid out for its own maximum degree,
// while the generating one builds a block sized for the call. They must lay
// out the same values in the same places.
TEST(GeneratingGrid, TruncatedCallsGenerateOnlyTheDegreesTheyUse) {
  constexpr auto gridDegree = std::ptrdiff_t{11};
  constexpr auto n = std::ptrdiff_t{1};

  auto stored = GenGrid(gridDegree, n, FFTWpp::Estimate);
  auto generated =
      GenGrid(gridDegree, n, FFTWpp::Estimate, Chunking::Automatic(),
              WignerValues::Generated());

  const auto field = GenField(static_cast<std::ptrdiff_t>(stored.FieldSize()));

  for (auto lMax = std::abs(n); lMax <= gridDegree; lMax++) {
    const auto coefficientSize =
        static_cast<std::ptrdiff_t>(stored.CoefficientSize(lMax, n));
    auto a = std::vector<GenComplex>(coefficientSize);
    auto b = std::vector<GenComplex>(coefficientSize);
    stored.ForwardTransformation(lMax, n, field, a);
    generated.ForwardTransformation(lMax, n, field, b);
    for (auto j = std::ptrdiff_t{0}; j < coefficientSize; j++) {
      EXPECT_EQ(a[j], b[j]) << "lMax = " << lMax << ", coefficient " << j;
    }
  }
}

// The reduced m >= 0 path of a real scalar grid, which indexes the generated
// block differently from the complex one.
TEST(GeneratingGrid, RealScalarGridAgreesExactlyToo) {
  constexpr auto lMax = std::ptrdiff_t{9};
  using RealGrid = GaussLegendreGrid<GenReal, NonNegative, All>;

  auto stored = RealGrid(lMax, 0, FFTWpp::Estimate);
  auto generated = RealGrid(lMax, 0, FFTWpp::Estimate, Chunking::Automatic(),
                            WignerValues::Generated());

  const auto fieldSize = static_cast<std::ptrdiff_t>(stored.FieldSize());

  // The reduced m >= 0 storage, which is the whole point of a real transform
  // and is not what CoefficientSize reports.
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(stored.RealCoefficientSize(lMax));

  auto field = std::vector<GenReal>(fieldSize);
  for (auto i = std::ptrdiff_t{0}; i < fieldSize; i++) {
    field[i] = std::cos(0.23 * i) + 0.5 * std::sin(0.11 * i);
  }

  auto a = std::vector<GenComplex>(coefficientSize);
  auto b = std::vector<GenComplex>(coefficientSize);
  stored.ForwardTransformation(lMax, 0, field, a);
  generated.ForwardTransformation(lMax, 0, field, b);
  for (auto j = std::ptrdiff_t{0}; j < coefficientSize; j++) {
    EXPECT_EQ(a[j], b[j]) << "coefficient " << j;
  }

  auto back = std::vector<GenReal>(fieldSize);
  auto backGenerated = std::vector<GenReal>(fieldSize);
  stored.InverseTransformation(lMax, 0, a, back);
  generated.InverseTransformation(lMax, 0, b, backGenerated);
  for (auto i = std::ptrdiff_t{0}; i < fieldSize; i++) {
    EXPECT_EQ(back[i], backGenerated[i]) << "sample " << i;
  }
}

// Batching and threading are where the per-thread scratch is actually
// exercised: every thread must generate its own block and no two may share
// one.
TEST(GeneratingGrid, BatchedAndParallelCallsAgreeExactly) {
  constexpr auto lMax = std::ptrdiff_t{10};
  constexpr auto n = std::ptrdiff_t{2};
  constexpr auto count = std::ptrdiff_t{5};

  auto stored = GenGrid(lMax, n, FFTWpp::Estimate, Chunking::Fixed(2));
  auto generated = GenGrid(lMax, n, FFTWpp::Estimate, Chunking::Fixed(2),
                           WignerValues::Generated());

  const auto fieldSize = static_cast<std::ptrdiff_t>(stored.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(stored.CoefficientSize(lMax, n));

  auto fields = std::vector<GenComplex>(count * fieldSize);
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    const auto one = GenField(fieldSize, k);
    std::copy(one.begin(), one.end(), fields.begin() + k * fieldSize);
  }

  // Like for like in both runs. The forward transform's parallel reduction
  // sums each thread's partials at the end rather than accumulating the
  // colatitudes in order, so a sequential run and a parallel one differ in
  // the last bits by design; comparing across that would be testing the
  // reduction, not the supplier.
  const auto Run = [&](auto& grid, Execution policy) {
    auto out = std::vector<GenComplex>(count * coefficientSize);
    grid.ForwardTransformation(lMax, n, fields,
                               Batch::Contiguous(count, fieldSize), out,
                               Batch::Contiguous(count, coefficientSize),
                               policy);
    return out;
  };

  for (auto policy : {Execution::Sequential(), Execution::Parallel(4)}) {
    const auto a = Run(stored, policy);
    const auto b = Run(generated, policy);
    for (auto j = std::ptrdiff_t{0}; j < count * coefficientSize; j++) {
      EXPECT_EQ(a[j], b[j])
          << "entry " << j << (policy.IsParallel() ? ", parallel" : "");
    }
  }
}

TEST(GeneratingGrid, BuildsNoTableAndForBandCarriesThePolicy) {
  constexpr auto lMax = std::ptrdiff_t{64};

  // Construction is where the difference shows: the stored grid pays for the
  // whole table and the generating one pays for two vectors of 2 lMax + 1.
  auto generated = GenGrid(lMax, 2, FFTWpp::Estimate, Chunking::Automatic(),
                           WignerValues::Generated());
  EXPECT_EQ(generated.MaxDegree(), lMax);

  // ForBand forwards the policy rather than silently resetting it to Stored.
  auto banded = GenGrid::ForBand(8, 2, 1.5, FFTWpp::Estimate,
                                 Chunking::Automatic(),
                                 WignerValues::Generated());
  EXPECT_EQ(banded.MaxDegree(), 12);

  auto stored = GenGrid(12, 2, FFTWpp::Estimate);
  const auto fieldSize = static_cast<std::ptrdiff_t>(stored.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(stored.CoefficientSize(12, 2));
  const auto field = GenField(fieldSize);
  auto a = std::vector<GenComplex>(coefficientSize);
  auto b = std::vector<GenComplex>(coefficientSize);
  stored.ForwardTransformation(12, 2, field, a);
  banded.ForwardTransformation(12, 2, field, b);
  for (auto j = std::ptrdiff_t{0}; j < coefficientSize; j++) {
    EXPECT_EQ(a[j], b[j]) << "coefficient " << j;
  }
}

// The two directions decompose differently, so they must not share a chunk
// rule. What the policy needs is the number of *copies* of the coefficient
// block that will be live at once: the forward transform gives every thread a
// private accumulator, the inverse gathers one shared read-only block.
//
// Serving both with the thread count starved the inverse -- a chunk of one at
// lMax = 256 and k = 8 on eight threads, where the whole batch fits, which
// measured 2.2x slower (core-plan.md section 10).
TEST(BatchedTransform, ChunkRuleCountsCopiesNotThreads) {
  // One field's coefficients at lMax = 256, n = 2.
  constexpr auto bytesPerField = std::ptrdiff_t{66045} * 16;

  const auto laptop = Chunking::ForCache(std::ptrdiff_t{16} << 20);
  EXPECT_GT(laptop.Count(bytesPerField, 1), laptop.Count(bytesPerField, 8))
      << "a shared block should be allowed a wider chunk than a private one";

  // Both of the anchors the formula was built on, and both are
  // forward-shaped: P2's optimum of eight was measured sequentially, so one
  // copy had the whole cache, and P8's prediction of two is sixty-four private
  // accumulators sharing 256 MiB.
  EXPECT_EQ(laptop.Count(bytesPerField, 1), 8);
  EXPECT_EQ((Chunking::ForCache(std::ptrdiff_t{256} << 20)
                 .Count(bytesPerField, 64)),
            2);

  // Fixed still defeats the heuristic from either side.
  EXPECT_EQ(Chunking::Fixed(5).Count(bytesPerField, 1), 5);
  EXPECT_EQ(Chunking::Fixed(5).Count(bytesPerField, 64), 5);
}

// -- The matrix kernel's Fourier stage (core-plan.md section 11, step M2).
//
// Section 11 asked for this to be checked against the loop kernel's own FFT
// stage. It is checked against a naive DFT written out here instead, which is
// strictly stronger: comparing two paths through the same FFTW plan would
// agree even if the plan were the wrong transform, and would say nothing at
// all about the sign of the exponent, the normalisation, or where negative
// orders live. An independent sum says all three.
namespace {

using StageReal = double;
using StageComplex = std::complex<StageReal>;
using StageGrid = GaussLegendreGrid<StageReal, All, All>;

// out[m][iTheta][k] by direct summation. FFTW's forward transform carries no
// normalisation and a negative exponent, and for a complex field its order m
// sits at index m for m >= 0 and at nPhi - |m| for m < 0 -- which is what
// running the index from 0 to nPhi - 1 says.
auto NaiveFourier(const std::vector<StageComplex>& fields, std::ptrdiff_t nTheta,
                  std::ptrdiff_t nPhi, std::ptrdiff_t count,
                  std::ptrdiff_t nFourier) {
  auto out = std::vector<StageComplex>(nFourier * nTheta * count);
  for (auto m = std::ptrdiff_t{0}; m < nFourier; m++) {
    for (auto iTheta = std::ptrdiff_t{0}; iTheta < nTheta; iTheta++) {
      for (auto k = std::ptrdiff_t{0}; k < count; k++) {
        auto sum = StageComplex{0, 0};
        for (auto iPhi = std::ptrdiff_t{0}; iPhi < nPhi; iPhi++) {
          const auto angle = -2 * std::numbers::pi_v<StageReal> *
                             static_cast<StageReal>(iPhi * m) /
                             static_cast<StageReal>(nPhi);
          sum += fields[k * nTheta * nPhi + iTheta * nPhi + iPhi] *
                 StageComplex{std::cos(angle), std::sin(angle)};
        }
        out[m * nTheta * count + iTheta * count + k] = sum;
      }
    }
  }
  return out;
}

}  // namespace

TEST(FourierStage, MatchesADirectTransform) {
  constexpr auto lMax = std::ptrdiff_t{5};
  constexpr auto count = std::ptrdiff_t{3};

  auto grid = StageGrid(lMax, 2, FFTWpp::Estimate);
  const auto nTheta = static_cast<std::ptrdiff_t>(grid.NumberOfCoLatitudes());
  const auto nPhi = static_cast<std::ptrdiff_t>(grid.NumberOfLongitudes());
  const auto fieldSize = nTheta * nPhi;

  auto fields = std::vector<StageComplex>(count * fieldSize);
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    const auto one = BatchField(fieldSize, k);
    std::copy(one.begin(), one.end(), fields.begin() + k * fieldSize);
  }

  const auto nFourier = grid.FourierSize<StageComplex>();
  EXPECT_EQ(nFourier, nPhi);

  auto stage = std::vector<StageComplex>(
      grid.ForwardFourierStageSize<StageComplex>(count));
  grid.ForwardFourierStage(fields, Batch::Contiguous(count, fieldSize), 0,
                           count, std::span<StageComplex>(stage));

  const auto expected = NaiveFourier(fields, nTheta, nPhi, count, nFourier);
  ASSERT_EQ(stage.size(), expected.size());
  for (std::size_t i = 0; i < stage.size(); i++) {
    EXPECT_NEAR(stage[i].real(), expected[i].real(), 1e-12) << "at " << i;
    EXPECT_NEAR(stage[i].imag(), expected[i].imag(), 1e-12) << "at " << i;
  }
}

// A real field's stage holds nPhi / 2 + 1 orders, the negative ones being the
// conjugates of the positive. The same direct sum, over the orders that are
// stored.
TEST(FourierStage, MatchesADirectTransformForARealField) {
  constexpr auto lMax = std::ptrdiff_t{5};
  constexpr auto count = std::ptrdiff_t{2};

  auto grid = StageGrid(lMax, 0, FFTWpp::Estimate);
  const auto nTheta = static_cast<std::ptrdiff_t>(grid.NumberOfCoLatitudes());
  const auto nPhi = static_cast<std::ptrdiff_t>(grid.NumberOfLongitudes());
  const auto fieldSize = nTheta * nPhi;

  auto fields = std::vector<StageReal>(count * fieldSize);
  auto asComplex = std::vector<StageComplex>(count * fieldSize);
  for (std::size_t i = 0; i < fields.size(); i++) {
    fields[i] = std::cos(0.11 * static_cast<StageReal>(i)) +
                0.5 * std::sin(0.37 * static_cast<StageReal>(i));
    asComplex[i] = StageComplex{fields[i], 0};
  }

  const auto nFourier = grid.FourierSize<StageReal>();
  EXPECT_EQ(nFourier, nPhi / 2 + 1);

  auto stage = std::vector<StageComplex>(
      grid.ForwardFourierStageSize<StageReal>(count));
  grid.ForwardFourierStage(fields, Batch::Contiguous(count, fieldSize), 0,
                           count, std::span<StageComplex>(stage));

  const auto expected = NaiveFourier(asComplex, nTheta, nPhi, count, nFourier);
  ASSERT_EQ(stage.size(), expected.size());
  for (std::size_t i = 0; i < stage.size(); i++) {
    EXPECT_NEAR(stage[i].real(), expected[i].real(), 1e-12) << "at " << i;
    EXPECT_NEAR(stage[i].imag(), expected[i].imag(), 1e-12) << "at " << i;
  }
}

// Blocking bounds the FFT workspace and must not otherwise exist. Bit-exact,
// not near: a block boundary that lost or duplicated a row would show up as a
// wrong number, and a block size that changed the answer at all would mean
// the copy back into place is reading the wrong run.
TEST(FourierStage, BlockingChangesNothing) {
  constexpr auto lMax = std::ptrdiff_t{7};
  constexpr auto count = std::ptrdiff_t{3};

  auto grid = StageGrid(lMax, 1, FFTWpp::Estimate);
  const auto nTheta = static_cast<std::ptrdiff_t>(grid.NumberOfCoLatitudes());
  const auto nPhi = static_cast<std::ptrdiff_t>(grid.NumberOfLongitudes());
  const auto fieldSize = nTheta * nPhi;

  auto fields = std::vector<StageComplex>(count * fieldSize);
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    const auto one = BatchField(fieldSize, k);
    std::copy(one.begin(), one.end(), fields.begin() + k * fieldSize);
  }

  const auto size = grid.ForwardFourierStageSize<StageComplex>(count);
  auto whole = std::vector<StageComplex>(size);
  grid.ForwardFourierStage(fields, Batch::Contiguous(count, fieldSize), 0,
                           count, std::span<StageComplex>(whole));

  for (auto block : {std::ptrdiff_t{1}, std::ptrdiff_t{2}, std::ptrdiff_t{3},
                     nTheta - 1, nTheta, nTheta + 5}) {
    auto blocked = std::vector<StageComplex>(size);
    grid.ForwardFourierStage(fields, Batch::Contiguous(count, fieldSize), 0,
                             count, std::span<StageComplex>(blocked), block);
    for (std::size_t i = 0; i < whole.size(); i++) {
      ASSERT_EQ(blocked[i], whole[i]) << "block " << block << ", at " << i;
    }
  }
}

// The caller's stride enters at PackRow and nowhere else, so an interleaved
// batch -- tensor components stored point by point -- must give the same
// answer as the same fields laid end to end. This is the [C9] path, and it is
// the one that would break silently if the pack were written against
// contiguity.
TEST(FourierStage, InterleavedBatchMatchesContiguous) {
  constexpr auto lMax = std::ptrdiff_t{5};
  constexpr auto count = std::ptrdiff_t{2};
  constexpr auto width = std::ptrdiff_t{5};

  auto grid = StageGrid(lMax, 1, FFTWpp::Estimate);
  const auto nTheta = static_cast<std::ptrdiff_t>(grid.NumberOfCoLatitudes());
  const auto nPhi = static_cast<std::ptrdiff_t>(grid.NumberOfLongitudes());
  const auto fieldSize = nTheta * nPhi;

  auto contiguous = std::vector<StageComplex>(count * fieldSize);
  auto interleaved = std::vector<StageComplex>(width * fieldSize, {0, 0});
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    const auto one = BatchField(fieldSize, k);
    std::copy(one.begin(), one.end(), contiguous.begin() + k * fieldSize);
    for (auto j = std::ptrdiff_t{0}; j < fieldSize; j++) {
      interleaved[j * width + k] = one[j];
    }
  }

  const auto size = grid.ForwardFourierStageSize<StageComplex>(count);
  auto fromContiguous = std::vector<StageComplex>(size);
  auto fromInterleaved = std::vector<StageComplex>(size);
  grid.ForwardFourierStage(contiguous, Batch::Contiguous(count, fieldSize), 0,
                           count, std::span<StageComplex>(fromContiguous));
  grid.ForwardFourierStage(interleaved, Batch::Interleaved(count, width), 0,
                           count, std::span<StageComplex>(fromInterleaved));

  for (std::size_t i = 0; i < fromContiguous.size(); i++) {
    EXPECT_EQ(fromInterleaved[i], fromContiguous[i]) << "at " << i;
  }
}

// A sub-range of a batch, which is what the chunking of step M3 will hand it.
TEST(FourierStage, TransformsASubRangeOfTheBatch) {
  constexpr auto lMax = std::ptrdiff_t{4};
  constexpr auto count = std::ptrdiff_t{4};

  auto grid = StageGrid(lMax, 1, FFTWpp::Estimate);
  const auto nTheta = static_cast<std::ptrdiff_t>(grid.NumberOfCoLatitudes());
  const auto nPhi = static_cast<std::ptrdiff_t>(grid.NumberOfLongitudes());
  const auto fieldSize = nTheta * nPhi;

  auto fields = std::vector<StageComplex>(count * fieldSize);
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    const auto one = BatchField(fieldSize, k);
    std::copy(one.begin(), one.end(), fields.begin() + k * fieldSize);
  }

  // Fields two and three of four, against the same two transformed alone.
  auto pair = std::vector<StageComplex>(2 * fieldSize);
  std::copy(fields.begin() + 2 * fieldSize, fields.begin() + 4 * fieldSize,
            pair.begin());

  const auto size = grid.ForwardFourierStageSize<StageComplex>(2);
  auto fromSubRange = std::vector<StageComplex>(size);
  auto fromOwnBatch = std::vector<StageComplex>(size);
  grid.ForwardFourierStage(fields, Batch::Contiguous(count, fieldSize), 2, 2,
                           std::span<StageComplex>(fromSubRange));
  grid.ForwardFourierStage(pair, Batch::Contiguous(2, fieldSize), 0, 2,
                           std::span<StageComplex>(fromOwnBatch));

  for (std::size_t i = 0; i < size; i++) {
    EXPECT_EQ(fromSubRange[i], fromOwnBatch[i]) << "at " << i;
  }
}

TEST(FourierStage, RejectsBadRequests) {
  constexpr auto lMax = std::ptrdiff_t{4};
  auto grid = StageGrid(lMax, 1, FFTWpp::Estimate);
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());

  auto fields = std::vector<StageComplex>(2 * fieldSize);
  auto out = std::vector<StageComplex>(
      grid.ForwardFourierStageSize<StageComplex>(2));
  const auto batch = Batch::Contiguous(2, fieldSize);

  // Fields outside the batch.
  EXPECT_THROW(grid.ForwardFourierStage(fields, batch, 1, 2,
                                        std::span<StageComplex>(out)),
               std::invalid_argument);
  EXPECT_THROW(grid.ForwardFourierStage(fields, batch, -1, 2,
                                        std::span<StageComplex>(out)),
               std::invalid_argument);
  // A count that is not a count.
  EXPECT_THROW(grid.ForwardFourierStage(fields, batch, 0, 0,
                                        std::span<StageComplex>(out)),
               std::invalid_argument);
  // An output buffer sized for the wrong number of fields.
  auto tooSmall = std::vector<StageComplex>(
      grid.ForwardFourierStageSize<StageComplex>(1));
  EXPECT_THROW(grid.ForwardFourierStage(fields, batch, 0, 2,
                                        std::span<StageComplex>(tooSmall)),
               std::invalid_argument);
}

// The aliasing guard shrinks a block whose output stride would be a power of
// two. That is a speed decision and must be invisible in the answers, so what
// is checked here is that a request which triggers it still agrees with one
// that does not: eight fields with a block of eight gives an output stride of
// 1024 bytes, which is the case measured at three times its neighbours'.
TEST(FourierStage, TheAliasingGuardChangesNoAnswers) {
  constexpr auto lMax = std::ptrdiff_t{6};
  constexpr auto count = std::ptrdiff_t{8};

  auto grid = StageGrid(lMax, 1, FFTWpp::Estimate);
  const auto nTheta = static_cast<std::ptrdiff_t>(grid.NumberOfCoLatitudes());
  const auto nPhi = static_cast<std::ptrdiff_t>(grid.NumberOfLongitudes());
  const auto fieldSize = nTheta * nPhi;

  auto fields = std::vector<StageComplex>(count * fieldSize);
  for (auto k = std::ptrdiff_t{0}; k < count; k++) {
    const auto one = BatchField(fieldSize, k);
    std::copy(one.begin(), one.end(), fields.begin() + k * fieldSize);
  }
  const auto batch = Batch::Contiguous(count, fieldSize);
  const auto size = grid.ForwardFourierStageSize<StageComplex>(count);

  auto reference = std::vector<StageComplex>(size);
  grid.ForwardFourierStage(fields, batch, 0, count,
                           std::span<StageComplex>(reference), 1);

  // 8 aliases at count = 8 and is shrunk; 16 and 32 likewise; 3 does not.
  for (auto block : {std::ptrdiff_t{3}, std::ptrdiff_t{8}, std::ptrdiff_t{16},
                     std::ptrdiff_t{32}}) {
    auto guarded = std::vector<StageComplex>(size);
    grid.ForwardFourierStage(fields, batch, 0, count,
                             std::span<StageComplex>(guarded), block);
    for (std::size_t i = 0; i < size; i++) {
      ASSERT_EQ(guarded[i], reference[i]) << "block " << block << ", at " << i;
    }
  }

  // And the default, which is what a caller who says nothing gets.
  auto byDefault = std::vector<StageComplex>(size);
  grid.ForwardFourierStage(fields, batch, 0, count,
                           std::span<StageComplex>(byDefault));
  for (std::size_t i = 0; i < size; i++) {
    ASSERT_EQ(byDefault[i], reference[i]) << "default, at " << i;
  }
}

// -- The matrix kernel (core-plan.md section 11, step M3).
//
// [C12] keeps both kernels permanently, and this is what that buys: the same
// inputs through two independent arrangements of the same sum. A GEMM sums in
// whatever order its kernel chooses, so this cannot be a bit comparison the
// way the batched-against-unbatched tests are -- but to a tolerance it checks
// the layout, the indexing, the FFT ordering and the accumulation together,
// which is very nearly everything the restructure can get wrong.
//
// The two are not independent in the d-values themselves: both read the same
// recursion. That half is pinned separately, by CheckWignerConvention against
// the l = 1 table and CheckLegendre against std::sph_legendre.
#ifdef GSHTRANS_HAVE_BLAS
namespace {

template <typename Grid, typename Scalar>
void CheckKernelsAgreeForward(std::ptrdiff_t lMax, std::ptrdiff_t gridDegree,
                              std::ptrdiff_t n, std::ptrdiff_t count,
                              double tolerance) {
  auto loop = Grid(gridDegree, std::abs(n), FFTWpp::Estimate);
  auto matrix = Grid(gridDegree, std::abs(n), FFTWpp::Estimate,
                     Chunking::Automatic(), WignerValues::Stored(),
                     TransformKernel::Matrix());

  const auto fieldSize = static_cast<std::ptrdiff_t>(loop.FieldSize());
  const auto coefficientSize = static_cast<std::ptrdiff_t>(
      std::is_same_v<Scalar, double> ? loop.RealCoefficientSize(lMax)
                                     : loop.CoefficientSize(lMax, n));

  auto fields = std::vector<Scalar>(count * fieldSize);
  for (std::size_t i = 0; i < fields.size(); i++) {
    if constexpr (std::is_same_v<Scalar, std::complex<double>>) {
      fields[i] = Scalar{std::cos(0.31 * static_cast<double>(i)),
                         std::sin(0.17 * static_cast<double>(i))};
    } else {
      fields[i] = std::cos(0.31 * static_cast<double>(i));
    }
  }

  auto fromLoop = std::vector<std::complex<double>>(count * coefficientSize);
  auto fromMatrix = std::vector<std::complex<double>>(count * coefficientSize);
  const auto fieldBatch = Batch::Contiguous(count, fieldSize);
  const auto coeffBatch = Batch::Contiguous(count, coefficientSize);

  loop.ForwardTransformation(lMax, n, fields, fieldBatch, fromLoop, coeffBatch);
  matrix.ForwardTransformation(lMax, n, fields, fieldBatch, fromMatrix,
                               coeffBatch);

  // Scaled to the largest coefficient, since the absolute size of these
  // depends on the data and a fixed tolerance would be a statement about the
  // test field rather than about the kernels.
  double scale = 0;
  for (const auto& z : fromLoop) scale = std::max(scale, std::abs(z));
  ASSERT_GT(scale, 0.0);

  for (std::size_t i = 0; i < fromLoop.size(); i++) {
    EXPECT_NEAR(fromMatrix[i].real(), fromLoop[i].real(), tolerance * scale)
        << "at " << i;
    EXPECT_NEAR(fromMatrix[i].imag(), fromLoop[i].imag(), tolerance * scale)
        << "at " << i;
  }
}

}  // namespace

TEST(MatrixKernel, ForwardAgreesWithTheLoopKernel) {
  using Grid = GaussLegendreGrid<double, All, All>;
  CheckKernelsAgreeForward<Grid, std::complex<double>>(8, 8, 0, 1, 1e-13);
  CheckKernelsAgreeForward<Grid, std::complex<double>>(8, 8, 2, 1, 1e-13);
  CheckKernelsAgreeForward<Grid, std::complex<double>>(8, 8, -2, 1, 1e-13);
}

TEST(MatrixKernel, ForwardAgreesOverABatch) {
  using Grid = GaussLegendreGrid<double, All, All>;
  CheckKernelsAgreeForward<Grid, std::complex<double>>(10, 10, 1, 5, 1e-13);
  CheckKernelsAgreeForward<Grid, std::complex<double>>(16, 16, 2, 8, 1e-13);
}

// A transform below the grid's own degree takes a prefix of each matrix's
// rows. Nothing is copied for it, so getting the prefix wrong would be
// invisible except here.
TEST(MatrixKernel, ForwardAgreesBelowTheGridDegree) {
  using Grid = GaussLegendreGrid<double, All, All>;
  CheckKernelsAgreeForward<Grid, std::complex<double>>(6, 12, 1, 3, 1e-13);
  CheckKernelsAgreeForward<Grid, std::complex<double>>(1, 12, 0, 2, 1e-13);
}

// The reduced m >= 0 storage, which is a different coefficient layout and a
// different FFT on the way in.
TEST(MatrixKernel, ForwardAgreesForARealField) {
  using Grid = GaussLegendreGrid<double, All, All>;
  CheckKernelsAgreeForward<Grid, double>(8, 8, 0, 1, 1e-13);
  CheckKernelsAgreeForward<Grid, double>(12, 12, 0, 4, 1e-13);
}

// A grid that stores only m >= 0 -- the real scalar grid of [C1] -- has a
// different Wigner layout again, since its matrices exist only at
// non-negative orders. Worth its own case because everything above runs on a
// grid holding all orders and merely declining to use half of them.
TEST(MatrixKernel, ForwardAgreesOnAScalarGrid) {
  using Grid = GaussLegendreGrid<double, NonNegative, All>;
  CheckKernelsAgreeForward<Grid, double>(8, 8, 0, 1, 1e-13);
  CheckKernelsAgreeForward<Grid, double>(11, 11, 0, 3, 1e-13);
}

TEST(MatrixKernel, RefusesWhatItCannotDo) {
  using Grid = GaussLegendreGrid<double, All, All>;
  using WideGrid = GaussLegendreGrid<long double, All, All>;

  // Generated values cannot be had one order at a time.
  EXPECT_THROW(Grid(8, 2, FFTWpp::Estimate, Chunking::Automatic(),
                    WignerValues::Generated(), TransformKernel::Matrix()),
               std::invalid_argument);

  // BLAS has no long double.
  EXPECT_THROW(WideGrid(8, 2, FFTWpp::Estimate, Chunking::Automatic(),
                        WignerValues::Stored(), TransformKernel::Matrix()),
               std::invalid_argument);

  // Both of those are fine on the loop kernel, which is the point of keeping
  // it: nothing the library could do before is withdrawn.
  EXPECT_NO_THROW(Grid(8, 2, FFTWpp::Estimate, Chunking::Automatic(),
                       WignerValues::Generated(), TransformKernel::Loop()));
  EXPECT_NO_THROW(WideGrid(8, 2, FFTWpp::Estimate, Chunking::Automatic(),
                           WignerValues::Stored(), TransformKernel::Loop()));
}
#endif  // GSHTRANS_HAVE_BLAS

#ifdef GSHTRANS_HAVE_BLAS
namespace {

template <typename Grid, typename Scalar>
void CheckKernelsAgreeInverse(std::ptrdiff_t lMax, std::ptrdiff_t gridDegree,
                              std::ptrdiff_t n, std::ptrdiff_t count,
                              double tolerance) {
  auto loop = Grid(gridDegree, std::abs(n), FFTWpp::Estimate);
  auto matrix = Grid(gridDegree, std::abs(n), FFTWpp::Estimate,
                     Chunking::Automatic(), WignerValues::Stored(),
                     TransformKernel::Matrix());

  const auto fieldSize = static_cast<std::ptrdiff_t>(loop.FieldSize());
  const auto coefficientSize = static_cast<std::ptrdiff_t>(
      std::is_same_v<Scalar, double> ? loop.RealCoefficientSize(lMax)
                                     : loop.CoefficientSize(lMax, n));

  // Coefficients rather than a field, so that the inverse is exercised on its
  // own rather than only as the right inverse of the forward.
  auto coefficients = std::vector<std::complex<double>>(count * coefficientSize);
  for (std::size_t i = 0; i < coefficients.size(); i++) {
    coefficients[i] = std::complex<double>{
        std::cos(0.23 * static_cast<double>(i)),
        std::sin(0.41 * static_cast<double>(i))};
  }

  auto fromLoop = std::vector<Scalar>(count * fieldSize);
  auto fromMatrix = std::vector<Scalar>(count * fieldSize);
  const auto coeffBatch = Batch::Contiguous(count, coefficientSize);
  const auto fieldBatch = Batch::Contiguous(count, fieldSize);

  loop.InverseTransformation(lMax, n, coefficients, coeffBatch, fromLoop,
                             fieldBatch);
  matrix.InverseTransformation(lMax, n, coefficients, coeffBatch, fromMatrix,
                               fieldBatch);

  double scale = 0;
  for (const auto& v : fromLoop) scale = std::max(scale, std::abs(v));
  ASSERT_GT(scale, 0.0);

  for (std::size_t i = 0; i < fromLoop.size(); i++) {
    if constexpr (std::is_same_v<Scalar, double>) {
      EXPECT_NEAR(fromMatrix[i], fromLoop[i], tolerance * scale) << "at " << i;
    } else {
      EXPECT_NEAR(fromMatrix[i].real(), fromLoop[i].real(), tolerance * scale)
          << "at " << i;
      EXPECT_NEAR(fromMatrix[i].imag(), fromLoop[i].imag(), tolerance * scale)
          << "at " << i;
    }
  }
}

}  // namespace

TEST(MatrixKernel, InverseAgreesWithTheLoopKernel) {
  using Grid = GaussLegendreGrid<double, All, All>;
  CheckKernelsAgreeInverse<Grid, std::complex<double>>(8, 8, 0, 1, 1e-13);
  CheckKernelsAgreeInverse<Grid, std::complex<double>>(8, 8, 2, 1, 1e-13);
  CheckKernelsAgreeInverse<Grid, std::complex<double>>(8, 8, -2, 1, 1e-13);
}

TEST(MatrixKernel, InverseAgreesOverABatch) {
  using Grid = GaussLegendreGrid<double, All, All>;
  CheckKernelsAgreeInverse<Grid, std::complex<double>>(10, 10, 1, 5, 1e-13);
  CheckKernelsAgreeInverse<Grid, std::complex<double>>(16, 16, 2, 8, 1e-13);
}

// Below the grid's degree the orders between lMax and nPhi - lMax carry no
// coefficient at all and must still be zero when the FFT reads them. The
// forward direction has no such band, so this is the inverse's own hazard.
TEST(MatrixKernel, InverseAgreesBelowTheGridDegree) {
  using Grid = GaussLegendreGrid<double, All, All>;
  CheckKernelsAgreeInverse<Grid, std::complex<double>>(6, 12, 1, 3, 1e-13);
  CheckKernelsAgreeInverse<Grid, std::complex<double>>(2, 16, 0, 2, 1e-13);
  CheckKernelsAgreeInverse<Grid, std::complex<double>>(1, 20, 1, 1, 1e-13);
}

TEST(MatrixKernel, InverseAgreesForARealField) {
  using Grid = GaussLegendreGrid<double, All, All>;
  CheckKernelsAgreeInverse<Grid, double>(8, 8, 0, 1, 1e-13);
  CheckKernelsAgreeInverse<Grid, double>(12, 12, 0, 4, 1e-13);
  CheckKernelsAgreeInverse<Grid, double>(5, 14, 0, 2, 1e-13);
}

TEST(MatrixKernel, InverseAgreesOnAScalarGrid) {
  using Grid = GaussLegendreGrid<double, NonNegative, All>;
  CheckKernelsAgreeInverse<Grid, double>(8, 8, 0, 1, 1e-13);
  CheckKernelsAgreeInverse<Grid, double>(9, 13, 0, 3, 1e-13);
}

// The round trip on a matrix grid, which the cross-kernel tests cannot see:
// they would both agree on a wrong answer if the two kernels were wrong the
// same way, and transforming a band-limited field there and back is an
// absolute check rather than a relative one.
TEST(MatrixKernel, RoundTripsOnItsOwn) {
  using Grid = GaussLegendreGrid<double, All, All>;
  constexpr auto lMax = std::ptrdiff_t{12};
  constexpr auto n = std::ptrdiff_t{1};
  constexpr auto count = std::ptrdiff_t{3};

  auto grid = Grid(lMax, std::abs(n), FFTWpp::Estimate, Chunking::Automatic(),
                   WignerValues::Stored(), TransformKernel::Matrix());
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(grid.CoefficientSize(lMax, n));

  auto coefficients = std::vector<std::complex<double>>(count * coefficientSize);
  for (std::size_t i = 0; i < coefficients.size(); i++) {
    coefficients[i] = std::complex<double>{
        std::cos(0.19 * static_cast<double>(i)),
        std::sin(0.53 * static_cast<double>(i))};
  }

  auto field = std::vector<std::complex<double>>(count * fieldSize);
  auto back = std::vector<std::complex<double>>(count * coefficientSize);
  const auto cb = Batch::Contiguous(count, coefficientSize);
  const auto fb = Batch::Contiguous(count, fieldSize);

  grid.InverseTransformation(lMax, n, coefficients, cb, field, fb);
  grid.ForwardTransformation(lMax, n, field, fb, back, cb);

  for (std::size_t i = 0; i < coefficients.size(); i++) {
    EXPECT_NEAR(back[i].real(), coefficients[i].real(), 1e-12) << "at " << i;
    EXPECT_NEAR(back[i].imag(), coefficients[i].imag(), 1e-12) << "at " << i;
  }
}
#endif  // GSHTRANS_HAVE_BLAS

#ifdef GSHTRANS_HAVE_BLAS
// Threading over orders (core-plan.md section 11, step M4).
//
// The orders write disjoint output, so a threaded run must give the *same*
// answer as a sequential one -- not a close one. Bit-exact is the right
// standard here even though the cross-kernel comparison is not: within one
// kernel nothing about the arithmetic changes, only which thread does it, and
// each order's product is a single BLAS call either way.
TEST(MatrixKernel, ThreadingChangesNothing) {
  using Grid = GaussLegendreGrid<double, All, All>;
  constexpr auto lMax = std::ptrdiff_t{20};
  constexpr auto n = std::ptrdiff_t{2};
  constexpr auto count = std::ptrdiff_t{4};

  auto grid = Grid(lMax, n, FFTWpp::Estimate, Chunking::Automatic(),
                   WignerValues::Stored(), TransformKernel::Matrix());
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(grid.CoefficientSize(lMax, n));

  auto fields = std::vector<std::complex<double>>(count * fieldSize);
  for (std::size_t i = 0; i < fields.size(); i++) {
    fields[i] = std::complex<double>{std::cos(0.31 * static_cast<double>(i)),
                                     std::sin(0.17 * static_cast<double>(i))};
  }
  const auto fb = Batch::Contiguous(count, fieldSize);
  const auto cb = Batch::Contiguous(count, coefficientSize);

  auto sequential = std::vector<std::complex<double>>(count * coefficientSize);
  auto threaded = std::vector<std::complex<double>>(count * coefficientSize);
  grid.ForwardTransformation(lMax, n, fields, fb, sequential, cb,
                             Execution::Sequential());
  grid.ForwardTransformation(lMax, n, fields, fb, threaded, cb,
                             Execution::Parallel(4));
  for (std::size_t i = 0; i < sequential.size(); i++) {
    ASSERT_EQ(threaded[i], sequential[i]) << "forward, at " << i;
  }

  auto backSequential = std::vector<std::complex<double>>(count * fieldSize);
  auto backThreaded = std::vector<std::complex<double>>(count * fieldSize);
  grid.InverseTransformation(lMax, n, sequential, cb, backSequential, fb,
                             Execution::Sequential());
  grid.InverseTransformation(lMax, n, sequential, cb, backThreaded, fb,
                             Execution::Parallel(4));
  for (std::size_t i = 0; i < backSequential.size(); i++) {
    ASSERT_EQ(backThreaded[i], backSequential[i]) << "inverse, at " << i;
  }
}

// Exactly one level threads: a transform asked to run in parallel from inside
// an existing parallel region runs sequentially instead. The matrix kernel
// has to keep that rule like everything else, and it is worth a test because
// its region is a new one.
TEST(MatrixKernel, DoesNotNestItsThreading) {
  using Grid = GaussLegendreGrid<double, All, All>;
  constexpr auto lMax = std::ptrdiff_t{12};
  constexpr auto n = std::ptrdiff_t{1};

  auto grid = Grid(lMax, n, FFTWpp::Estimate, Chunking::Automatic(),
                   WignerValues::Stored(), TransformKernel::Matrix());
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(grid.CoefficientSize(lMax, n));

  constexpr auto slices = std::ptrdiff_t{4};
  auto fields = std::vector<std::complex<double>>(slices * fieldSize);
  for (std::size_t i = 0; i < fields.size(); i++) {
    fields[i] = std::complex<double>{std::cos(0.07 * static_cast<double>(i)),
                                     std::sin(0.29 * static_cast<double>(i))};
  }
  auto outer = std::vector<std::complex<double>>(slices * coefficientSize);
  auto inner = std::vector<std::complex<double>>(slices * coefficientSize);

  // A caller threading over slices, each slice transforming "in parallel".
#pragma omp parallel for num_threads(4)
  for (std::ptrdiff_t s = 0; s < slices; s++) {
    auto one = std::span<const std::complex<double>>(
        fields.data() + s * fieldSize, fieldSize);
    auto target = std::span<std::complex<double>>(
        outer.data() + s * coefficientSize, coefficientSize);
    grid.ForwardTransformation(lMax, n, one, target, Execution::Parallel(4));
  }

  for (std::ptrdiff_t s = 0; s < slices; s++) {
    auto one = std::span<const std::complex<double>>(
        fields.data() + s * fieldSize, fieldSize);
    auto target = std::span<std::complex<double>>(
        inner.data() + s * coefficientSize, coefficientSize);
    grid.ForwardTransformation(lMax, n, one, target, Execution::Sequential());
  }

  for (std::size_t i = 0; i < outer.size(); i++) {
    ASSERT_EQ(outer[i], inner[i]) << "at " << i;
  }
}
#endif  // GSHTRANS_HAVE_BLAS
