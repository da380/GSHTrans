#include <gtest/gtest.h>

#include <GSHTrans/All>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>

// Reality relations between canonical components, and the removal of
// real-valued transforms at nonzero upper index (core-plan.md step A).
//
// The relation this file exercises is theory note eq:basiclevel,
//
//     (conj f)^{-N}_{l,-m} = (-1)^{m-N} conj(f^N_{lm}),
//
// which is true at every N and is the engine of phase 4's reality reduction.
// It relates two *different* fields, f and conj(f), and it survives.
//
// What does not survive is the self-relation
//
//     f^N_{l,-m} = (-1)^{m-N} conj(f^N_{lm}),
//
// which holds only when f is its own conjugate, i.e. only at N = 0. That
// self-relation is what the reduced m >= 0 coefficient storage assumes, so
// real-valued transforms now exist only at n = 0. The tests below keep the
// relation at n = 2, move the storage tests to n = 0, and pin the rejection.
//
// Written against raw coefficient buffers rather than the
// CanonicalComponentExpansion classes, which are superseded and are deleted at
// step 7 of the field-algebra plan; this file is meant to outlive them.

namespace {

using namespace GSHTrans;

using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;
using ScalarGrid = GaussLegendreGrid<Real, NonNegative, All>;

constexpr Int lMax = 5;
constexpr Int nSpin = 2;
constexpr Real tolerance = 2.0e-11;

Real Phase(Int exponent) {
  return exponent % 2 == 0 ? static_cast<Real>(1) : static_cast<Real>(-1);
}

void ExpectNear(Complex actual, Complex expected) {
  EXPECT_NEAR(actual.real(), expected.real(), tolerance);
  EXPECT_NEAR(actual.imag(), expected.imag(), tolerance);
}

// A real-valued field on the grid, in both real and complex storage. The
// cos(lMax * phi) term deliberately excites the orders m = +-lMax, which the
// grid cannot separate while nPhi = 2 * lMax.
template <typename AnyGrid>
void MakeRealSamples(const AnyGrid& grid, FFTWpp::vector<Real>& realSamples,
                     FFTWpp::vector<Complex>& complexSamples) {
  auto i = Int{0};
  for (auto [theta, phi] : grid.Points()) {
    const auto value = 1.25 + 0.4 * std::cos(theta) +
                       0.3 * std::sin(theta) * std::cos(phi) -
                       0.2 * std::sin(2.0 * theta) * std::sin(2.0 * phi) +
                       0.15 * std::cos(3.0 * theta) * std::cos(3.0 * phi) +
                       0.45 * (1.0 + 0.3 * std::cos(theta)) *
                           std::cos(lMax * phi);
    realSamples[i] = value;
    complexSamples[i] = Complex{value, 0.0};
    ++i;
  }
}

// Full (all orders) coefficients of a complex-valued field.
template <typename AnyGrid>
auto FullCoefficients(const AnyGrid& grid, Int n,
                      const FFTWpp::vector<Complex>& samples) {
  auto coefficients = FFTWpp::vector<Complex>(grid.CoefficientSize(lMax, n));
  std::ranges::fill(coefficients, Complex{});
  grid.ForwardTransformation(lMax, n, samples, coefficients);
  return coefficients;
}

// Reduced (m >= 0) coefficients of a real-valued field. Upper index zero only.
template <typename AnyGrid>
auto ReducedCoefficients(const AnyGrid& grid,
                         const FFTWpp::vector<Real>& samples) {
  auto coefficients = FFTWpp::vector<Complex>(grid.RealCoefficientSize(lMax));
  std::ranges::fill(coefficients, Complex{});
  grid.ForwardTransformation(lMax, 0, samples, coefficients);
  return coefficients;
}

auto FullIndices(Int n) { return GSHIndices<All>(lMax, lMax, n); }
auto ReducedIndices() { return GSHIndices<NonNegative>(lMax, lMax, 0); }

// Returns the message of the std::invalid_argument thrown by f, or an empty
// string if it did not throw. Used to check *why* a call was rejected: an
// out-of-range upper index throws the same type, so testing the type alone
// would let these tests pass for the wrong reason.
template <typename Callable>
std::string RejectionMessage(Callable&& f) {
  try {
    f();
  } catch (const std::invalid_argument& error) {
    return error.what();
  }
  return {};
}

bool Mentions(const std::string& message, const std::string& fragment) {
  return message.find(fragment) != std::string::npos;
}

// A complex-valued transform needs all orders, so on a scalar grid it is not
// merely rejected at run time, it does not compile. Both directions are
// asserted: without the positive case a typo in the expression would make the
// negative one vacuously true.
template <typename AnyGrid>
concept AdmitsComplexForward =
    requires(AnyGrid grid, FFTWpp::vector<Complex> field,
             FFTWpp::vector<Complex> coefficients) {
      grid.ForwardTransformation(lMax, 0, field, coefficients);
    };

static_assert(AdmitsComplexForward<Grid>);
static_assert(!AdmitsComplexForward<ScalarGrid>);

}  // namespace

// The reduced storage at n = 0 holds exactly the non-negative orders of the
// full storage.
TEST(RealFieldSymmetry, ReducedSpectrumMatchesFullNonNegativeOrders) {
  auto grid = Grid(lMax, nSpin, FFTWpp::Estimate);
  auto realSamples = FFTWpp::vector<Real>(grid.FieldSize());
  auto complexSamples = FFTWpp::vector<Complex>(grid.FieldSize());
  MakeRealSamples(grid, realSamples, complexSamples);

  const auto reduced = ReducedCoefficients(grid, realSamples);
  const auto full = FullCoefficients(grid, 0, complexSamples);
  const auto reducedIndices = ReducedIndices();
  const auto fullIndices = FullIndices(0);

  for (auto [l, m] : reducedIndices.Indices()) {
    if (l == lMax && m == lMax) continue;  // see the Nyquist test below
    ExpectNear(reduced[reducedIndices.Index(l, m)],
               full[fullIndices.Index(l, m)]);
  }

  // The zonal coefficients of a real field are real.
  for (auto l : reducedIndices.Degrees()) {
    EXPECT_NEAR(reduced[reducedIndices.Index(l, 0)].imag(), 0.0, tolerance);
  }
}

// The self-relation f^0_{l,-m} = (-1)^m conj(f^0_{lm}) for a real field. This
// is what makes the reduced m >= 0 storage lossless -- and it is exactly what
// fails at nonzero upper index, which is why that storage is gone there.
TEST(RealFieldSymmetry, RealFieldIsSelfConjugateAtUpperIndexZero) {
  auto grid = Grid(lMax, nSpin, FFTWpp::Estimate);
  auto realSamples = FFTWpp::vector<Real>(grid.FieldSize());
  auto complexSamples = FFTWpp::vector<Complex>(grid.FieldSize());
  MakeRealSamples(grid, realSamples, complexSamples);

  const auto full = FullCoefficients(grid, 0, complexSamples);
  const auto indices = FullIndices(0);

  for (auto [l, m] : indices.Indices()) {
    if (l == lMax && std::abs(m) == lMax) continue;  // unresolvable, see below
    ExpectNear(full[indices.Index(l, -m)],
               Phase(m) * std::conj(full[indices.Index(l, m)]));
  }
}

// The orders m = +-lMax are the same discrete mode while nPhi = 2 * lMax, so
// the complex forward transform zeroes (lMax, lMax) and the real transform
// does not. These expectations record the current aliasing behaviour rather
// than any mathematics, and they are expected to change when core-plan.md
// step D (task T3) makes nPhi large enough to resolve both orders.
TEST(RealFieldSymmetry, UnresolvableOrdersFollowTheCurrentAliasingRule) {
  auto grid = Grid(lMax, nSpin, FFTWpp::Estimate);
  auto realSamples = FFTWpp::vector<Real>(grid.FieldSize());
  auto complexSamples = FFTWpp::vector<Complex>(grid.FieldSize());
  MakeRealSamples(grid, realSamples, complexSamples);

  const auto reduced = ReducedCoefficients(grid, realSamples);
  const auto full = FullCoefficients(grid, 0, complexSamples);
  const auto nyquist = reduced[ReducedIndices().Index(lMax, lMax)];
  const auto indices = FullIndices(0);

  EXPECT_NEAR(nyquist.imag(), 0.0, tolerance);
  EXPECT_GT(std::abs(nyquist), tolerance);
  ExpectNear(full[indices.Index(lMax, lMax)], Complex{});
  ExpectNear(full[indices.Index(lMax, -lMax)], Phase(lMax) * std::conj(nyquist));
}

// eq:basiclevel at n = 2, through complex transforms of a real-valued field:
// the coefficients of conj(f) at -N are determined by those of f at +N. This
// relation is true at every N and is the oracle phase 4 reuses.
TEST(RealFieldSymmetry, FullTransformsSatisfyCrossUpperIndexIdentity) {
  auto grid = Grid(lMax, nSpin, FFTWpp::Estimate);
  auto realSamples = FFTWpp::vector<Real>(grid.FieldSize());
  auto complexSamples = FFTWpp::vector<Complex>(grid.FieldSize());
  MakeRealSamples(grid, realSamples, complexSamples);

  const auto plus = FullCoefficients(grid, nSpin, complexSamples);
  const auto minus = FullCoefficients(grid, -nSpin, complexSamples);
  const auto plusIndices = FullIndices(nSpin);
  const auto minusIndices = FullIndices(-nSpin);

  for (auto [l, m] : plusIndices.Indices()) {
    if (l == lMax && std::abs(m) == lMax) continue;
    ExpectNear(minus[minusIndices.Index(l, -m)],
               Phase(m - nSpin) * std::conj(plus[plusIndices.Index(l, m)]));
  }
}

// The reduced inverse at n = 0 reconstructs the same field as the full inverse
// applied to the Hermitian pair the reduced storage stands for.
TEST(RealFieldSymmetry, ReducedInverseUsesTheImplicitHermitianPair) {
  auto grid = Grid(lMax, nSpin, FFTWpp::Estimate);
  auto realSamples = FFTWpp::vector<Real>(grid.FieldSize());
  auto complexSamples = FFTWpp::vector<Complex>(grid.FieldSize());
  MakeRealSamples(grid, realSamples, complexSamples);

  const auto reduced = ReducedCoefficients(grid, realSamples);
  const auto reducedIndices = ReducedIndices();
  const auto fullIndices = FullIndices(0);

  auto expanded = FFTWpp::vector<Complex>(grid.CoefficientSize(lMax, 0));
  std::ranges::fill(expanded, Complex{});
  for (auto [l, m] : reducedIndices.Indices()) {
    const auto value = reduced[reducedIndices.Index(l, m)];
    expanded[fullIndices.Index(l, m)] = value;
    if (m > 0 && m < lMax) {
      expanded[fullIndices.Index(l, -m)] = Phase(m) * std::conj(value);
    }
  }

  auto reducedField = FFTWpp::vector<Real>(grid.FieldSize());
  auto expandedField = FFTWpp::vector<Complex>(grid.FieldSize());
  grid.InverseTransformation(lMax, 0, reduced, reducedField);
  grid.InverseTransformation(lMax, 0, expanded, expandedField);

  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    EXPECT_NEAR(expandedField[i].imag(), 0.0, tolerance);
    EXPECT_NEAR(reducedField[i], expandedField[i].real(), tolerance);
  }
}

// The removal itself.
TEST(RealFieldSymmetry, RealTransformsAreRejectedAtNonzeroUpperIndex) {
  auto grid = Grid(lMax, nSpin, FFTWpp::Estimate);
  auto realField = FFTWpp::vector<Real>(grid.FieldSize());
  auto coefficients =
      FFTWpp::vector<Complex>(GSHIndices<NonNegative>(lMax, lMax, nSpin).Size());

  for (auto n : {nSpin, -nSpin}) {
    const auto forward = RejectionMessage([&] {
      grid.ForwardTransformation(lMax, n, realField, coefficients);
    });
    const auto inverse = RejectionMessage([&] {
      grid.InverseTransformation(lMax, n, coefficients, realField);
    });
    EXPECT_TRUE(Mentions(forward, "upper index zero")) << forward;
    EXPECT_TRUE(Mentions(inverse, "upper index zero")) << inverse;
  }

  // Control: an upper index outside the grid's range is rejected too, but for
  // a different reason. If these two messages ever coincide, the test above
  // has stopped discriminating.
  const auto outOfRange = RejectionMessage(
      [&] { grid.ForwardTransformation(lMax, lMax, realField, coefficients); });
  EXPECT_FALSE(Mentions(outOfRange, "upper index zero")) << outOfRange;
  EXPECT_FALSE(outOfRange.empty());

  // The complex transforms at the same upper index are untouched.
  auto complexField = FFTWpp::vector<Complex>(grid.FieldSize());
  auto full = FFTWpp::vector<Complex>(grid.CoefficientSize(lMax, nSpin));
  EXPECT_NO_THROW(
      grid.ForwardTransformation(lMax, nSpin, complexField, full));
}

// An MRange = NonNegative grid is a real scalar grid: it serves real transforms
// at upper index zero and nothing else, so a nonzero maximum upper index leaves
// it able to serve nothing (core-plan.md [C1]).
TEST(RealFieldSymmetry, ScalarGridRejectsNonzeroMaximumUpperIndex) {
  for (auto nMax : {1, static_cast<int>(nSpin)}) {
    const auto message =
        RejectionMessage([&] { ScalarGrid(lMax, nMax, FFTWpp::Estimate); });
    EXPECT_TRUE(Mentions(message, "maximum upper index must be zero"))
        << message;
  }
  EXPECT_NO_THROW(ScalarGrid(lMax, 0, FFTWpp::Estimate));
}

// And that a real scalar grid actually works at n = 0.
TEST(RealFieldSymmetry, ScalarGridRoundTripsAtUpperIndexZero) {
  auto grid = ScalarGrid(lMax, 0, FFTWpp::Estimate);
  auto realSamples = FFTWpp::vector<Real>(grid.FieldSize());
  auto complexSamples = FFTWpp::vector<Complex>(grid.FieldSize());
  MakeRealSamples(grid, realSamples, complexSamples);

  const auto coefficients = ReducedCoefficients(grid, realSamples);
  auto reconstructed = FFTWpp::vector<Real>(grid.FieldSize());
  grid.InverseTransformation(lMax, 0, coefficients, reconstructed);

  // The (lMax, lMax) mode is not resolvable at nPhi = 2 * lMax, so the round
  // trip is checked against a field with that mode projected out: forward,
  // inverse, forward again must reproduce the first spectrum.
  auto again = FFTWpp::vector<Complex>(grid.RealCoefficientSize(lMax));
  std::ranges::fill(again, Complex{});
  grid.ForwardTransformation(lMax, 0, reconstructed, again);
  for (auto i = Int{0}; i < static_cast<Int>(coefficients.size()); ++i) {
    ExpectNear(again[i], coefficients[i]);
  }

  // The comparison above is idempotence, which an all-zero spectrum would
  // satisfy. Check the spectrum is actually carrying the field.
  const auto largest = std::ranges::max(
      coefficients | std::ranges::views::transform(
                         [](auto c) { return std::abs(c); }));
  EXPECT_GT(largest, 1.0);
}
