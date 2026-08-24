#include <gtest/gtest.h>

#include <GSHTrans/All>

#include <complex>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

#include "TestRandom.h"

namespace {

using namespace GSHTrans;

using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

constexpr auto pi = std::numbers::pi_v<Real>;

// A field's worth of arbitrary but reproducible samples.
auto Samples(Int n, int seed) {
  auto values = std::vector<Complex>();
  values.reserve(static_cast<std::size_t>(n));
  for (auto i = Int{0}; i < n; i++) {
    seed = seed * 1103515245 + 12345;
    const auto a = static_cast<Real>((seed >> 8) % 1000) / 1000;
    seed = seed * 1103515245 + 12345;
    const auto b = static_cast<Real>((seed >> 8) % 1000) / 1000;
    values.emplace_back(a, b);
  }
  return values;
}

//--------------------------------------------------------------------------//
//         P1: the padded grid, which needs no expansion to check            //
//--------------------------------------------------------------------------//

class PaddedGrid : public ::testing::Test {
 protected:
  Grid grid{8, 2};
  Int nTheta = static_cast<Int>(grid.NumberOfCoLatitudes());
  Int nPhi = static_cast<Int>(grid.NumberOfLongitudes());

  auto Build() const {
    const auto values = Samples(nTheta * nPhi, 11);
    const auto north = Samples(nPhi, 22);
    const auto south = Samples(nPhi, 33);
    return std::tuple{values, north, south,
                      InterpolateDetails::Pad<Grid, Complex>(
                          grid, values, north, south)};
  }
};

TEST_F(PaddedGrid, CoversTheClosedSphere) {
  const auto [values, north, south, padded] = Build();

  EXPECT_EQ(padded.Rows(), static_cast<std::size_t>(nTheta + 2));
  EXPECT_EQ(padded.Columns(), static_cast<std::size_t>(nPhi + 1));
  EXPECT_EQ(padded.values.size(), padded.Rows() * padded.Columns());

  EXPECT_DOUBLE_EQ(padded.theta.front(), 0.0);
  EXPECT_DOUBLE_EQ(padded.theta.back(), pi);
  EXPECT_DOUBLE_EQ(padded.phi.front(), 0.0);
  EXPECT_DOUBLE_EQ(padded.phi.back(), 2 * pi);
}

// Upstream refuses an axis that is not strictly increasing, so this is the
// property that decides whether the padded grid is usable at all.
TEST_F(PaddedGrid, AxesAreStrictlyIncreasing) {
  const auto [values, north, south, padded] = Build();

  for (std::size_t i = 1; i < padded.theta.size(); i++) {
    EXPECT_GT(padded.theta[i], padded.theta[i - 1]) << "at colatitude " << i;
  }
  for (std::size_t j = 1; j < padded.phi.size(); j++) {
    EXPECT_GT(padded.phi[j], padded.phi[j - 1]) << "at longitude " << j;
  }
}

// The interior is the field, unmoved: the whole point of the layout agreeing
// is that padding is an insertion rather than a repack.
TEST_F(PaddedGrid, ReproducesTheFieldAtEveryOriginalNode) {
  const auto [values, north, south, padded] = Build();

  for (auto i = Int{0}; i < nTheta; i++) {
    EXPECT_DOUBLE_EQ(padded.theta[static_cast<std::size_t>(i + 1)],
                     grid.CoLatitudes()[i]);
    for (auto j = Int{0}; j < nPhi; j++) {
      EXPECT_EQ(padded.At(static_cast<std::size_t>(i + 1),
                          static_cast<std::size_t>(j)),
                values[static_cast<std::size_t>(i * nPhi + j)])
          << "at (" << i << ", " << j << ")";
    }
  }
}

TEST_F(PaddedGrid, CarriesThePolarRowsAsGiven) {
  const auto [values, north, south, padded] = Build();

  for (auto j = Int{0}; j < nPhi; j++) {
    EXPECT_EQ(padded.At(0, static_cast<std::size_t>(j)),
              north[static_cast<std::size_t>(j)]);
    EXPECT_EQ(padded.At(padded.Rows() - 1, static_cast<std::size_t>(j)),
              south[static_cast<std::size_t>(j)]);
  }
}

// phi = 2 pi is phi = 0. Every row closes on itself, the polar rows included,
// which is what makes the wrap exact rather than an approximation.
TEST_F(PaddedGrid, WrapColumnIsColumnZero) {
  const auto [values, north, south, padded] = Build();

  for (std::size_t i = 0; i < padded.Rows(); i++) {
    EXPECT_EQ(padded.At(i, padded.Columns() - 1), padded.At(i, 0))
        << "at colatitude " << i;
  }
}

TEST_F(PaddedGrid, RefusesInputsThatDoNotFitTheGrid) {
  const auto values = Samples(nTheta * nPhi, 11);
  const auto north = Samples(nPhi, 22);
  const auto shortRow = Samples(nPhi - 1, 44);
  const auto shortField = Samples(nTheta * nPhi - 1, 55);

  EXPECT_THROW((InterpolateDetails::Pad<Grid, Complex>(grid, shortField, north,
                                                       north)),
               std::invalid_argument);
  EXPECT_THROW(
      (InterpolateDetails::Pad<Grid, Complex>(grid, values, shortRow, north)),
      std::invalid_argument);
  EXPECT_THROW(
      (InterpolateDetails::Pad<Grid, Complex>(grid, values, north, shortRow)),
      std::invalid_argument);
}

// A real field pads exactly as a complex one does; the scalar type is carried
// through rather than promoted, which is what [I9] asks for.
TEST_F(PaddedGrid, PadsARealFieldWithoutPromoting) {
  auto values = std::vector<Real>(
      static_cast<std::size_t>(nTheta * nPhi), Real{2});
  auto row = std::vector<Real>(static_cast<std::size_t>(nPhi), Real{5});

  const auto padded =
      InterpolateDetails::Pad<Grid, Real>(grid, values, row, row);

  static_assert(std::same_as<decltype(padded.values)::value_type, Real>);
  EXPECT_DOUBLE_EQ(padded.At(1, 0), 2.0);
  EXPECT_DOUBLE_EQ(padded.At(0, 0), 5.0);
}

//--------------------------------------------------------------------------//
//              P2: the spectral interpolant, which is the reference          //
//--------------------------------------------------------------------------//

// Fill an expansion with arbitrary but reproducible coefficients.
template <typename Expansion>
void Fill(Expansion& e, int seed) {
  for (auto l : e.Degrees()) {
    for (auto m : e.Orders(l)) {
      seed = seed * 1103515245 + 12345;
      const auto a = static_cast<Real>((seed >> 8) % 1000) / 1000;
      seed = seed * 1103515245 + 12345;
      const auto b = static_cast<Real>((seed >> 8) % 1000) / 1000;
      // A real field's m = 0 coefficient is real, by its own reality
      // condition; giving it an imaginary part would make the expansion
      // describe no field at all.
      if constexpr (std::same_as<typename Expansion::Value, RealValued>) {
        e[l, m] = m == 0 ? Complex(a, 0) : Complex(a, b);
      } else {
        e[l, m] = Complex(a, b);
      }
    }
  }
}

// The interpolant must be copyable, because ProjectFunction takes its
// callable by value and copies it into a lambda -- which is [I8]'s whole
// point, and the property most easily broken by a change of storage.
static_assert(std::copy_constructible<SpectralInterpolant<2, Grid>>);
static_assert(std::copy_constructible<
              SpectralInterpolant<0, Grid, RealValued>>);

// The decisive agreement: the same numbers as the transform, at every point
// the transform produces. This is what says the direct sum of section 22.1 is
// the synthesis rather than something like it.
TEST(SpectralInterpolant, MatchesEvaluateAtEveryGridPoint) {
  constexpr Int N = 2;
  const Int lMax = 8;
  auto grid = Grid(lMax, 2);
  auto e = SpinExpansion<N, Grid>(grid, lMax);
  Fill(e, 101);

  const auto field = Evaluate(e);
  const auto at = Interpolate(e);

  auto worst = Real{0};
  auto iTheta = Int{0};
  for (auto theta : grid.CoLatitudes()) {
    auto iPhi = Int{0};
    for (auto phi : grid.Longitudes()) {
      worst = std::max(worst, std::abs(at(theta, phi) - field[iTheta, iPhi]));
      iPhi++;
    }
    iTheta++;
  }
  EXPECT_LT(worst, 1e-13) << "worst difference " << worst;
}

// The same, through the reduced m >= 0 storage, where the sum has to supply
// the Hermitian partner f_{l,-m} = (-1)^m conj(f_{lm}) itself.
TEST(SpectralInterpolant, MatchesEvaluateForARealField) {
  const Int lMax = 8;
  auto grid = Grid(lMax, 0);
  auto e = SpinExpansion<0, Grid, RealValued>(grid, lMax);
  Fill(e, 202);

  const auto field = Evaluate(e);
  const auto at = Interpolate(e);

  static_assert(std::same_as<decltype(at(0.0, 0.0)), Real>);

  auto worst = Real{0};
  auto iTheta = Int{0};
  for (auto theta : grid.CoLatitudes()) {
    auto iPhi = Int{0};
    for (auto phi : grid.Longitudes()) {
      worst = std::max(worst, std::abs(at(theta, phi) - field[iTheta, iPhi]));
      iPhi++;
    }
    iTheta++;
  }
  EXPECT_LT(worst, 1e-13) << "worst difference " << worst;
}

// Off the grid, where agreement with Evaluate says nothing, the oracle is a
// closed form. The l = 1 generalised Legendre functions of Dahlen & Tromp
// (C.115) are the table field-algebra-plan.md section 2 uses to pin the
// convention, so this is the same oracle at a different point.
TEST(SpectralInterpolant, IsExactOnALowDegreeHarmonic) {
  constexpr Int N = 1;
  auto grid = Grid(1, 1);
  auto e = SpinExpansion<N, Grid>(grid, 1);
  const auto c = std::array{Complex(0.3, -0.7), Complex(-1.1, 0.4),
                            Complex(0.9, 0.2)};
  e[1, -1] = c[0];
  e[1, 0] = c[1];
  e[1, 1] = c[2];

  const auto at = Interpolate(e);
  const auto norm = std::sqrt(3 / (4 * pi));

  for (auto theta : {0.31, 1.0, 2.4, 3.0}) {
    for (auto phi : {0.0, 0.9, 4.7}) {
      // P^{1}_{1,-1} = (1 - cos)/2,  P^{1}_{1,0} = sin/sqrt(2),
      // P^{1}_{1,+1} = (1 + cos)/2.
      const auto want =
          norm * (c[0] * ((1 - std::cos(theta)) / 2) *
                      std::exp(Complex(0, -phi)) +
                  c[1] * (std::sin(theta) / std::sqrt(2.0)) +
                  c[2] * ((1 + std::cos(theta)) / 2) *
                      std::exp(Complex(0, phi)));
      EXPECT_NEAR(std::abs(at(theta, phi) - want), 0.0, 1e-14)
          << "at theta = " << theta << ", phi = " << phi;
    }
  }
}

// The poles are inside the domain and are where the whole padding question
// comes from, so the reference must answer there. Section 22.1's rule:
// the order m = +N survives at the north and m = -N at the south, the latter
// with a sign alternating in the degree.
TEST(SpectralInterpolant, AnswersAtThePolesByTheStatedRule) {
  constexpr Int N = 2;
  const Int lMax = 6;
  auto grid = Grid(lMax, 2);
  auto e = SpinExpansion<N, Grid>(grid, lMax);
  Fill(e, 303);

  const auto at = Interpolate(e);

  auto north = Complex{};
  auto south = Complex{};
  for (auto l = N; l <= lMax; l++) {
    const auto norm = std::sqrt((2 * static_cast<Real>(l) + 1) / (4 * pi));
    north += e[l, N] * norm;
    south += e[l, -N] * norm * ((l - N) % 2 == 0 ? 1.0 : -1.0);
  }

  for (auto phi : {0.0, 1.3, 5.5}) {
    EXPECT_NEAR(std::abs(at(0.0, phi) - north * std::exp(Complex(0, N * phi))),
                0.0, 1e-13)
        << "north pole at phi = " << phi;
    EXPECT_NEAR(
        std::abs(at(pi, phi) - south * std::exp(Complex(0, -N * phi))), 0.0,
        1e-13)
        << "south pole at phi = " << phi;
  }
}

// A pole value depends on phi at N != 0, which is the frame ambiguity rather
// than a defect, and is the fact section 9 of thoughts.md got wrong. Pinned
// because a polar row built as a constant would pass every other test here.
TEST(SpectralInterpolant, PoleValueVariesWithLongitudeAtNonzeroUpperIndex) {
  constexpr Int N = 2;
  auto grid = Grid(6, 2);
  auto e = SpinExpansion<N, Grid>(grid, 6);
  Fill(e, 404);
  const auto at = Interpolate(e);

  const auto a = at(0.0, 0.0);
  const auto b = at(0.0, pi / (2 * N));  // a quarter turn of exp(i N phi)
  EXPECT_GT(std::abs(a - b), 1e-3) << "the north pole looks constant in phi";
  EXPECT_NEAR(std::abs(a), std::abs(b), 1e-13) << "only the phase should move";
}

//--------------------------------------------------------------------------//
//        P3: the polar rows, and Interpolate on a field                      //
//--------------------------------------------------------------------------//

// Every scheme here passes through its own nodes, so this is the identity
// test -- the same check RadialResample's turned on, and the one that says
// the padding's indices line up with the field's.
template <typename Field, typename Interpolant>
void ExpectReproducesNodes(const Field& field, const Interpolant& at,
                           Real tolerance) {
  const auto& grid = field.Grid();
  auto worst = Real{0};
  auto iTheta = Int{0};
  for (auto theta : grid.CoLatitudes()) {
    auto iPhi = Int{0};
    for (auto phi : grid.Longitudes()) {
      worst = std::max(worst, std::abs(at(theta, phi) - field[iTheta, iPhi]));
      iPhi++;
    }
    iTheta++;
  }
  EXPECT_LT(worst, tolerance) << "worst difference at a node " << worst;
}

TEST(FieldInterpolant, EverySchemeReproducesTheFieldAtEveryNode) {
  constexpr Int N = 2;
  const Int lMax = 8;
  auto grid = Grid(lMax, 2);
  auto e = SpinExpansion<N, Grid>(grid, lMax);
  Fill(e, 505);
  const auto field = Evaluate(e);

  ExpectReproducesNodes(field, Interpolate(field, Scheme::Spectral()), 1e-12);
#ifdef GSHTRANS_HAVE_INTERPOLATION
  ExpectReproducesNodes(field, Interpolate(field, Scheme::Bilinear()), 1e-13);
  ExpectReproducesNodes(field, Interpolate(field, Scheme::Bicubic()), 1e-12);
#endif
}

TEST(FieldInterpolant, EverySchemeReproducesARealFieldAtEveryNode) {
  const Int lMax = 8;
  auto grid = Grid(lMax, 0);
  auto e = SpinExpansion<0, Grid, RealValued>(grid, lMax);
  Fill(e, 606);
  const auto field = Evaluate(e);

  ExpectReproducesNodes(field, Interpolate(field, Scheme::Spectral()), 1e-12);
#ifdef GSHTRANS_HAVE_INTERPOLATION
  ExpectReproducesNodes(field, Interpolate(field, Scheme::Bilinear()), 1e-13);
  ExpectReproducesNodes(field, Interpolate(field, Scheme::Bicubic()), 1e-12);
#endif
}

#ifdef GSHTRANS_HAVE_INTERPOLATION

// The poles are the reason the padding exists, and a polar row is exact where
// it is a node -- which is at the grid's own longitudes and nowhere else.
// Between them the row is interpolated along phi like any other row, because
// exp(i N phi) is not what a local scheme reproduces. Measured: 2e-15 at the
// nodes against 6e-3 (bicubic) and 1e-1 (bilinear) between them, the latter
// being no worse than the same schemes manage in the interior.
//
// So this asserts the exactness where it is claimed and not where it is not,
// which is the distinction section 22's P5 exists to measure.
TEST(FieldInterpolant, LocalSchemesAreExactAtThePolarNodes) {
  constexpr Int N = 2;
  const Int lMax = 8;
  auto grid = Grid(lMax, 2);
  auto e = SpinExpansion<N, Grid>(grid, lMax);
  Fill(e, 707);
  const auto field = Evaluate(e);

  const auto reference = Interpolate(field, Scheme::Spectral());
  const auto bilinear = Interpolate(field, Scheme::Bilinear());
  const auto bicubic = Interpolate(field, Scheme::Bicubic());

  for (auto phi : grid.Longitudes()) {
    for (auto theta : {0.0, pi}) {
      EXPECT_NEAR(std::abs(bilinear(theta, phi) - reference(theta, phi)), 0.0,
                  1e-13)
          << "bilinear at theta = " << theta << ", phi = " << phi;
      EXPECT_NEAR(std::abs(bicubic(theta, phi) - reference(theta, phi)), 0.0,
                  1e-13)
          << "bicubic at theta = " << theta << ", phi = " << phi;
    }
  }
}

// A constant polar row would satisfy every node test above and be wrong at
// every N != 0, so the phase is pinned separately -- on a local scheme, since
// that is where a constant row would have been written. At the row's own
// nodes the ratio between two longitudes is exactly exp(i N (phi - phi')).
TEST(FieldInterpolant, PolarRowsCarryTheFramePhase) {
  constexpr Int N = 2;
  auto grid = Grid(8, 2);
  auto e = SpinExpansion<N, Grid>(grid, 8);
  Fill(e, 808);
  const auto field = Evaluate(e);
  const auto at = Interpolate(field, Scheme::Bicubic());

  const auto base = at(0.0, 0.0);
  ASSERT_GT(std::abs(base), 1e-6) << "degenerate test data";

  auto moved = false;
  for (auto phi : grid.Longitudes()) {
    const auto want = base * std::exp(Complex(0, N * phi));
    EXPECT_NEAR(std::abs(at(0.0, phi) - want), 0.0, 1e-13)
        << "north pole at phi = " << phi;
    if (std::abs(at(0.0, phi) - base) > 1e-3) moved = true;
  }
  EXPECT_TRUE(moved) << "the north pole looks constant in phi";
}

// The wrap: the last cell used to interpolate against nothing. Just below
// 2 pi the answer must approach the value at zero, which it cannot do without
// the extra column.
TEST(FieldInterpolant, TheLastLongitudeCellClosesOnZero) {
  const Int lMax = 8;
  auto grid = Grid(lMax, 0);
  auto e = SpinExpansion<0, Grid, RealValued>(grid, lMax);
  Fill(e, 909);
  const auto field = Evaluate(e);
  const auto at = Interpolate(field, Scheme::Bilinear());

  EXPECT_NEAR(at(1.0, 2 * pi - 1e-9), at(1.0, 0.0), 1e-7);
  // And a query past 2 pi is the same point, reduced. Not bit-identical:
  // fmod(2 pi + 0.3, 2 pi) is 0.3 in exact arithmetic and one ulp away in
  // floating point, so the two land in the same cell at slightly different
  // places. That is a property of the reduction, not of the interpolant.
  EXPECT_NEAR(at(1.0, 2 * pi + 0.3), at(1.0, 0.3), 1e-12);
}

// Interpolate is a template over the scheme tag, so a build without the
// dependency does not have these two factories at all. Guarded here so the
// test file compiles either way, which is what says the guard is real.
static_assert(requires { Scheme::Bilinear(); });
static_assert(requires { Scheme::Bicubic(); });

#endif  // GSHTRANS_HAVE_INTERPOLATION

//--------------------------------------------------------------------------//
//     P4: the interpolant as a function on the sphere, and remeshing         //
//--------------------------------------------------------------------------//

// [I8]: modelling ScalarFunctionS2 is the point of the feature rather than a
// bonus, because it is what makes remeshing one line. Asserted because it is
// the property most easily broken by a change of signature.
static_assert(
    ScalarFunctionS2<SpectralInterpolant<2, Grid>, Real, Complex>);
static_assert(ScalarFunctionS2<SpectralInterpolant<0, Grid, RealValued>, Real,
                               Real>);

// Remeshing, end to end: a band-limited field sampled on one grid, evaluated
// on another. Spectral interpolation is exact for such a field, so the
// remeshed samples must agree with the expansion evaluated on the second grid
// directly -- which is an independent route to the same numbers.
TEST(FieldInterpolant, RemeshesOntoAnotherGridExactly) {
  constexpr Int N = 2;
  const Int band = 8;
  auto coarse = Grid(band, 2);
  auto fine = Grid(band + 4, 2);

  auto e = SpinExpansion<N, Grid>(coarse, band);
  Fill(e, 1111);
  const auto field = Evaluate(e);

  // The same coefficients on the finer grid, evaluated there.
  auto eFine = SpinExpansion<N, Grid>(fine, band);
  for (auto l : e.Degrees())
    for (auto m : e.Orders(l)) eFine[l, m] = e[l, m];
  const auto want = Evaluate(eFine);

  // And by handing the interpolant to the field constructor, which is the
  // one-line remesh [I8] promises.
  const auto got = SpinField<N, Grid>(fine, Interpolate(field));

  auto worst = Real{0};
  for (auto iTheta : fine.CoLatitudeIndices())
    for (auto iPhi : fine.LongitudeIndices())
      worst = std::max(worst,
                       std::abs(got[iTheta, iPhi] - want[iTheta, iPhi]));
  EXPECT_LT(worst, 1e-12) << "worst difference " << worst;
}

// ProjectFunction takes its callable by value and copies it into a lambda, so
// this exercises the copyability that [I1]'s shared state exists to provide.
TEST(FieldInterpolant, SurvivesProjectFunctionWhichCopiesIt) {
  const Int lMax = 6;
  auto grid = Grid(lMax, 0);
  auto e = SpinExpansion<0, Grid, RealValued>(grid, lMax);
  Fill(e, 1212);
  const auto field = Evaluate(e);

  const auto at = Interpolate(field);
  auto worst = Real{0};
  auto iPoint = Int{0};
  for (auto value : grid.ProjectFunction(at)) {
    worst = std::max(worst, std::abs(value - field.Data()[iPoint]));
    iPoint++;
  }
  EXPECT_EQ(iPoint, grid.FieldSize());
  EXPECT_LT(worst, 1e-12) << "worst difference " << worst;
}

TEST(FieldInterpolant, RefusesAColatitudeOffTheSphere) {
  auto grid = Grid(4, 0);
  auto e = SpinExpansion<0, Grid, RealValued>(grid, 4);
  Fill(e, 1010);
  const auto at = Interpolate(e);

  EXPECT_THROW(at(-0.1, 0.0), std::invalid_argument);
  EXPECT_THROW(at(pi + 0.1, 0.0), std::invalid_argument);
  EXPECT_NO_THROW(at(0.0, 0.0));
  EXPECT_NO_THROW(at(pi, 0.0));
}

}  // namespace
