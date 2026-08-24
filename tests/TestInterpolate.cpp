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

}  // namespace
