#include <gtest/gtest.h>

#include <GSHTrans/All>

#include <cmath>
#include <complex>
#include <cstddef>
#include <numeric>
#include <span>
#include <vector>

namespace {

using namespace GSHTrans;
using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

auto Radii(Int nR, Real inner = 0.5, Real outer = 1.0) {
  auto r = std::vector<Real>(nR);
  for (Int i = 0; i < nR; i++) {
    r[i] = inner + (outer - inner) * i / static_cast<Real>(nR - 1);
  }
  return RadialGrid<Real>(r);
}

}  // namespace

//--------------------------------------------------------------------------//
//                              The radial grid                              //
//--------------------------------------------------------------------------//

TEST(RadialGrid, CarriesNodesWeightsAndIdentity) {
  auto r = std::vector<Real>{1.0, 2.0, 3.0};
  auto w = std::vector<Real>{0.5, 1.0, 0.5};
  auto grid = RadialGrid<Real>(r, w);

  EXPECT_EQ(grid.NumberOfRadii(), 3);
  EXPECT_EQ(grid.Radius(1), 2.0);
  EXPECT_EQ(grid.Weight(1), 1.0);
  EXPECT_TRUE(grid.HasWeights());

  // Points alone are allowed: a caller who never integrates radially need not
  // invent weights.
  auto bare = RadialGrid<Real>(r);
  EXPECT_FALSE(bare.HasWeights());

  // Identity, not structure, exactly as for the angular grid.
  auto copy = grid;
  auto twin = RadialGrid<Real>(r, w);
  EXPECT_EQ(copy.Identity(), grid.Identity());
  EXPECT_NE(twin.Identity(), grid.Identity());
}

TEST(RadialGrid, RejectsWhatItCannotHold) {
  EXPECT_THROW(RadialGrid<Real>(std::vector<Real>{}), std::invalid_argument);
  EXPECT_THROW((RadialGrid<Real>({1.0, 2.0}, {1.0})), std::invalid_argument);
  EXPECT_THROW((RadialGrid<Real>({2.0, 1.0})), std::invalid_argument);
  EXPECT_THROW((RadialGrid<Real>({-1.0, 1.0})), std::invalid_argument);
}

//--------------------------------------------------------------------------//
//                            The layered field                              //
//--------------------------------------------------------------------------//

TEST(LayeredSpinField, IsAStackOfAngularSlices) {
  constexpr auto nR = Int{5};
  auto grid = Grid(8, 2, FFTWpp::Estimate);
  auto radial = Radii(nR);

  auto f = LayeredSpinField<1, Grid>(radial, grid);
  EXPECT_EQ(f.NumberOfRadii(), nR);
  EXPECT_EQ(f.Size(), nR * f.FieldSize());

  // A slice is an ordinary phase-1 node, not a new kind of object.
  auto slice = f.Slice(2);
  static_assert(SpinWeighted<decltype(slice)>);
  static_assert(decltype(slice)::UpperIndex == 1);

  // Writing through it writes the stack, radius-major.
  slice[0, 0] = Complex{7.0, -1.0};
  EXPECT_EQ(f.Data()[2 * f.FieldSize()], (Complex{7.0, -1.0}));

  // And the slices are disjoint.
  EXPECT_EQ((f.Slice(1)[0, 0]), Complex{});
  EXPECT_EQ((f.Slice(3)[0, 0]), Complex{});

  EXPECT_THROW(f.Slice(nR), std::invalid_argument);
  EXPECT_THROW(f.Slice(-1), std::invalid_argument);
}

TEST(LayeredSpinField, SlicesComposeWithTheAngularAlgebra) {
  auto grid = Grid(8, 2, FFTWpp::Estimate);
  auto radial = Radii(4);

  auto f = LayeredSpinField<2, Grid>(radial, grid);
  auto g = LayeredSpinField<2, Grid>(radial, grid);
  for (auto i : f.RadiusIndices()) {
    auto a = f.Slice(i);
    auto b = g.Slice(i);
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        a[iTheta, iPhi] = Complex{std::cos(0.3 * iTheta + i), 0.2 * iPhi};
        b[iTheta, iPhi] = Complex{std::sin(0.4 * iPhi), 0.1 * i};
      }
    }
  }

  // The whole index algebra lifts: conj reverses, the product lands at zero,
  // and only there can it be integrated. Nothing about the stack is involved.
  const auto& stack = f;
  const auto& other = g;
  auto pairing = Integrate(conj(stack.Slice(1)) * other.Slice(1));
  static_assert(std::same_as<decltype(pairing), Complex>);

  auto expected = Integrate(conj(f.Slice(1)) * g.Slice(1));
  EXPECT_EQ(pairing, expected);
}

TEST(LayeredSpinField, BroadcastLiftsAnExpressionToEveryRadius) {
  auto grid = Grid(8, 2, FFTWpp::Estimate);
  auto radial = Radii(4);

  auto u = SpinField<1, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::cos(theta), std::sin(phi)};
  });

  // Any spin-weighted node, so a lazy expression needs no materialising.
  const auto& node = u;
  auto stack = Broadcast(radial, 2.0 * node);
  static_assert(decltype(stack)::UpperIndex == 1);
  EXPECT_EQ(stack.NumberOfRadii(), 4);

  for (auto i : stack.RadiusIndices()) {
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        EXPECT_EQ((stack.Slice(i)[iTheta, iPhi]),
                  2.0 * (u[iTheta, iPhi]));
      }
    }
  }
}

//--------------------------------------------------------------------------//
//                     The radial axis is the batch axis                     //
//--------------------------------------------------------------------------//

// The whole stack transforms in one call, and the answer must be exactly what
// transforming each radius separately gives -- batching widens the inner loop
// without reordering any sum.
TEST(LayeredSpinField, TransformsAsOneBatchAndAgreesExactly) {
  constexpr auto lMax = Int{6};
  constexpr auto nR = Int{7};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = Radii(nR);

  auto f = LayeredSpinField<2, Grid>(radial, grid);
  for (auto i : f.RadiusIndices()) {
    auto slice = f.Slice(i);
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        slice[iTheta, iPhi] =
            Complex{std::cos(0.31 * iTheta + i), std::sin(0.17 * iPhi - i)};
      }
    }
  }

  auto stacked = Expand(f, lMax);
  EXPECT_EQ(stacked.NumberOfRadii(), nR);
  EXPECT_EQ(stacked.Size(), nR * stacked.CoefficientSize());

  for (auto i : f.RadiusIndices()) {
    auto one = SpinField<2, Grid>(grid);
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        one[iTheta, iPhi] = f.Slice(i)[iTheta, iPhi];
      }
    }
    auto alone = Expand(one, lMax);
    for (auto l : alone.Degrees()) {
      for (auto m : alone.Orders(l)) {
        EXPECT_EQ((stacked[i, l, m]), (alone[l, m]))
            << "radius " << i << ", l = " << l << ", m = " << m;
      }
    }
  }
}

TEST(LayeredSpinField, RoundTripsThroughTheSpectralDomain) {
  constexpr auto lMax = Int{6};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = Radii(5);

  auto e = LayeredSpinExpansion<1, Grid>(radial, grid, lMax);
  for (auto i : e.RadiusIndices()) {
    for (auto l : e.Degrees()) {
      for (auto m : e.Orders(l)) {
        e[i, l, m] = Complex{std::cos(0.3 * l + m + i), std::sin(0.2 * l - i)};
      }
    }
  }

  auto field = Evaluate(e);
  auto back = Expand(field, lMax);

  for (auto i : e.RadiusIndices()) {
    for (auto l : e.Degrees()) {
      for (auto m : e.Orders(l)) {
        EXPECT_NEAR((back[i, l, m]).real(), (e[i, l, m]).real(), 1.0e-11);
        EXPECT_NEAR((back[i, l, m]).imag(), (e[i, l, m]).imag(), 1.0e-11);
      }
    }
  }
}

TEST(LayeredSpinField, ThreadingIsTheCallersChoiceAndChangesNothing) {
  constexpr auto lMax = Int{8};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = Radii(9);

  auto e = LayeredSpinExpansion<0, Grid>(radial, grid, lMax);
  for (auto i : e.RadiusIndices()) {
    for (auto l : e.Degrees()) {
      for (auto m : e.Orders(l)) {
        e[i, l, m] = Complex{std::sin(0.4 * l + m), std::cos(0.1 * i)};
      }
    }
  }

  auto sequential = Evaluate(e);
  auto parallel = Evaluate(e, Execution::Parallel(4));

  // The inverse transform's colatitudes write disjoint rows, so threading
  // reorders nothing and the answer is bit-identical.
  for (auto j = Int{0}; j < sequential.Size(); j++) {
    EXPECT_EQ(sequential.Data()[j], parallel.Data()[j]) << "at " << j;
  }
}
