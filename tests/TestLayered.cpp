#include <gtest/gtest.h>

#include <GSHTrans/All>

#include <cmath>
#include <complex>
#include <cstddef>
#include <numeric>
#include <ranges>
#include <span>
#include <stdexcept>
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

// A radial grid with the trapezium rule's weights, so that the quadrature
// tests have something exact to check against.
auto WeightedRadii(Int nR, Real inner = 0.5, Real outer = 1.0) {
  auto r = std::vector<Real>(nR);
  auto w = std::vector<Real>(nR);
  const auto h = (outer - inner) / static_cast<Real>(nR - 1);
  for (Int i = 0; i < nR; i++) {
    r[i] = inner + h * i;
    w[i] = (i == 0 || i == nR - 1) ? h / 2 : h;
  }
  return RadialGrid<Real>(r, w);
}

// A three-point centred difference on a uniform mesh, one-sided at the ends.
// Stands in here for whatever the application actually supplies.
template <typename Scalar>
struct CentredDifference {
  Real h;

  void operator()(std::span<const Scalar> in, std::span<Scalar> out) const {
    const auto n = static_cast<Int>(in.size());
    for (Int i = 0; i < n; i++) {
      if (i == 0) {
        out[i] = (in[1] - in[0]) / h;
      } else if (i == n - 1) {
        out[i] = (in[n - 1] - in[n - 2]) / h;
      } else {
        out[i] = (in[i + 1] - in[i - 1]) / (2 * h);
      }
    }
  }
};

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

//--------------------------------------------------------------------------//
//                              The radial seam                              //
//--------------------------------------------------------------------------//

TEST(RadialOperator, AppliesAlongTheRadialAxisAndNowhereElse) {
  constexpr auto lMax = Int{6};
  constexpr auto nR = Int{33};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  const auto h = 0.5 / (nR - 1);
  auto radial = Radii(nR);

  // f(r, theta, phi) = r^2 * g(theta, phi), so d/dr is 2r * g and the answer
  // is known at every angular point independently.
  auto stack = LayeredSpinField<0, Grid>(radial, grid);
  auto reference = std::vector<Complex>(static_cast<std::size_t>(stack.Size()));
  for (auto i : stack.RadiusIndices()) {
    const auto r = radial.Radius(i);
    auto slice = stack.Slice(i);
    for (Int j = 0; j < slice.Size(); j++) {
      const auto g = Complex{std::sin(0.3 * j), 0.25 * j};
      slice.Data()[j] = r * r * g;
      reference[static_cast<std::size_t>(i * stack.FieldSize() + j)] = 2 * r * g;
    }
  }

  auto derivative = ApplyRadially(stack, CentredDifference<Complex>{h});

  // A centred difference is exact on a quadratic in the interior; the
  // one-sided ends are not, so they are checked to their own order.
  for (auto i : stack.RadiusIndices()) {
    const auto interior = i > 0 && i < nR - 1;
    for (Int j = 0; j < stack.FieldSize(); j++) {
      const auto k = static_cast<std::size_t>(i * stack.FieldSize() + j);
      EXPECT_NEAR(derivative.Data()[k].real(), reference[k].real(),
                  interior ? 1.0e-12 : 2.0e-2)
          << "at radius " << i << ", point " << j;
      EXPECT_NEAR(derivative.Data()[k].imag(), reference[k].imag(),
                  interior ? 1.0e-11 : 1.0);
    }
  }
}

TEST(RadialOperator, ActsTheSameOnEitherSideOfTheTransform) {
  // A radial operator touches only r, and the angular transform touches only
  // (theta, phi), so the two commute. That is the property that lets the model
  // workflow do its radial work in the spectral domain, and it is worth
  // pinning down rather than assuming.
  constexpr auto lMax = Int{8};
  constexpr auto nR = Int{17};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  const auto h = 0.5 / (nR - 1);
  auto radial = Radii(nR);

  auto e = LayeredSpinExpansion<0, Grid>(radial, grid, lMax);
  for (auto i : e.RadiusIndices()) {
    const auto r = radial.Radius(i);
    for (auto l : e.Degrees()) {
      for (auto m : e.Orders(l)) {
        e[i, l, m] = r * r * r * Complex{std::sin(0.4 * l + m), 0.2 * m - 0.1};
      }
    }
  }

  const auto op = CentredDifference<Complex>{h};

  auto spectralFirst = Evaluate(ApplyRadially(e, op));
  auto spatialFirst = ApplyRadially(Evaluate(e), op);

  for (Int j = 0; j < spectralFirst.Size(); j++) {
    EXPECT_NEAR(spectralFirst.Data()[j].real(), spatialFirst.Data()[j].real(),
                1.0e-12);
    EXPECT_NEAR(spectralFirst.Data()[j].imag(), spatialFirst.Data()[j].imag(),
                1.0e-12);
  }
}

TEST(RadialOperator, ThreadingChangesNothing) {
  constexpr auto lMax = Int{8};
  constexpr auto nR = Int{21};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = Radii(nR);

  auto stack = LayeredSpinField<0, Grid>(radial, grid);
  for (Int j = 0; j < stack.Size(); j++) {
    stack.Data()[j] = Complex{std::cos(0.017 * j), std::sin(0.023 * j)};
  }

  const auto op = CentredDifference<Complex>{0.5 / (nR - 1)};
  auto sequential = ApplyRadially(stack, op);
  auto parallel = ApplyRadially(stack, op, Execution::Parallel(4));

  // Each line is independent, so this is a partition of the work and not a
  // reordering of any sum: bit-identical, not merely close.
  for (Int j = 0; j < sequential.Size(); j++) {
    EXPECT_EQ(sequential.Data()[j], parallel.Data()[j]) << "at " << j;
  }
}

TEST(RadialOperator, RefusesWhatItCannotDo) {
  constexpr auto lMax = Int{4};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  const auto op = CentredDifference<Complex>{0.1};

  auto stack = LayeredSpinField<0, Grid>(Radii(6), grid);

  // In place would silently read a half-written line.
  EXPECT_THROW(ApplyRadially(stack, stack, op), std::invalid_argument);

  // A different radial grid is a different radial axis, whatever its length.
  auto other = LayeredSpinField<0, Grid>(Radii(6), grid);
  EXPECT_THROW(ApplyRadially(stack, other, op), std::invalid_argument);
}

TEST(RadialOperator, IntegratesWithTheGridsOwnWeights) {
  constexpr auto lMax = Int{4};
  constexpr auto nR = Int{201};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = WeightedRadii(nR);

  // f(r, x) = r * g(x); the exact integral over [0.5, 1] is 0.375 * g(x).
  auto stack = LayeredSpinField<0, Grid>(radial, grid);
  for (auto i : stack.RadiusIndices()) {
    auto slice = stack.Slice(i);
    for (Int j = 0; j < slice.Size(); j++) {
      slice.Data()[j] = radial.Radius(i) * Complex{1.0 + j, 0.5};
    }
  }

  auto integral = IntegrateRadially(stack);
  ASSERT_EQ(static_cast<Int>(integral.size()), stack.FieldSize());
  for (Int j = 0; j < stack.FieldSize(); j++) {
    EXPECT_NEAR(integral[static_cast<std::size_t>(j)].real(),
                0.375 * (1.0 + j), 1.0e-10);
    EXPECT_NEAR(integral[static_cast<std::size_t>(j)].imag(), 0.375 * 0.5,
                1.0e-12);
  }

  // Weights are optional, and a grid without them says so rather than
  // inventing a rule.
  auto unweighted = LayeredSpinField<0, Grid>(Radii(nR), grid);
  EXPECT_THROW(IntegrateRadially(unweighted), std::invalid_argument);
}

//--------------------------------------------------------------------------//
//                            Layered tensors                                //
//--------------------------------------------------------------------------//

namespace {

// d/dr acting exactly on c * r^power, so that a gradient test measures the
// angular algebra and the seam and not a difference formula's truncation.
struct PowerDerivative {
  Real power;
  std::vector<Real> r;

  void operator()(std::span<const Complex> in, std::span<Complex> out) const {
    for (std::size_t i = 0; i < in.size(); i++) {
      out[i] = power * in[i] / r[i];
    }
  }
};

const auto TestRadii = std::vector<Real>{0.5, 0.75, 1.0};

}  // namespace

TEST(LayeredTensorField, HoldsOneRadiusMajorStackPerStoredComponent) {
  constexpr auto lMax = Int{6};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = RadialGrid<Real>(TestRadii);

  using V = LayeredVectorField<Grid>;
  auto v = V(radial, grid);

  // A real vector stores two components: the pinned-real one at upper index
  // zero and one complex one, which is three real degrees of freedom.
  EXPECT_EQ(V::StoredComponents, 2);

  auto& stack = v.ComponentStack<0>();
  EXPECT_EQ(stack.NumberOfRadii(), 3);
  EXPECT_EQ(stack.FieldSize(), grid.FieldSize());

  // A stored component's stack is an ordinary LayeredSpinField, so writing
  // through a slice writes the tensor.
  auto slice = v.Component<0>(1);
  slice.Data()[0] = 2.5;
  EXPECT_EQ(stack.Slice(1).Data()[0], 2.5);
}

TEST(LayeredTensorField, TransformsComponentByComponentAsOneBatch) {
  constexpr auto lMax = Int{8};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = RadialGrid<Real>(TestRadii);

  auto e = LayeredVectorExpansion<Grid, ComplexTensor>(radial, grid, lMax);
  for (auto i : e.RadiusIndices()) {
    for (auto l : e.ComponentStack<1>().Degrees()) {
      for (auto m : e.ComponentStack<1>().Orders(l)) {
        e.ComponentStack<1>()[i, l, m] = Complex{0.1 * l + m, 0.2 - 0.05 * m};
      }
    }
    for (auto l : e.ComponentStack<0>().Degrees()) {
      for (auto m : e.ComponentStack<0>().Orders(l)) {
        e.ComponentStack<0>()[i, l, m] = Complex{std::cos(0.3 * l), 0.1 * m};
      }
    }
    for (auto l : e.ComponentStack<-1>().Degrees()) {
      for (auto m : e.ComponentStack<-1>().Orders(l)) {
        e.ComponentStack<-1>()[i, l, m] = Complex{0.4, std::sin(0.2 * l + m)};
      }
    }
  }

  auto field = Evaluate(e);
  auto back = Expand(field, lMax);

  for (auto i : e.RadiusIndices()) {
    for (auto l : e.ComponentStack<1>().Degrees()) {
      for (auto m : e.ComponentStack<1>().Orders(l)) {
        const auto was = e.ComponentStack<1>()[i, l, m];
        const auto is = back.ComponentStack<1>()[i, l, m];
        EXPECT_NEAR(is.real(), was.real(), 1.0e-11);
        EXPECT_NEAR(is.imag(), was.imag(), 1.0e-11);
      }
    }
  }
}

//--------------------------------------------------------------------------//
//                              The gradient                                 //
//--------------------------------------------------------------------------//

TEST(LayeredGradient, HasTheKnownCoefficientsOnAScalar) {
  constexpr auto lMax = Int{8};
  constexpr auto l = Int{4};
  constexpr auto m = Int{2};
  const auto a = Real{3};

  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = RadialGrid<Real>(TestRadii);

  // f_{lm}(r) = r^a in one (l, m) and zero elsewhere.
  auto f = LayeredScalarExpansion<Grid, ComplexTensor>(radial, grid, lMax);
  for (auto i : radial.RadiusIndices()) {
    f.ComponentStack<>()[i, l, m] = std::pow(TestRadii[i], a);
  }

  auto g = Gradient(f, PowerDerivative{a, TestRadii});

  // (grad f)^0 = df/dr, and (grad f)^{+-} = r^{-1} Omega^{0}_l f with
  // Omega^{0}_l = sqrt(l(l+1)/2), the same for both signs because a scalar
  // has upper index zero.
  const auto omega = std::sqrt(l * (l + 1.0) / 2);
  for (auto i : radial.RadiusIndices()) {
    const auto r = TestRadii[i];
    EXPECT_NEAR((g.Coefficient<0>(i, l, m)).real(), a * std::pow(r, a - 1),
                1.0e-12);
    EXPECT_NEAR((g.Coefficient<1>(i, l, m)).real(), omega * std::pow(r, a - 1),
                1.0e-12);
    EXPECT_NEAR((g.Coefficient<-1>(i, l, m)).real(), omega * std::pow(r, a - 1),
                1.0e-12);
  }
}

TEST(LayeredGradient, TwoGradientsContractToTheLaplacian) {
  // The decisive test of the whole construction. Contracting grad grad f with
  // the metric must give
  //
  //     lap (r^a Y_{lm}) = [a(a+1) - l(l+1)] r^{a-2},
  //
  // and it is decisive because the *second* gradient acts on a vector whose
  // e_0 component is not zero. Getting the right answer needs the r^{-1}, the
  // radial block, and the connection terms that move a slot between e_0 and
  // e_{+-} -- the last of which nothing before this exercised at all.
  constexpr auto lMax = Int{10};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = RadialGrid<Real>(TestRadii);

  struct Case {
    Real a;
    Int l;
    Int m;
  };
  const auto cases = std::vector<Case>{{3, 4, 2},  {1, 1, 0},  {-2, 6, -5},
                                       {0, 2, 1},  {2.5, 0, 0}, {-1, 10, 10}};

  for (const auto& c : cases) {
    auto f = LayeredScalarExpansion<Grid, ComplexTensor>(radial, grid, lMax);
    for (auto i : radial.RadiusIndices()) {
      f.ComponentStack<>()[i, c.l, c.m] = std::pow(TestRadii[i], c.a);
    }

    auto g = Gradient(f, PowerDerivative{c.a, TestRadii});
    auto h = Gradient(g, PowerDerivative{c.a - 1, TestRadii});

    // g_{alpha beta} = (-1)^alpha delta_{alpha + beta, 0}.
    for (auto i : radial.RadiusIndices()) {
      const auto trace = h.Coefficient<0, 0>(i, c.l, c.m) -
                         h.Coefficient<1, -1>(i, c.l, c.m) -
                         h.Coefficient<-1, 1>(i, c.l, c.m);
      const auto expected = (c.a * (c.a + 1) - c.l * (c.l + 1.0)) *
                            std::pow(TestRadii[i], c.a - 2);
      EXPECT_NEAR(trace.real(), expected, 1.0e-10 * (1 + std::abs(expected)))
          << "a = " << c.a << ", l = " << c.l << ", m = " << c.m
          << ", r = " << TestRadii[i];
      EXPECT_NEAR(trace.imag(), 0.0, 1.0e-10);
    }
  }
}

TEST(LayeredGradient, HoldsOnARealScalarThroughTheReducedStorage) {
  // The same identity with the reality condition switched on, so that the
  // scalar's block holds only m >= 0 and the gradient's negative-index
  // components are derived rather than stored.
  constexpr auto lMax = Int{8};
  constexpr auto l = Int{5};
  const auto a = Real{2};

  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = RadialGrid<Real>(TestRadii);

  auto f = LayeredScalarExpansion<Grid>(radial, grid, lMax);
  EXPECT_EQ(decltype(f)::StoredComponents, 1);
  for (auto i : radial.RadiusIndices()) {
    f.ComponentStack<>()[i, l, 3] = std::pow(TestRadii[i], a);
  }

  auto g = Gradient(f, PowerDerivative{a, TestRadii});
  auto h = Gradient(g, PowerDerivative{a - 1, TestRadii});

  for (auto i : radial.RadiusIndices()) {
    for (auto m : {Int{3}, Int{-3}}) {
      const auto trace = h.Coefficient<0, 0>(i, l, m) -
                         h.Coefficient<1, -1>(i, l, m) -
                         h.Coefficient<-1, 1>(i, l, m);
      const auto f0 = f.Coefficient<>(i, l, m);
      const auto expected =
          (a * (a + 1) - l * (l + 1.0)) * std::pow(TestRadii[i], a - 2) * f0 /
          std::pow(TestRadii[i], a);
      EXPECT_NEAR(trace.real(), expected.real(), 1.0e-10);
      EXPECT_NEAR(trace.imag(), expected.imag(), 1.0e-10);
    }
  }
}

TEST(LayeredGradient, RefusesTheOrigin) {
  constexpr auto lMax = Int{4};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  const auto radii = std::vector<Real>{0.0, 0.5, 1.0};
  auto radial = RadialGrid<Real>(radii);

  auto f = LayeredScalarExpansion<Grid>(radial, grid, lMax);

  // r^{-1} is in the operator, not in the field: the basis is singular at the
  // origin and saying so beats returning an infinity.
  EXPECT_THROW(Gradient(f, PowerDerivative{1, radii}), std::invalid_argument);
}

//--------------------------------------------------------------------------//
//                            The other layout                               //
//--------------------------------------------------------------------------//

TEST(RadialMajor, IsTheSameDataWithTheAxesExchanged) {
  constexpr auto lMax = Int{10};
  constexpr auto nR = Int{37};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = Radii(nR);

  auto e = LayeredSpinExpansion<0, Grid>(radial, grid, lMax);
  for (auto i : e.RadiusIndices()) {
    for (auto l : e.Degrees()) {
      for (auto m : e.Orders(l)) {
        e[i, l, m] = Complex{0.3 * i + l, 0.1 * m - 0.7 * i};
      }
    }
  }

  auto major = RadialMajor(e);
  EXPECT_EQ(major.NumberOfRadii(), nR);
  EXPECT_EQ(major.NumberOfLines(), e.CoefficientSize());

  // The line at a degree and order is that coefficient over every radius, and
  // it is contiguous, which is the whole reason the layout exists.
  for (auto l : e.Degrees()) {
    for (auto m : e.Orders(l)) {
      const auto line = major.Line(e.CoefficientIndex(l, m));
      ASSERT_EQ(static_cast<Int>(line.size()), nR);
      for (auto i : e.RadiusIndices()) {
        EXPECT_EQ(line[static_cast<std::size_t>(i)], (e[i, l, m]))
            << "at l = " << l << ", m = " << m << ", radius " << i;
      }
    }
  }

  // Round trip, exactly: a repack moves values and does not compute with them.
  auto back = e.SameShape();
  major.CopyInto(back);
  for (Int j = 0; j < e.Size(); j++) {
    EXPECT_EQ(back.Data()[j], e.Data()[j]) << "at " << j;
  }
}

TEST(RadialMajor, TransposesTheSameWhetherThreadedOrNot) {
  constexpr auto lMax = Int{12};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto e = LayeredSpinExpansion<0, Grid>(Radii(29), grid, lMax);
  for (Int j = 0; j < e.Size(); j++) {
    e.Data()[j] = Complex{std::cos(0.013 * j), std::sin(0.019 * j)};
  }

  auto sequential = RadialMajor(e);
  auto parallel = RadialMajor(e, Execution::Parallel(4));
  for (Int j = 0; j < sequential.Size(); j++) {
    EXPECT_EQ(sequential.Data()[j], parallel.Data()[j]) << "at " << j;
  }
}

TEST(RadialMajor, GivesTheSameAnswerAsApplyingThroughTheGather) {
  // The two routes to a radial operator: gather each line out of the
  // radius-major stack, or transpose once and work on contiguous lines. They
  // must agree exactly, since neither changes the arithmetic -- only when the
  // copying happens.
  constexpr auto lMax = Int{10};
  constexpr auto nR = Int{33};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = Radii(nR);
  const auto op = CentredDifference<Complex>{0.5 / (nR - 1)};

  auto e = LayeredSpinExpansion<0, Grid>(radial, grid, lMax);
  for (Int j = 0; j < e.Size(); j++) {
    e.Data()[j] = Complex{std::sin(0.007 * j), 0.4 - 0.011 * j};
  }

  auto throughGather = ApplyRadially(e, op);

  auto major = RadialMajor(e);
  auto applied = RadialMajor(e);
  ApplyToLines(major, applied, op, Execution::Parallel(4));
  auto throughTranspose = e.SameShape();
  applied.CopyInto(throughTranspose);

  for (Int j = 0; j < e.Size(); j++) {
    EXPECT_EQ(throughGather.Data()[j], throughTranspose.Data()[j])
        << "at " << j;
  }
}

TEST(RadialMajor, RefusesAShapeItDidNotComeFrom) {
  constexpr auto lMax = Int{6};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto e = LayeredSpinExpansion<0, Grid>(Radii(11), grid, lMax);
  auto major = RadialMajor(e);

  auto shorter = LayeredSpinExpansion<0, Grid>(Radii(9), grid, lMax);
  EXPECT_THROW(major.CopyInto(shorter), std::invalid_argument);
}

TEST(RadialMajor, RefillsWithoutAllocating) {
  constexpr auto lMax = Int{8};
  constexpr auto nR = Int{19};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = Radii(nR);

  auto first = LayeredSpinExpansion<0, Grid>(radial, grid, lMax);
  auto second = LayeredSpinExpansion<0, Grid>(radial, grid, lMax);
  for (Int j = 0; j < first.Size(); j++) {
    first.Data()[j] = Complex{0.5 * j, 1.0};
    second.Data()[j] = Complex{-0.25 * j, 2.0};
  }

  auto major = RadialMajor(first);
  const auto* before = major.Data().data();

  major.CopyFrom(second);
  EXPECT_EQ(major.Data().data(), before) << "refilling must not reallocate";

  auto back = second.SameShape();
  major.CopyInto(back);
  for (Int j = 0; j < second.Size(); j++) {
    EXPECT_EQ(back.Data()[j], second.Data()[j]) << "at " << j;
  }

  // A scratch buffer of the same shape, transposed from nothing.
  auto scratch = major.SameShape();
  EXPECT_EQ(scratch.NumberOfRadii(), major.NumberOfRadii());
  EXPECT_EQ(scratch.NumberOfLines(), major.NumberOfLines());

  auto other = LayeredSpinExpansion<0, Grid>(Radii(nR + 2), grid, lMax);
  EXPECT_THROW(major.CopyFrom(other), std::invalid_argument);
}

//--------------------------------------------------------------------------//
//                        Tangential layered tensors                         //
//--------------------------------------------------------------------------//

namespace {

template <typename Reality>
using TangentialLayered =
    LayeredTensorField<2, NoSymmetry<2>, Reality, Grid, TangentialSlots>;

template <typename Reality>
using TangentialLayeredExpansion =
    LayeredTensorExpansion<2, NoSymmetry<2>, Reality, Grid, TangentialSlots>;

// Whether the layered gradient will take an operand at all, asked as a
// concept so the negative case is an unsatisfied requirement.
template <typename E>
concept SurfaceDifferentiable = requires(const E& e) { SurfaceGradient(e); };

}  // namespace

// The layered type adds a radial axis and nothing else, so the alphabet only
// has to be passed through: which components exist is the flat type's
// question, and the counts here are the flat ones.
TEST(LayeredTensorField, ATangentialStackHasOneStackPerTangentialComponent) {
  constexpr auto lMax = Int{6};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = RadialGrid<Real>(TestRadii);

  using T = TangentialLayered<RealTensor>;
  static_assert(std::same_as<T::SlotSet, TangentialSlots>);
  static_assert(T::Components == 4);
  static_assert(T::StoredComponents == 2);

  auto t = T(radial, grid);
  auto& stack = t.ComponentStack<-1, 1>();
  EXPECT_EQ(stack.NumberOfRadii(), 3);
  EXPECT_EQ(stack.FieldSize(), grid.FieldSize());

  // Writing through a slice writes the tensor, as on the general type.
  auto slice = t.Component<-1, 1>(1);
  slice.Data()[0] = Complex{2.5, -1.5};
  EXPECT_EQ(stack.Slice(1).Data()[0], (Complex{2.5, -1.5}));

  // And a radial index is refused here too, the accessors being conditioned
  // on the flat type's Represents.
  static_assert(!T::Represents<0, 1>);
  static_assert(!T::Writable<0, 1>);
}

// Both bridges, over a tangential alphabet: one batched transform per stored
// component, and the coefficients come back as they went in.
TEST(LayeredTensorField, ATangentialStackRoundTripsThroughBothBridges) {
  constexpr auto lMax = Int{8};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = RadialGrid<Real>(TestRadii);

  auto e = TangentialLayeredExpansion<RealTensor>(radial, grid, lMax);
  const auto fill = [&](auto&& stack, Real tag) {
    for (auto i : e.RadiusIndices()) {
      for (auto l : stack.Degrees()) {
        for (auto m : stack.Orders(l)) {
          stack[i, l, m] = Complex{tag + std::cos(0.3 * l), 0.1 * m - tag};
        }
      }
    }
  };
  fill(e.ComponentStack<-1, -1>(), 1.0);
  fill(e.ComponentStack<-1, 1>(), 2.0);

  auto field = Evaluate(e);
  static_assert(std::same_as<decltype(field)::SlotSet, TangentialSlots>);
  auto back = Expand(field, lMax);
  static_assert(std::same_as<decltype(back)::SlotSet, TangentialSlots>);

  const auto compare = [&](auto&& was, auto&& is) {
    for (auto i : e.RadiusIndices()) {
      for (auto l : was.Degrees()) {
        for (auto m : was.Orders(l)) {
          // Named first: a multi-argument subscript inside a macro reads as
          // three macro arguments.
          const auto before = was[i, l, m];
          const auto after = is[i, l, m];
          EXPECT_NEAR(after.real(), before.real(), 1.0e-11);
          EXPECT_NEAR(after.imag(), before.imag(), 1.0e-11);
        }
      }
    }
  };
  compare(e.ComponentStack<-1, -1>(), back.ComponentStack<-1, -1>());
  compare(e.ComponentStack<-1, 1>(), back.ComponentStack<-1, 1>());
}

// [D9] again, on the layered side: grad_1 moves slots between e_0 and e_+-,
// so it does not close on the tangential bundle and takes no operand from it.
TEST(LayeredGradient, TakesNoTangentialOperand) {
  using General = LayeredTensorExpansion<2, NoSymmetry<2>, ComplexTensor, Grid>;
  using Tangential = TangentialLayeredExpansion<ComplexTensor>;

  static_assert(SurfaceDifferentiable<General>);
  static_assert(!SurfaceDifferentiable<Tangential>);
  SUCCEED();
}

#ifdef GSHTRANS_HAVE_INTERPOLATION
#include <Interpolation/CubicSpline.hpp>
#endif

//--------------------------------------------------------------------------//
//                         Ready-made radial derivatives                     //
//--------------------------------------------------------------------------//

namespace {

// Deliberately unequally spaced, and not close to uniform: a rule that
// happened to assume equal spacing would pass on a uniform grid and fail
// here, which is the point of choosing these.
const auto UnevenRadii =
    std::vector<Real>{0.40, 0.55, 0.70, 1.00, 1.30, 1.45};

Real Monomial(Real r, Int degree) { return std::pow(r, degree); }
Real MonomialSlope(Real r, Int degree) {
  return degree == 0 ? Real{0} : degree * std::pow(r, degree - 1);
}

// Apply an operator to one line, which is what a radial operator is.
template <typename Op>
auto Line(const Op& op, const std::vector<Real>& in) {
  auto out = std::vector<Real>(in.size());
  op(std::span<const Real>(in), std::span<Real>(out));
  return out;
}

}  // namespace

// The property that makes a finite-difference rule what it is: exact for
// polynomials up to its order, at *every* node including the ends, where the
// stencil is one-sided. A rule that was only centred would have nothing to
// say at the first and last radii -- which are exactly the radii a boundary
// condition is applied at.
TEST(RadialDerivatives, FiniteDifferencesAreExactToTheirOrder) {
  auto radial = RadialGrid<Real>(UnevenRadii);

  for (auto order : {Int{1}, Int{2}, Int{3}, Int{4}}) {
    const auto d = FiniteDifferenceDerivative<Real>(radial, order);
    EXPECT_EQ(d.Order(), order);

    for (auto degree = Int{0}; degree <= order; degree++) {
      auto values = std::vector<Real>{};
      for (auto r : UnevenRadii) values.push_back(Monomial(r, degree));
      const auto got = Line(d, values);
      for (std::size_t i = 0; i < UnevenRadii.size(); i++) {
        EXPECT_NEAR(got[i], MonomialSlope(UnevenRadii[i], degree), 1.0e-11)
            << "order " << order << ", degree " << degree << ", node " << i;
      }
    }

    // And one degree higher is not exact, which is what says the order means
    // something rather than the test being satisfied by any weights at all.
    auto tooHigh = std::vector<Real>{};
    for (auto r : UnevenRadii) tooHigh.push_back(Monomial(r, order + 1));
    const auto got = Line(d, tooHigh);
    auto worst = Real{0};
    for (std::size_t i = 0; i < UnevenRadii.size(); i++) {
      worst = std::max(worst, std::abs(got[i] - MonomialSlope(UnevenRadii[i],
                                                              order + 1)));
    }
    EXPECT_GT(worst, 1.0e-9) << "order " << order;
  }
}

TEST(RadialDerivatives, FiniteDifferencesRefuseWhatTheyCannotDo) {
  auto radial = RadialGrid<Real>(UnevenRadii);
  EXPECT_THROW((FiniteDifferenceDerivative<Real>(radial, 0)),
               std::invalid_argument);
  // A stencil wider than the grid.
  EXPECT_THROW((FiniteDifferenceDerivative<Real>(radial, 6)),
               std::invalid_argument);
  EXPECT_NO_THROW((FiniteDifferenceDerivative<Real>(radial, 5)));
}

// The differentiation matrix is exact to the highest degree any operator on
// these nodes could be, which is one less than their number.
TEST(RadialDerivatives, TheDifferentiationMatrixIsExactToTheNodeCount) {
  auto radial = RadialGrid<Real>(UnevenRadii);
  const auto d = LagrangeDerivative<Real>(radial);
  const auto n = static_cast<Int>(UnevenRadii.size());

  for (auto degree = Int{0}; degree < n; degree++) {
    auto values = std::vector<Real>{};
    for (auto r : UnevenRadii) values.push_back(Monomial(r, degree));
    const auto got = Line(d, values);
    for (std::size_t i = 0; i < UnevenRadii.size(); i++) {
      EXPECT_NEAR(got[i], MonomialSlope(UnevenRadii[i], degree), 1.0e-10)
          << "degree " << degree << ", node " << i;
    }
  }

  // A constant differentiates to zero to rounding, and the diagonal is what
  // makes that so: it is imposed as minus the row sum rather than evaluated
  // from the formula, which is the negative-sum trick. What it buys is that
  // the error here is the rounding of *this* sum -- a few epsilon, and the
  // summation order is why it is not identically zero -- rather than the
  // accuracy of a closed form for the diagonal, which is what would otherwise
  // grow with the number of nodes.
  auto ones = std::vector<Real>(UnevenRadii.size(), Real{1});
  const auto zero = Line(d, ones);
  for (auto value : zero) EXPECT_NEAR(value, Real{0}, 1.0e-14);
}

// Radial lines are complex in the spectral domain, which is where the model
// application does its radial work, so an operator that only took real lines
// would be useless where it is most wanted.
TEST(RadialDerivatives, ActOnComplexLinesAsReadilyAsRealOnes) {
  auto radial = RadialGrid<Real>(UnevenRadii);
  const auto fd = FiniteDifferenceDerivative<Real>(radial, 2);
  const auto lagrange = LagrangeDerivative<Real>(radial);

  auto values = std::vector<Complex>{};
  for (auto r : UnevenRadii) values.push_back(Complex{r * r, 3.0 * r});

  for (const auto& apply : {std::function<void(std::span<const Complex>,
                                               std::span<Complex>)>(
                                [&](auto in, auto out) { fd(in, out); }),
                            std::function<void(std::span<const Complex>,
                                               std::span<Complex>)>(
                                [&](auto in, auto out) { lagrange(in, out); })}) {
    auto got = std::vector<Complex>(values.size());
    apply(std::span<const Complex>(values), std::span<Complex>(got));
    for (std::size_t i = 0; i < UnevenRadii.size(); i++) {
      EXPECT_NEAR(got[i].real(), 2.0 * UnevenRadii[i], 1.0e-10) << i;
      EXPECT_NEAR(got[i].imag(), 3.0, 1.0e-10) << i;
    }
  }
}

// The contract of RadialOperator.h: one operator, shared const, called from
// every thread, with whatever scratch it needs in thread_local storage. This
// is the pattern the pre-built operators are written to and that a caller
// writing their own has to follow, so it is pinned rather than described.
TEST(RadialOperator, OneOperatorServesEveryThread) {
  constexpr auto lMax = Int{8};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = RadialGrid<Real>(UnevenRadii);

  // An operator that genuinely needs scratch, so that a mutable member would
  // be visibly wrong rather than accidentally right.
  struct ScratchOperator {
    void operator()(std::span<const Complex> in,
                    std::span<Complex> out) const {
      thread_local auto work = std::vector<Complex>{};
      if (work.size() < in.size()) work.resize(in.size());
      for (std::size_t i = 0; i < in.size(); i++) {
        work[i] = in[i] * Complex{0.0, 1.0};
      }
      for (std::size_t i = 0; i < in.size(); i++) {
        out[i] = work[in.size() - 1 - i];
      }
    }
  };

  auto f = LayeredSpinField<0, Grid, ComplexValued>(radial, grid);
  for (Int i = 0; i < f.NumberOfRadii() * f.SliceSize(); i++) {
    f.Data()[i] = Complex{std::cos(0.03 * i), std::sin(0.07 * i)};
  }

  const auto op = ScratchOperator{};
  const auto sequential = ApplyRadially(f, op, Execution::Sequential());
  const auto threaded = ApplyRadially(f, op, Execution::Parallel());

  ASSERT_EQ(sequential.Data().size(), threaded.Data().size());
  for (std::size_t i = 0; i < sequential.Data().size(); i++) {
    EXPECT_EQ(sequential.Data()[i], threaded.Data()[i]) << "at " << i;
  }

  // And the ready-made ones are usable the same way, which is the point.
  const auto d = FiniteDifferenceDerivative<Real>(radial, 2);
  const auto one = ApplyRadially(f, d, Execution::Sequential());
  const auto many = ApplyRadially(f, d, Execution::Parallel());
  for (std::size_t i = 0; i < one.Data().size(); i++) {
    EXPECT_EQ(one.Data()[i], many.Data()[i]) << "at " << i;
  }
}

// What the whole exercise is for: Gradient now runs without the caller
// writing a differentiation matrix first. The identity is the Laplacian one,
// at a power both operators integrate exactly -- r^2 is degree two, and both
// a three-point rule and a three-node matrix are exact there, so the answer
// is machine precision rather than a truncation error to be tolerated.
TEST(LayeredGradient, RunsWithAReadyMadeRadialDerivative) {
  constexpr auto lMax = Int{8};
  constexpr auto l = Int{3};
  constexpr auto m = Int{-2};
  const auto a = Real{2};

  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = RadialGrid<Real>(TestRadii);

  auto f = LayeredScalarExpansion<Grid, ComplexTensor>(radial, grid, lMax);
  for (auto i : radial.RadiusIndices()) {
    f.ComponentStack<>()[i, l, m] = std::pow(TestRadii[i], a);
  }

  const auto fd = FiniteDifferenceDerivative<Real>(radial, 2);
  const auto lagrange = LagrangeDerivative<Real>(radial);

  const auto check = [&](const auto& op, const char* what) {
    auto g = Gradient(f, op);
    auto h = Gradient(g, op);
    for (auto i : radial.RadiusIndices()) {
      // `template`, because inside a generic lambda `h` is dependent.
      const auto trace = h.template Coefficient<0, 0>(i, l, m) -
                         h.template Coefficient<1, -1>(i, l, m) -
                         h.template Coefficient<-1, 1>(i, l, m);
      const auto expected =
          (a * (a + 1) - l * (l + 1.0)) * std::pow(TestRadii[i], a - 2);
      EXPECT_NEAR(trace.real(), expected, 1.0e-10 * (1 + std::abs(expected)))
          << what << " at r = " << TestRadii[i];
      EXPECT_NEAR(trace.imag(), 0.0, 1.0e-10) << what;
    }
  };

  check(fd, "finite differences");
  check(lagrange, "the differentiation matrix");
}

// A natural spline is exact on straight lines and on nothing else, since the
// end conditions force the second derivative to vanish where a curve's does
// not. So the property test is the linear one, and everything beyond it is
// the oracle below.
TEST(RadialDerivatives, TheSplineIsExactOnStraightLines) {
  auto radial = RadialGrid<Real>(UnevenRadii);
  const auto d = SplineDerivative<Real>(radial);

  for (auto degree = Int{0}; degree <= 1; degree++) {
    auto values = std::vector<Real>{};
    for (auto r : UnevenRadii) values.push_back(Monomial(r, degree));
    const auto got = Line(d, values);
    for (std::size_t i = 0; i < UnevenRadii.size(); i++) {
      EXPECT_NEAR(got[i], MonomialSlope(UnevenRadii[i], degree), 1.0e-12)
          << "degree " << degree << ", node " << i;
    }
  }
}

TEST(RadialDerivatives, TheSplineRefusesARepeatedRadius) {
  // A repeated radius is how a two-sided material interface is written, and a
  // single spline through it is not what is wanted there. Refusing says so.
  auto interface = std::vector<Real>{0.4, 0.7, 0.7, 1.0};
  auto radial = RadialGrid<Real>(interface);
  EXPECT_THROW((SplineDerivative<Real>(radial)), std::invalid_argument);
}

#ifdef GSHTRANS_HAVE_INTERPOLATION
// The oracle, and the reason the operator is allowed to be hand-written at
// all: an independent implementation of the same spline, by another author in
// another library, must agree with it to rounding. Anything this file could
// asssert about the operator on its own terms would be weaker.
TEST(RadialDerivatives, TheSplineAgreesWithAnIndependentImplementation) {
  auto radial = RadialGrid<Real>(UnevenRadii);
  const auto d = SplineDerivative<Real>(radial);

  // Data with no polynomial structure, so that agreement is not an accident
  // of both being exact on it.
  auto values = std::vector<Real>{};
  for (auto r : UnevenRadii) values.push_back(std::exp(r) * std::sin(3.0 * r));

  const auto got = Line(d, values);
  const Interpolation::CubicSpline oracle{UnevenRadii, values};
  for (std::size_t i = 0; i < UnevenRadii.size(); i++) {
    EXPECT_NEAR(got[i], oracle.Evaluate<1>(UnevenRadii[i]), 1.0e-12)
        << "node " << i;
  }
}

// And on complex lines, which is what the spectral domain hands it. The
// oracle takes complex ordinates directly, so the comparison is like for
// like rather than a real-and-imaginary split on either side.
TEST(RadialDerivatives, TheSplineAgreesOnComplexLinesToo) {
  auto radial = RadialGrid<Real>(UnevenRadii);
  const auto d = SplineDerivative<Real>(radial);

  auto values = std::vector<Complex>{};
  for (auto r : UnevenRadii) {
    values.push_back(Complex{std::exp(r) * std::sin(3.0 * r), std::cos(2.0 * r)});
  }

  auto got = std::vector<Complex>(values.size());
  d(std::span<const Complex>(values), std::span<Complex>(got));

  const Interpolation::CubicSpline oracle{UnevenRadii, values};
  for (std::size_t i = 0; i < UnevenRadii.size(); i++) {
    const auto expected = oracle.Evaluate<1>(UnevenRadii[i]);
    EXPECT_NEAR(got[i].real(), expected.real(), 1.0e-12) << "node " << i;
    EXPECT_NEAR(got[i].imag(), expected.imag(), 1.0e-12) << "node " << i;
  }
}
#endif

// The spline serves the gradient like the other two, and one operator serves
// every thread.
TEST(RadialDerivatives, TheSplineIsAnOrdinaryRadialOperator) {
  constexpr auto lMax = Int{6};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto radial = RadialGrid<Real>(UnevenRadii);

  auto f = LayeredSpinField<0, Grid, ComplexValued>(radial, grid);
  for (Int i = 0; i < f.NumberOfRadii() * f.SliceSize(); i++) {
    f.Data()[i] = Complex{std::cos(0.05 * i), std::sin(0.11 * i)};
  }

  const auto d = SplineDerivative<Real>(radial);
  const auto one = ApplyRadially(f, d, Execution::Sequential());
  const auto many = ApplyRadially(f, d, Execution::Parallel());
  for (std::size_t i = 0; i < one.Data().size(); i++) {
    EXPECT_EQ(one.Data()[i], many.Data()[i]) << "at " << i;
  }
}
