// The conventions, pinned to things a user can check.
//
// docs/gshtrans-reference.tex §2.4 used to list two statements the library
// could not settle from the inside: the sign in
//
//     e_{+-} = -+ (1/sqrt2)(theta-hat +- i phi-hat),   e_0 = r-hat        (*)
//
// and the overall signs of eth and eth-bar. Neither was observed by anything:
// every test worked in canonical components throughout, and a convention that
// nothing reads back into the physical frame is not observable.
//
// This file reads it back. The claim under test is the one a user actually
// depends on -- that the gradient of a function is the gradient of that
// function -- and it is stated where the answer is known independently of
// this library, in the orthonormal (theta-hat, phi-hat) frame. From (*),
//
//     v = v^+ e_+ + v^- e_- + v^0 e_0
//        =>  v_theta = (v^- - v^+)/sqrt2
//            v_phi   = -i (v^+ + v^-)/sqrt2
//            v_r     = v^0
//
// and inverting, v^{+-} = (-+ v_theta + i v_phi)/sqrt2.
//
// What that settles, and what it does not. Flipping (*) and the sign of the
// eth pair *together* is still a relabelling that nothing can see -- it
// renames which component is called +1. What is no longer free is the
// relation between them: given (*) as written, the operator signs are
// determined, and given the operator signs, (*) is. The two were previously
// independent unknowns and are now one convention, fixed here, with the
// physical gradient as the thing that fixes it. Every test below fails by an
// O(1) amount, not by a tolerance, if either sign is flipped alone.
//
// One trap, met while writing these and worth recording. A field written as
// a(theta, phi) theta-hat is almost never smooth on the sphere, because
// theta-hat is not: at the pole its direction depends on the azimuth from
// which the pole is approached. a = sin(theta) cos(phi) looks harmless and is
// not differentiable there, and its spin-weighted expansion does not
// converge -- the divergence check below came out wrong by 0.18 at every
// truncation from lMax = 8 to 128, which is what a non-convergent expansion
// looks like rather than what a bug looks like. The tangential fields here
// are therefore built as grad f + r-hat x grad g from smooth scalars, which
// is smooth by construction.

#include <gtest/gtest.h>

#include <GSHTrans/All>

#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>

namespace {

using namespace GSHTrans;
using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

constexpr auto RootTwo = std::numbers::sqrt2_v<Real>;

// The physical components of a vector, from its canonical contravariant ones.
struct Physical {
  Complex theta, phi, r;
};

template <typename Vector>
Physical ToPhysicalFrame(const Vector& v, Int iTheta, Int iPhi) {
  const auto plus = v.template Component<1>()[iTheta, iPhi];
  const auto minus = v.template Component<-1>()[iTheta, iPhi];
  const auto zero = v.template Component<0>()[iTheta, iPhi];
  return {(minus - plus) / RootTwo, Complex{0, -1} * (plus + minus) / RootTwo,
          zero};
}

// And back again, which is how a tangential field is built below.
struct Canonical {
  Complex plus, minus;
};

Canonical FromPhysicalFrame(Real vTheta, Real vPhi) {
  return {(Complex{-vTheta, 0} + Complex{0, vPhi}) / RootTwo,
          (Complex{vTheta, 0} + Complex{0, vPhi}) / RootTwo};
}

}  // namespace

//--------------------------------------------------------------------------//
//              Rank 0: the gradient of a scalar is the gradient             //
//--------------------------------------------------------------------------//

// grad f = theta-hat d_theta f + phi-hat (1/sin theta) d_phi f, and both
// components are checked pointwise. The second field carries the azimuthal
// dependence the first does not, so a sign error in the phi component cannot
// hide behind a zero.
//
// Every scalar here is a real spherical harmonic at low degree -- P_1, the
// real part of Y_1^1, the real part of Y_2^2 -- which is deliberate and is
// the general rule for this file: a harmonic is smooth at the poles by
// construction, so its expansion converges and the check measures the
// operator rather than the truncation. Trigonometry that merely looks smooth
// is where the hours go.
TEST(Conventions, TheSurfaceGradientOfAScalarIsTheGradient) {
  constexpr auto lMax = Int{12};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  struct Case {
    const char* name;
    Real (*f)(Real, Real);
    Real (*dTheta)(Real, Real);
    Real (*dPhiOverSin)(Real, Real);
  };

  const Case cases[] = {
      {"cos(theta)", [](Real t, Real) { return std::cos(t); },
       [](Real t, Real) { return -std::sin(t); },
       [](Real, Real) { return 0.0; }},
      {"sin(theta) cos(phi)",
       [](Real t, Real p) { return std::sin(t) * std::cos(p); },
       [](Real t, Real p) { return std::cos(t) * std::cos(p); },
       [](Real, Real p) { return -std::sin(p); }},
      // Degree two, and with azimuthal order two, because degree one is
      // special in enough ways that a check resting on it alone is thin.
      {"sin^2(theta) cos(2 phi)",
       [](Real t, Real p) { return std::sin(t) * std::sin(t) * std::cos(2 * p); },
       [](Real t, Real p) { return 2 * std::sin(t) * std::cos(t) * std::cos(2 * p); },
       [](Real t, Real p) { return -2 * std::sin(t) * std::sin(2 * p); }},
  };

  for (const auto& c : cases) {
    auto scalar = TensorField<0, NoSymmetry<0>, ComplexTensor, Grid>(grid);
    for (auto iTheta : grid.CoLatitudeIndices()) {
      const auto theta = grid.CoLatitudes()[iTheta];
      for (auto iPhi : grid.LongitudeIndices()) {
        const auto phi = grid.Longitudes()[iPhi];
        scalar.Component<>()[iTheta, iPhi] = Complex{c.f(theta, phi), 0};
      }
    }

    auto gradient = Evaluate(SurfaceGradient(Expand(scalar, lMax)));

    for (auto iTheta : grid.CoLatitudeIndices()) {
      const auto theta = grid.CoLatitudes()[iTheta];
      for (auto iPhi : grid.LongitudeIndices()) {
        const auto phi = grid.Longitudes()[iPhi];
        const auto v = ToPhysicalFrame(gradient, iTheta, iPhi);

        EXPECT_NEAR(v.theta.real(), c.dTheta(theta, phi), 1.0e-12)
            << c.name << " at theta = " << theta << ", phi = " << phi;
        EXPECT_NEAR(v.phi.real(), c.dPhiOverSin(theta, phi), 1.0e-12)
            << c.name << " at theta = " << theta << ", phi = " << phi;

        // A gradient is real where the scalar is, and the surface gradient
        // has no radial part at all.
        EXPECT_NEAR(v.theta.imag(), 0.0, 1.0e-12);
        EXPECT_NEAR(v.phi.imag(), 0.0, 1.0e-12);
        EXPECT_NEAR(std::abs(v.r), 0.0, 1.0e-12);
      }
    }
  }
}

//--------------------------------------------------------------------------//
//        Rank 1: the divergence of a vector is the divergence               //
//--------------------------------------------------------------------------//

// The metric trace of the surface gradient of a tangential field, against
//
//     div v = (1/sin t) d_t (sin t v_theta) + (1/sin t) d_phi v_phi.
//
// This is the rank-1 statement, and it is where the connection terms of the
// contravariant derivative would show up if they were wrong -- although not
// in the trace itself, where for a tangential field they contribute 2 v^0 and
// so cancel. What it does test, and the rank-0 case cannot, is that a vector
// assembled from physical components transforms, differentiates and contracts
// back to the right physical scalar.
//
// The field is v = grad f + r-hat x grad g, so it is smooth at the poles and
// its divergence is the Laplacian of f alone: the second term is
// divergence-free.
TEST(Conventions, TheTraceOfTheSurfaceGradientIsTheSurfaceDivergence) {
  constexpr auto lMax = Int{16};
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);

  struct Case {
    const char* name;
    Real (*vTheta)(Real, Real);
    Real (*vPhi)(Real, Real);
    Real (*divergence)(Real, Real);
  };

  const Case cases[] = {
      // f = Re Y_1^1 ~ sin(t) cos(p), g = P_1 ~ cos(t). div v = -2 f.
      //
      //   v_theta = d_t f - (1/s) d_p g,   v_phi = (1/s) d_p f + d_t g
      {"degree one",
       [](Real t, Real p) { return std::cos(t) * std::cos(p); },
       [](Real t, Real p) { return -std::sin(p) - std::sin(t); },
       [](Real t, Real p) { return -2 * std::sin(t) * std::cos(p); }},
      // f = P_2 ~ (3 cos^2 t - 1)/2, g = Re Y_2^2 ~ sin^2(t) cos(2p).
      // div v = -6 f.
      {"degree two",
       [](Real t, Real p) {
         return -3 * std::cos(t) * std::sin(t) + 2 * std::sin(t) * std::sin(2 * p);
       },
       [](Real t, Real p) { return 2 * std::sin(t) * std::cos(t) * std::cos(2 * p); },
       [](Real t, Real) {
         return -3 * (3 * std::cos(t) * std::cos(t) - 1);
       }},
  };

  for (const auto& c : cases) {
    auto v = TensorField<1, NoSymmetry<1>, ComplexTensor, Grid>(grid);
    for (auto iTheta : grid.CoLatitudeIndices()) {
      const auto t = grid.CoLatitudes()[iTheta];
      for (auto iPhi : grid.LongitudeIndices()) {
        const auto p = grid.Longitudes()[iPhi];
        const auto canonical = FromPhysicalFrame(c.vTheta(t, p), c.vPhi(t, p));
        v.Component<1>()[iTheta, iPhi] = canonical.plus;
        v.Component<-1>()[iTheta, iPhi] = canonical.minus;
        v.Component<0>()[iTheta, iPhi] = Complex{0, 0};
      }
    }

    // Round-trip the construction first, so that a failure below is about the
    // gradient rather than about the frame conversion.
    for (auto iTheta : grid.CoLatitudeIndices()) {
      const auto t = grid.CoLatitudes()[iTheta];
      for (auto iPhi : grid.LongitudeIndices()) {
        const auto p = grid.Longitudes()[iPhi];
        const auto back = ToPhysicalFrame(v, iTheta, iPhi);
        EXPECT_NEAR(back.theta.real(), c.vTheta(t, p), 1.0e-14) << c.name;
        EXPECT_NEAR(back.phi.real(), c.vPhi(t, p), 1.0e-14) << c.name;
      }
    }

    auto gradient = Evaluate(SurfaceGradient(Expand(v, lMax)));

    // g_{ab} = (-1)^a delta_{a+b,0}, so the trace is -T^{-+} + T^{00} - T^{+-}.
    for (auto iTheta : grid.CoLatitudeIndices()) {
      const auto t = grid.CoLatitudes()[iTheta];
      for (auto iPhi : grid.LongitudeIndices()) {
        const auto p = grid.Longitudes()[iPhi];
        const auto trace = -gradient.Component<-1, 1>()[iTheta, iPhi] +
                           gradient.Component<0, 0>()[iTheta, iPhi] -
                           gradient.Component<1, -1>()[iTheta, iPhi];
        EXPECT_NEAR(trace.real(), c.divergence(t, p), 1.0e-11)
            << c.name << " at theta = " << t << ", phi = " << p;
        EXPECT_NEAR(trace.imag(), 0.0, 1.0e-11) << c.name;
      }
    }
  }
}
