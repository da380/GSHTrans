#ifndef GSH_TRANS_ETH_GUARD_H
#define GSH_TRANS_ETH_GUARD_H

#include <cmath>
#include <complex>
#include <cstddef>

#include "../Concepts.h"
#include "../Utility.h"
#include "SpinExpansion.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                        Raising and lowering the index                     //
//--------------------------------------------------------------------------//

// The operators that connect different upper indices; see section 6 of the
// theory note, docs/canonical-components.tex.
//
// They are the reason the library has two representations rather than one.
// Everything in the field algebra is local in (theta, phi); these are local in
// (l, m) and not in position, so "the gradient of a product" is necessarily
// evaluated in both -- the product in one, the gradient in the other. Nothing
// here can be a spin-weighted node.
//
// In the spectral domain each is a multiplication of the coefficients by an
// l-dependent factor with no coupling between different (l, m):
//
//   eth      Y^N_{lm} = -sqrt((l - N)(l + N + 1)) Y^{N+1}_{lm}
//   eth-bar  Y^N_{lm} = +sqrt((l + N)(l - N + 1)) Y^{N-1}_{lm}
//
// **The overall signs used to be an open convention, and are not any more.**
// Two identities are checked and pin everything but the common sign: eth-bar
// eth is the surface Laplacian on a scalar, the factors multiplying to
// -l(l+1), and the commutator is -2N. Both survive flipping the pair, which
// is why the sign was long carried as a question for Phinney & Burridge.
//
// It is not their question. P&B do not use eth as such, and what has to be
// right is that gradients of functions come out as gradients -- which is
// checkable here, and is checked, by reading a gradient back into the
// physical (theta-hat, phi-hat) frame through the stated e_{+-} convention
// and requiring it to equal d_theta f and (sin theta)^{-1} d_phi f pointwise.
// That fixes the sign here against the sign in e_{+-}: flipping either alone
// fails by an O(1) amount, flipping both together is a relabelling of which
// component is called +1. tests/TestConventions.cpp is where it is observed,
// and section 2.4 of docs/gshtrans-reference.tex states it.
//
// The degree ranges look after themselves. Raising a field at N >= 0 gives a
// field starting at l = N + 1, and the factor at l = N is sqrt(0) = 0, so the
// coefficient that has nowhere to go was going to be zero anyway. Lowering
// gives a field starting at a *lower* degree, whose extra coefficient has no
// source; it is zero, which is what a field with no content there means.
//
// **These are not the tensor gradient.** Phinney & Burridge, and Dahlen &
// Tromp after them, work with a contravariant derivative d^sigma taking the
// coefficients of a rank-q tensor to those of the rank-(q+1) tensor grad T.
// For a *scalar* it agrees with these operators up to a factor of sqrt(2) --
// the normalisation of e_{+-} -- but for rank one and above it is not a
// multiplication at all: it also subtracts the components with one slot index
// shifted, once per slot, because the canonical basis vectors themselves vary
// over the sphere and differentiating a tensor field differentiates its basis.
//
// So applying Raise or Lower to each component of a tensor expansion does not
// give grad T. Nothing in the types says so, because a tensor's components are
// ordinary spin expansions. Use SurfaceGradient, which is d^sigma; section 7
// of docs/gshtrans-reference.tex gives it in full.
//
// **But there is a case where component-wise eth is exactly right**, and it
// would be a mistake to read the paragraph above as forbidding it. There are
// two derivatives on the sphere. d^sigma is the *ambient* one: it
// differentiates the tensor and its basis, and the basis leaves the tangent
// plane. The *intrinsic* covariant derivative -- the Levi-Civita connection of
// the induced metric, which is what surface differential geometry means by
// differentiation -- is closed on tangential tensors, those with no radial
// slot. On such a tensor every slot shift in d^sigma either leaves the
// alphabet or lands on a radial slot the tensor does not have, so the
// connection terms all vanish and what is left is pure Omega multiplication:
//
//   on a tangential tensor, the intrinsic covariant derivative is Raise and
//   Lower applied component by component, up to sqrt(2), and it is closed.
//
// The difference between the two is algebraic, not differential: for a
// tangential T the radial-slot components of d^sigma T are -T with that slot
// replaced, which is the extrinsic curvature of the sphere. So eth is not a
// poor relation of the contravariant derivative; it is the intrinsic
// derivative, exact on exactly the objects -- tangential, spin-weighted --
// that it was invented for.

namespace EthDetails {

using Int = std::ptrdiff_t;

template <RealFloatingPoint Real>
constexpr Real RaisingFactor(Int l, Int n) {
  const auto a = static_cast<Real>(l - n);
  const auto b = static_cast<Real>(l + n + 1);
  return -std::sqrt(a * b);
}

template <RealFloatingPoint Real>
constexpr Real LoweringFactor(Int l, Int n) {
  const auto a = static_cast<Real>(l + n);
  const auto b = static_cast<Real>(l - n + 1);
  return std::sqrt(a * b);
}

// Omega^{s}_l = sqrt((l + s)(l - s + 1) / 2), so that Omega^{+N} is
// Omega(l, N) and Omega^{-N} is Omega(l, -N).
//
// D&T's coefficient rather than eth's, and it lives here beside eth's own two
// because both derivatives that use it -- the contravariant one and the
// intrinsic one -- are built from the same raising and lowering, and this is
// where a reader looking for the factor will come.
template <RealFloatingPoint Real>
Real Omega(Int l, Int s) {
  const auto a = static_cast<Real>(l + s);
  const auto b = static_cast<Real>(l - s + 1);
  return std::sqrt(a * b / 2);
}

}  // namespace EthDetails

// A coefficient of an expansion, or zero where it has none.
//
// Lowering reaches degrees below the source's minimum, where the source has
// no coefficient because the harmonic does not exist. That is not a missing
// value to be looked up but a zero: a field of upper index N has no content
// below degree |N|.
//
// It also fills in the negative orders a real field does not store, by
// f_{l,-m} = (-1)^m conj(f_{lm}), so that raising or lowering a real scalar
// gives the right thing rather than reading off the end.
template <std::ptrdiff_t N, AngularGrid Grid, RealOrComplexValued Value>
auto Coefficient(const SpinExpansion<N, Grid, Value>& expansion,
                 std::ptrdiff_t l, std::ptrdiff_t m) {
  using Complex = std::complex<typename Grid::Real>;
  if (l < expansion.MinDegree() || l > expansion.MaxDegree()) return Complex{};
  if (m < -l || m > l) return Complex{};
  if constexpr (std::same_as<Value, RealValued>) {
    if (m < 0) {
      return static_cast<typename Grid::Real>(MinusOneToPower(m)) *
             std::conj(expansion[l, -m]);
    }
  }
  return expansion[l, m];
}

// Raise the upper index by one. The result is an expansion at N + 1 over the
// same degrees, less the one the raised field cannot carry.
template <std::ptrdiff_t N, AngularGrid Grid, RealOrComplexValued Value>
auto Raise(const SpinExpansion<N, Grid, Value>& expansion) {
  using Real = typename Grid::Real;
  auto raised =
      SpinExpansion<N + 1, Grid, ComplexValued>(expansion.Grid(),
                                                expansion.MaxDegree());
  for (auto l : raised.Degrees()) {
    for (auto m : raised.Orders(l)) {
      raised[l, m] = EthDetails::RaisingFactor<Real>(l, N) *
                     Coefficient(expansion, l, m);
    }
  }
  return raised;
}

// Lower it by one.
template <std::ptrdiff_t N, AngularGrid Grid, RealOrComplexValued Value>
auto Lower(const SpinExpansion<N, Grid, Value>& expansion) {
  using Real = typename Grid::Real;
  auto lowered =
      SpinExpansion<N - 1, Grid, ComplexValued>(expansion.Grid(),
                                                expansion.MaxDegree());
  for (auto l : lowered.Degrees()) {
    for (auto m : lowered.Orders(l)) {
      lowered[l, m] = EthDetails::LoweringFactor<Real>(l, N) *
                      Coefficient(expansion, l, m);
    }
  }
  return lowered;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_ETH_GUARD_H
