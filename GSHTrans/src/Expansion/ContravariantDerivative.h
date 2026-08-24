#ifndef GSH_TRANS_CONTRAVARIANT_DERIVATIVE_GUARD_H
#define GSH_TRANS_CONTRAVARIANT_DERIVATIVE_GUARD_H

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <utility>

#include "../Concepts.h"
#include "../Tensor/MultiIndex.h"
#include "TensorExpansion.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                        The contravariant derivative                       //
//--------------------------------------------------------------------------//

// The surface gradient of a tensor field, in the formalism of Phinney &
// Burridge and of Dahlen & Tromp (C.151)-(C.153).
//
// This is the derivative a tensor field actually wants, and it is *not* eth
// applied component by component. For a scalar the two agree up to a factor
// of sqrt(2) -- the normalisation of e_{+-} -- but for rank one and above the
// operator is not a multiplication at all: it also subtracts the components
// with one slot index shifted, once per slot, because the canonical basis
// vectors vary over the sphere and differentiating a tensor field
// differentiates its basis too:
//
//   [+-d_theta + i(sin theta)^{-1} d_phi] e_a = a cot(theta) e_a
//                                               - sqrt(2) e_{a +- 1}.
//
// With Omega^{+-N}_l = sqrt((l +- N)(l -+ N + 1) / 2), the operator takes the
// coefficients of a rank-q tensor to those of the rank-(q+1) tensor grad T:
//
//   (grad T)^{sigma a_1...a_q}_{lm} = Omega^{-+N}_l T^{a_1...a_q}_{lm}
//                                     - sum_i T^{a_1...(a_i + sigma)...a_q}_{lm}
//
// for sigma = +-1, the upper sign going with sigma = -1. The operator
// *prepends* a slot, and the result carries upper index sigma + N, which is
// eq:N applied to the lengthened multi-index.
//
// The radial part is not here. D&T's full gradient is
// [e_0 d_r + r^{-1} grad_1], and the sigma = 0 component is dT/dr, which an
// angular library cannot supply; the *surface* gradient grad_1 has no e_0
// component at all (C.145), so this operator is complete without it and the
// radial derivative is the user's to supply. What that means for the result
// is that its sigma = 0 block is stored and zero -- a rank-(q+1) tensor whose
// e_0 slot vanishes is an ordinary tensor, and keeping it as one leaves the
// result composable with contraction, further gradients and the transform.

namespace EthDetails {

// Omega^{s}_l = sqrt((l + s)(l - s + 1) / 2), so that Omega^{+N} is
// Omega(l, N) and Omega^{-N} is Omega(l, -N).
template <RealFloatingPoint Real>
Real Omega(std::ptrdiff_t l, std::ptrdiff_t s) {
  const auto a = static_cast<Real>(l + s);
  const auto b = static_cast<Real>(l - s + 1);
  return std::sqrt(a * b / 2);
}

}  // namespace EthDetails

namespace ContravariantDetails {

using Int = std::ptrdiff_t;

// The multi-index of the operand that a component of the result draws on:
// the result's indices with the leading one dropped.
template <Int Rank, std::size_t Size>
constexpr auto DropFirst(const std::array<Int, Size>& indices) {
  auto rest = std::array<Int, Rank>{};
  for (auto i = Int{0}; i < Rank; i++) rest[i] = indices[i + 1];
  return rest;
}

// The same, with slot `which` shifted by `by`. Returns whether it stayed
// inside {-1, 0, 1}: a shifted index that leaves the range names no
// component, and D&T stipulate that its coefficient is zero.
template <Int Rank>
constexpr auto Shifted(const std::array<Int, Rank>& indices, Int which,
                       Int by) {
  auto shifted = indices;
  shifted[which] += by;
  const auto inside = shifted[which] >= -1 && shifted[which] <= 1;
  return std::pair(shifted, inside);
}

// Read one term, by expanding a compile-time multi-index into the pack that
// TensorExpansion::Coefficient wants.
template <auto Indices, typename Expansion, std::size_t... I>
auto TermOf(const Expansion& operand, Int l, Int m, std::index_sequence<I...>) {
  return operand.template Coefficient<Indices[I]...>(l, m);
}

template <auto Indices, typename Expansion>
auto Term(const Expansion& operand, Int l, Int m) {
  return TermOf<Indices>(operand, l, m,
                         std::make_index_sequence<Indices.size()>{});
}

// Reach a component of the result, expanding a compile-time multi-index into
// the pack the accessor wants.
template <auto Indices, typename Result, std::size_t... I>
auto ComponentOfArray(Result& result, std::index_sequence<I...>) {
  return result.template Component<Indices[I]...>();
}

// One component of the result, over every (l, m) it carries.
//
// Only the components the result *stores* are filled: on a real tensor the
// rest follow from eq:reality, and the gradient of a real tensor is real.
//
// `scale` multiplies every coefficient. It is 1 for the surface gradient and
// r^{-1} for the angular part of the full three-dimensional one, which is the
// only difference between them: D&T's gradient is [e_0 d_r + r^{-1} grad_1],
// so the same formula serves both and the mathematics is written once.
template <auto Indices, typename Result, typename Expansion>
void FillComponent(Result& result, const Expansion& operand,
                   typename Result::Real scale = 1) {
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.
  using Real = typename Result::Real;  ///< The precision.
  using Complex = typename Result::Complex;  ///< `std::complex` over the precision.

  constexpr auto Rank = static_cast<Int>(Indices.size()) - 1;
  constexpr auto sigma = Indices[0];

  // grad_1 has no e_0 component (D&T C.145). The block was built zero and
  // stays that way, which keeps the result an ordinary rank-(q+1) tensor.
  if constexpr (sigma == 0) {
    (void)result;
    (void)operand;
    (void)scale;
    return;
  } else {
    constexpr auto source = DropFirst<Rank>(Indices);
    constexpr auto sourceN = MultiIndex<Rank>(source).UpperIndex();
    constexpr auto flat = MultiIndex<Rank + 1>(Indices).Flat();
    constexpr auto constraint = Result::Orbits.constraint[flat];

    auto block = ComponentOfArray<Indices>(
        result, std::make_index_sequence<static_cast<std::size_t>(Rank + 1)>{});

    for (auto l : block.Degrees()) {
      for (auto m : block.Orders(l)) {
        // The Omega term draws on the operand's own upper index; the shift
        // terms are already at the result's.
        auto value =
            EthDetails::Omega<Real>(l, sigma == -1 ? sourceN : -sourceN) *
            Term<source>(operand, l, m);

        [&]<std::size_t... Slot>(std::index_sequence<Slot...>) {
          (
              [&] {
                constexpr auto shifted = Shifted<Rank>(source, Slot, sigma);
                if constexpr (shifted.second) {
                  value -= Term<shifted.first>(operand, l, m);
                }
              }(),
              ...);
        }(std::make_index_sequence<static_cast<std::size_t>(Rank)>{});

        // A component the reality condition pins as purely imaginary stores
        // the real field whose i-multiple it is, so what goes into the block
        // is the value divided by i. A pinned-real one stores itself.
        //
        // The imaginary branch is unreachable as things stand and is here to
        // stay correct if that changes: the result always has NoSymmetry, and
        // under negation alone the only self-paired multi-index is the
        // all-zero one, which is pinned real. An imaginary pinning needs a
        // permutation sign, so it can only arise in an *operand*.
        value *= scale;
        if constexpr (constraint == ComponentConstraint::Imaginary) {
          block[l, m] = Complex{0, -1} * value;
        } else {
          block[l, m] = value;
        }
      }
    }
  }
}

}  // namespace ContravariantDetails

// The surface gradient of a tensor expansion: rank q in, rank q + 1 out.
//
// The result has no symmetry, because grad_1 of a symmetric tensor is not
// symmetric in the new slot against the old ones, and inferring a symmetry
// from an operator is the problem Materialise declined to solve. It keeps the
// operand's reality, since the gradient of a real tensor is real.
//
// And it takes a general tensor, deliberately. grad_1 is an operator on
// sections of the general bundle: it moves slots between e_0 and e_+-, so it
// does not close on a tangential tensor and no signature over that alphabet
// would be honest. A tangential operand is embedded first and then
// differentiated, in that order and visibly, which is why the deduction below
// simply does not match one. What
// *is* closed on a tangential tensor is the intrinsic derivative, and that
// has a name of its own.
template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid>
auto SurfaceGradient(
    const TensorExpansion<Rank, Symmetry, Reality, Grid>& operand) {
  using Result =
      TensorExpansion<Rank + 1, NoSymmetry<Rank + 1>, Reality, Grid>;

  auto result = Result(operand.Grid(), operand.MaxDegree());
  [&]<std::size_t... Slot>(std::index_sequence<Slot...>) {
    (
        [&] {
          constexpr auto flat = Result::ComponentLayout.flatOfSlot[Slot];
          constexpr auto indices =
              MultiIndex<Rank + 1>::FromFlat(flat).Slots();
          ContravariantDetails::FillComponent<indices>(result, operand);
        }(),
        ...);
  }(std::make_index_sequence<static_cast<std::size_t>(
        Result::StoredComponents)>{});
  return result;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_CONTRAVARIANT_DERIVATIVE_GUARD_H
