#ifndef GSH_TRANS_INTRINSIC_DERIVATIVE_GUARD_H
#define GSH_TRANS_INTRINSIC_DERIVATIVE_GUARD_H

#include <cstddef>
#include <utility>

#include "../Concepts.h"
#include "../Tensor/MultiIndex.h"
#include "ContravariantDerivative.h"
#include "Eth.h"
#include "TensorExpansion.h"

namespace GSHTrans {

// The intrinsic covariant derivative of a tangential tensor: the Levi-Civita
// connection of the metric the sphere induces, and the operator surface
// differential geometry means by differentiation.
//
// There are two derivatives here and they must be kept apart.
//
// SurfaceGradient is D&T's grad_1, the angular part of the three-dimensional
// gradient. It differentiates the tensor *and* its basis, and the basis leaves
// the tangent plane -- so it does not close on a tangential tensor, and it
// takes no operand from that bundle at all.
//
// This is the other one. The two are related by the Gauss formula: grad_1 is D
// plus a term in the second fundamental form, which for the unit sphere is
// algebraic rather than differential. Splitting grad_1 of a tangential T by
// whether the inherited slot is radial,
//
//   (grad_1 T)^{sigma a_1...a_q}      = Omega^{-+N}_l T^{a_1...a_q}
//   (grad_1 T)^{sigma a_1...0_j...a_q} = -T^{a_1...sigma_j...a_q}
//
// the first line is the tangential block and *is* D; the second is the
// extrinsic curvature, which for a unit sphere is exactly minus the field with
// the slot replaced.
//
// So the first line has no connection terms at all. Every shift a_i + sigma
// either leaves {-1, 0, +1} or lands on a radial slot, and a tangential tensor
// has no radial slot to land on -- which is why the operator is a pure Omega
// multiplication and why it is closed. That makes it eth applied component by
// component, up to the sqrt(2) that is eth's own normalisation.
//
// **The normalisation is forced, not chosen.** Gauss says D is the tangential
// block of grad_1, exactly, so it carries D&T's Omega. Adopting eth's instead
// would make this operator differ by sqrt(2) from the operator it is defined
// to be. The sqrt(2) stays where it belongs, as a fact about eth.

namespace IntrinsicDetails {

using Int = std::ptrdiff_t;

// One component of the result: the operand's component with the leading slot
// dropped, scaled by Omega at every degree.
//
// The Omega term draws on the *operand's* upper index, as it does in the
// surface gradient, and the upper sign goes with sigma = -1.
template <auto Indices, typename Result, typename Expansion>
void FillComponent(Result& result, const Expansion& operand) {
  using Real = typename Result::Real;

  constexpr auto Rank = static_cast<Int>(Indices.size()) - 1;
  constexpr auto sigma = Indices[0];
  constexpr auto source = ContravariantDetails::DropFirst<Rank>(Indices);
  constexpr auto sourceN =
      MultiIndex<Rank, TangentialSlots>(source).UpperIndex();

  // The result is tangential, of rank at least one, and carries no symmetry.
  // Negation therefore has no fixed point among its multi-indices -- the
  // all-zero one does not exist over this alphabet -- and no permutation can
  // pin anything either, so no component of the result is constrained. That
  // is where knowing the orbits removes a branch rather than a buffer: the
  // surface gradient needs a pinned-imaginary case here and this does not.
  constexpr auto flat = MultiIndex<Rank + 1, TangentialSlots>(Indices).Flat();
  static_assert(Result::Orbits.constraint[flat] == ComponentConstraint::None,
                "A tangential tensor with no symmetry has no pinned "
                "components, so this operator never has to store one");

  auto block = ContravariantDetails::ComponentOfArray<Indices>(
      result, std::make_index_sequence<static_cast<std::size_t>(Rank + 1)>{});

  for (auto l : block.Degrees()) {
    for (auto m : block.Orders(l)) {
      block[l, m] =
          EthDetails::Omega<Real>(l, sigma == -1 ? sourceN : -sourceN) *
          ContravariantDetails::Term<source>(operand, l, m);
    }
  }
}

}  // namespace IntrinsicDetails

// D T for a tangential tensor expansion: rank q in, rank q + 1 out, and the
// result is tangential too -- which is the whole point of the operator having
// a name of its own.
//
// The result has no symmetry, for the reason the surface gradient's has none:
// D of a symmetric tensor is not symmetric in the new slot against the old
// ones, and inferring a symmetry from an operator is the problem Materialise
// declined to solve. It keeps the operand's reality, since D of a real tensor
// is real.
template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid>
auto IntrinsicDerivative(const TensorExpansion<Rank, Symmetry, Reality, Grid,
                                               TangentialSlots>& operand) {
  using Result = TensorExpansion<Rank + 1, NoSymmetry<Rank + 1>, Reality, Grid,
                                 TangentialSlots>;

  auto result = Result(operand.Grid(), operand.MaxDegree());
  [&]<std::size_t... Slot>(std::index_sequence<Slot...>) {
    (
        [&] {
          constexpr auto flat = Result::ComponentLayout.flatOfSlot[Slot];
          constexpr auto indices =
              MultiIndex<Rank + 1, TangentialSlots>::FromFlat(flat).Slots();
          IntrinsicDetails::FillComponent<indices>(result, operand);
        }(),
        ...);
  }(std::make_index_sequence<static_cast<std::size_t>(
        Result::StoredComponents)>{});
  return result;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_INTRINSIC_DERIVATIVE_GUARD_H
