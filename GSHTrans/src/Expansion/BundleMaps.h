#ifndef GSH_TRANS_EXPANSION_BUNDLE_MAPS_GUARD_H
#define GSH_TRANS_EXPANSION_BUNDLE_MAPS_GUARD_H

#include <complex>
#include <concepts>
#include <cstddef>
#include <utility>

#include "../Concepts.h"
#include "../Tensor/MultiIndex.h"
#include "../Tensor/Orbits.h"
#include "ContravariantDerivative.h"
#include "TensorExpansion.h"

namespace GSHTrans {

// The same two maps as GSHTrans/src/Tensor/BundleMaps.h, in the spectral
// domain: Embed widens a tangential expansion into the general bundle, and
// Tangential projects one back.
//
// **These copy, where the spatial pair does not, and that is a decision rather
// than an oversight** (field-algebra-plan.md section 18.2 [D13]). Three things
// make it the right one. There are no lazy nodes on this side at all, and by
// design -- section 15.4 records that a derived component is not offered
// spectrally because deriving one reverses the order index. Making Embed lazy
// would mean re-signing SurfaceGradient against a concept that could not carry
// the one member distinguishing this type from the layered one, since
// Coefficient is a template whose pack depends on the rank and a concept
// cannot require it. And the saving would be small: SurfaceGradient allocates
// a rank-(q+1) result whatever its operand is, which is larger than the
// embedded operand it would have avoided.
//
// What the pair is for is the identity that closes the split:
//
//   Tangential(SurfaceGradient(Embed(T))) == IntrinsicDerivative(T)
//
// which says the intrinsic derivative is the tangential block of the ambient
// one, and is what a caller writes when they want the ambient gradient of a
// tangential tensor at all.

namespace BundleDetails {

// Whether a multi-index given as a compile-time array has a radial slot.
template <auto Indices>
constexpr bool AnyRadialIn() {
  for (auto alpha : Indices) {
    if (alpha == 0) return true;
  }
  return false;
}

// Copy one component's block across, at every degree and order it carries.
//
// The two bundles agree about this component's orbit -- permutation and
// negation both preserve whether a multi-index has a radial slot, so an orbit
// of tangential indices contains only tangential indices and has the same
// representative and the same constraint on either side. So the value read is
// the value to store, up to the one adjustment a pinned component needs.
template <auto Indices, typename Result, typename Operand>
void CopyComponent(Result& result, const Operand& operand) {
  using Complex = typename Result::Complex;  ///< `std::complex` over the precision.
  constexpr auto Rank = static_cast<std::ptrdiff_t>(Indices.size());
  constexpr auto flat =
      MultiIndex<Rank, typename Result::SlotSet>(Indices).Flat();
  constexpr auto constraint = Result::Orbits.constraint[flat];

  auto block = ContravariantDetails::ComponentOfArray<Indices>(
      result, std::make_index_sequence<static_cast<std::size_t>(Rank)>{});

  for (auto l : block.Degrees()) {
    for (auto m : block.Orders(l)) {
      const auto value = ContravariantDetails::Term<Indices>(operand, l, m);
      // A component the reality condition pins as purely imaginary stores the
      // real field whose i-multiple it is, exactly as in the surface gradient.
      if constexpr (constraint == ComponentConstraint::Imaginary) {
        block[l, m] = Complex{0, -1} * value;
      } else {
        block[l, m] = value;
      }
    }
  }
}

// Fill every stored component of the result that the operand can supply,
// leaving the rest as the buffer was built: zero.
template <typename Result, typename Operand, bool SkipRadial>
void FillFrom(Result& result, const Operand& operand) {
  [&]<std::size_t... Slot>(std::index_sequence<Slot...>) {
    (
        [&] {
          constexpr auto flat = Result::ComponentLayout.flatOfSlot[Slot];
          constexpr auto indices =
              MultiIndex<Result::Rank,
                         typename Result::SlotSet>::FromFlat(flat)
                  .Slots();
          if constexpr (!SkipRadial || !AnyRadialIn<indices>()) {
            CopyComponent<indices>(result, operand);
          }
        }(),
        ...);
  }(std::make_index_sequence<static_cast<std::size_t>(
        Result::StoredComponents)>{});
}

}  // namespace BundleDetails

// Tangential into the general. Rank, symmetry, reality and grid unchanged --
// a permutation preserves the slot sum and maps radial slots to radial slots,
// so an embedded symmetric tensor is symmetric and an embedded real one is
// real. Every component with a radial slot is left at zero, which is what the
// embedded tensor has there.
template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid>
auto Embed(const TensorExpansion<Rank, Symmetry, Reality, Grid,
                                 TangentialSlots>& operand) {
  using Operand =
      TensorExpansion<Rank, Symmetry, Reality, Grid, TangentialSlots>;
  using Result = TensorExpansion<Rank, Symmetry, Reality, Grid, AllSlots>;
  auto result = Result(operand.Grid(), operand.MaxDegree());
  BundleDetails::FillFrom<Result, Operand, true>(result, operand);
  return result;
}

// General onto the tangential. Every component the result stores is one the
// operand has, so nothing is skipped on this side.
template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid>
auto Tangential(
    const TensorExpansion<Rank, Symmetry, Reality, Grid, AllSlots>& operand) {
  using Operand = TensorExpansion<Rank, Symmetry, Reality, Grid, AllSlots>;
  using Result =
      TensorExpansion<Rank, Symmetry, Reality, Grid, TangentialSlots>;
  auto result = Result(operand.Grid(), operand.MaxDegree());
  BundleDetails::FillFrom<Result, Operand, false>(result, operand);
  return result;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_EXPANSION_BUNDLE_MAPS_GUARD_H
