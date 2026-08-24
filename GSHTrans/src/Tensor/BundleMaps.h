#ifndef GSH_TRANS_BUNDLE_MAPS_GUARD_H
#define GSH_TRANS_BUNDLE_MAPS_GUARD_H

#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "../Concepts.h"
#include "MultiIndex.h"
#include "TensorExpr.h"

namespace GSHTrans {

// The two maps between the tangential bundle and the general one.
//
// A tangential tensor lives in another bundle, so it is another object
// (field-algebra-plan.md section 18.2 [D8]). What connects the two is the
// inclusion of the tangent plane in the ambient space, and its adjoint:
//
//   Embed(T)       tangential -> general, the components with a radial slot
//                  being ones the embedded tensor does not have
//   Tangential(T)  general -> tangential, dropping those same components
//
// Neither is implicit. A product of operands from different alphabets does not
// compile ([D10]) and grad_1 takes no tangential operand ([D9]); both are
// written by embedding at the call site, where it can be read.
//
// Both are lazy, and neither owns storage. That is what makes crossing bundles
// cost a node rather than a buffer: Embed allocates nothing and writes no
// zeros, because a component with a radial slot is one it declines to
// represent -- the same signal an antisymmetric tensor's diagonal already
// emits, which every compile-time traversal, Contract and Materialise already
// handle. See [D13] for why the spectral pair is not built this way.

namespace BundleDetails {

using Int = std::ptrdiff_t;

// Whether a multi-index given as a pack has a radial slot.
template <Int... Alphas>
constexpr bool AnyRadial() {
  return ((Alphas == 0) or ...);
}

}  // namespace BundleDetails

//--------------------------------------------------------------------------//
//                        Tangential into the general                        //
//--------------------------------------------------------------------------//

// The inclusion. Rank and grid unchanged; the alphabet widens.
//
// A component with a radial slot is not represented rather than zero: the
// embedded tensor genuinely has no content there, and saying so is both
// cheaper than storing a zero and more informative to a traversal, which can
// then skip the component rather than evaluate it. It is the same statement
// the orbit table makes about a component an antisymmetry annihilates.
template <typename Operand>
class EmbedNode {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.
  using OperandType = std::remove_cvref_t<Operand>;

  static constexpr Int Rank = OperandType::Rank;
  using GridType = typename OperandType::GridType;  ///< The angular grid this is defined on.
  using SlotSet = AllSlots;  ///< The alphabet the slots are drawn from.

  static_assert(std::same_as<typename OperandType::SlotSet, TangentialSlots>,
                "Embed widens a tangential tensor; a general one is already "
                "where it is going");

  explicit EmbedNode(Operand&& operand)
      : _operand{std::forward<Operand>(operand)} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _operand.Grid(); }

  template <Int... Alphas>
  static constexpr bool RepresentsFn() {
    if constexpr (sizeof...(Alphas) != static_cast<std::size_t>(Rank)) {
      return false;
    } else if constexpr (!AreSlotLetters<SlotSet, Alphas...>()) {
      return false;
    } else if constexpr (BundleDetails::AnyRadial<Alphas...>()) {
      return false;
    } else {
      return OperandType::template Represents<Alphas...>;
    }
  }

  template <Int... Alphas>
  static constexpr bool Represents = RepresentsFn<Alphas...>();

  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const {
    return _operand.template Component<Alphas...>();
  }

 private:
  OperandStorage<Operand> _operand;
};

template <typename T>
requires TensorExpr<std::remove_cvref_t<T>> &&
         std::same_as<typename std::remove_cvref_t<T>::SlotSet,
                      TangentialSlots>
auto Embed(T&& tensor) {
  return EmbedNode<T>(std::forward<T>(tensor));
}

//--------------------------------------------------------------------------//
//                        General onto the tangential                        //
//--------------------------------------------------------------------------//

// The projection, which is how a caller comes back down: the tangential part
// of a three-dimensional strain, or the tangential block of a gradient.
//
// A pure relabelling like Permute, and for the same reason harmless to phase
// 1's aliasing theorem: it selects tensor slots and does not touch grid
// points.
template <typename Operand>
class TangentialNode {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.
  using OperandType = std::remove_cvref_t<Operand>;

  static constexpr Int Rank = OperandType::Rank;
  using GridType = typename OperandType::GridType;  ///< The angular grid this is defined on.
  using SlotSet = TangentialSlots;  ///< The alphabet the slots are drawn from.

  static_assert(std::same_as<typename OperandType::SlotSet, AllSlots>,
                "Tangential projects a general tensor; a tangential one is "
                "already where it is going");

  explicit TangentialNode(Operand&& operand)
      : _operand{std::forward<Operand>(operand)} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _operand.Grid(); }

  template <Int... Alphas>
  static constexpr bool RepresentsFn() {
    if constexpr (sizeof...(Alphas) != static_cast<std::size_t>(Rank)) {
      return false;
    } else if constexpr (!AreSlotLetters<SlotSet, Alphas...>()) {
      return false;
    } else {
      return OperandType::template Represents<Alphas...>;
    }
  }

  template <Int... Alphas>
  static constexpr bool Represents = RepresentsFn<Alphas...>();

  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const {
    return _operand.template Component<Alphas...>();
  }

 private:
  OperandStorage<Operand> _operand;
};

template <typename T>
requires TensorExpr<std::remove_cvref_t<T>> &&
         std::same_as<typename std::remove_cvref_t<T>::SlotSet, AllSlots>
auto Tangential(T&& tensor) {
  return TangentialNode<T>(std::forward<T>(tensor));
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_BUNDLE_MAPS_GUARD_H
