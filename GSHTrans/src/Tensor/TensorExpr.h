#ifndef GSH_TRANS_TENSOR_EXPR_GUARD_H
#define GSH_TRANS_TENSOR_EXPR_GUARD_H

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "../Concepts.h"
#include "../SpinField/SpinWeighted.h"
#include "MultiIndex.h"
#include "TensorField.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                            The tensor-level node                          //
//--------------------------------------------------------------------------//

// A tensor expression: a rank, a grid, and a component accessor.
//
// The whole of the tensor algebra is that an operation on tensors is an
// operation on that accessor. Every node below returns a *phase-1* node from
// Component<...>(), so the index algebra of upper indices is already enforced
// and there is no second expression system to keep consistent with the first.
// What this layer has to get right is only *which* components it asks for.
//
// The concept cannot require Component<...>() itself: the multi-index is a
// pack whose admissible values depend on the tensor's symmetry -- the all-zero
// component of an antisymmetric tensor does not exist -- so there is no single
// instantiation every tensor must provide. It requires the parts that are
// universal, and `Represents` is what a caller consults before asking.
template <typename T>
concept TensorExpr = requires {
  requires std::same_as<std::remove_cv_t<decltype(T::Rank)>, std::ptrdiff_t>;
  requires T::Rank >= 0;
  typename T::GridType;
  requires AngularGrid<typename T::GridType>;
  { T::template Represents<> } -> std::convertible_to<bool>;
} and requires(const T& tensor) {
  { tensor.Grid() } -> std::convertible_to<const typename T::GridType&>;
};

template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid, TensorLayout Layout>
struct IsTerminalTrait<TensorField<Rank, Symmetry, Reality, Grid, Layout>>
    : std::true_type {};

namespace TensorDetails {

using Int = std::ptrdiff_t;

// Ask an operand for the component named by a compile-time array, by
// expanding it back into a template pack.
//
// This is the piece of machinery every node below needs, because a tensor
// operation computes a multi-index and then has to *use* it as template
// arguments. An index array is a structural type, so it can be a template
// parameter; turning it back into a pack takes an index_sequence.
template <std::array<Int, 0> Indices, typename Operand>
auto ComponentOf(const Operand& operand, std::index_sequence<>) {
  return operand.template Component<>();
}

template <auto Indices, typename Operand, std::size_t... I>
auto ComponentOf(const Operand& operand, std::index_sequence<I...>) {
  return operand.template Component<Indices[I]...>();
}

template <auto Indices, typename Operand>
auto ComponentOf(const Operand& operand) {
  return ComponentOf<Indices>(
      operand, std::make_index_sequence<Indices.size()>{});
}

// Whether an operand represents the component named by an array.
template <auto Indices, typename Operand, std::size_t... I>
constexpr bool RepresentsImpl(std::index_sequence<I...>) {
  return Operand::template Represents<Indices[I]...>;
}

template <auto Indices, typename Operand>
constexpr bool Represents() {
  return RepresentsImpl<Indices, Operand>(
      std::make_index_sequence<Indices.size()>{});
}

}  // namespace TensorDetails

//--------------------------------------------------------------------------//
//                            Permuting the slots                            //
//--------------------------------------------------------------------------//

// A tensor with its slots relabelled: (Permute T)^{alpha} = T^{pi(alpha)}.
//
// No arithmetic at all -- the component accessor forwards to the operand with
// the indices reordered -- which makes this the step that pins the machinery
// the rest of the layer uses.
//
// It permutes *tensor slots* and not grid points, so it does not touch the
// pointwise-and-index-preserving invariant that phase 1's aliasing theorem
// rests on. That is worth stating because a re-indexing view is exactly the
// thing that would break it, and this is the closest the library comes to one.
template <auto Image, typename Operand>
class PermuteNode {
 public:
  using Int = std::ptrdiff_t;
  using OperandType = std::remove_cvref_t<Operand>;

  static constexpr Int Rank = OperandType::Rank;
  using GridType = typename OperandType::GridType;

  static_assert(Image.size() == static_cast<std::size_t>(Rank),
                "A slot permutation needs one image per tensor slot");

  explicit PermuteNode(Operand&& operand)
      : _operand{std::forward<Operand>(operand)} {}

  const GridType& Grid() const { return _operand.Grid(); }

  // Slot i of this tensor is slot Image[i] of the operand, which is the same
  // convention MultiIndex::Permuted uses.
  template <Int... Alphas>
  static constexpr auto Source =
      MultiIndex<Rank>(std::array<Int, Rank>{Alphas...}).Permuted(Image).Slots();

  template <Int... Alphas>
  static constexpr bool Represents =
      sizeof...(Alphas) == static_cast<std::size_t>(Rank) and
      TensorDetails::Represents<Source<Alphas...>, OperandType>();

  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const {
    return TensorDetails::ComponentOf<Source<Alphas...>>(_operand);
  }

 private:
  OperandStorage<Operand> _operand;
};

// Relabel a tensor's slots. The image is given as a std::array, so that
// Permute<std::array{1, 0}>(t) reads as "slot 0 of the result is slot 1 of t".
template <auto Image, typename T>
requires TensorExpr<std::remove_cvref_t<T>>
auto Permute(T&& tensor) {
  return PermuteNode<Image, T>(std::forward<T>(tensor));
}

// The rank-2 case, which is the one with a name.
template <typename T>
requires TensorExpr<std::remove_cvref_t<T>> &&
         (std::remove_cvref_t<T>::Rank == 2)
auto Transpose(T&& tensor) {
  return Permute<std::array<std::ptrdiff_t, 2>{1, 0}>(std::forward<T>(tensor));
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_TENSOR_EXPR_GUARD_H
