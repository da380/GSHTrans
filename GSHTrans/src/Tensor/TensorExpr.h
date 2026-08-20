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

// The first N and last N entries of an index array, for splitting a product's
// multi-index between its two operands.
template <Int N, std::size_t Size>
constexpr auto Head(const std::array<Int, Size>& indices) {
  auto head = std::array<Int, N>{};
  for (auto i = Int{0}; i < N; i++) head[i] = indices[i];
  return head;
}

template <Int N, std::size_t Size>
constexpr auto Tail(const std::array<Int, Size>& indices) {
  auto tail = std::array<Int, N>{};
  for (auto i = Int{0}; i < N; i++) {
    tail[i] = indices[static_cast<Int>(Size) - N + i];
  }
  return tail;
}

// The operand multi-index of a contraction: the surviving indices in order,
// with `alpha` at slot J and `-alpha` at slot K.
template <Int J, Int K, Int Rank, std::size_t Size>
constexpr auto Insert(const std::array<Int, Size>& surviving, Int alpha) {
  auto full = std::array<Int, Rank>{};
  full[J] = alpha;
  full[K] = -alpha;
  auto next = Int{0};
  for (auto slot = Int{0}; slot < Rank; slot++) {
    if (slot == J || slot == K) continue;
    full[slot] = surviving[next++];
  }
  return full;
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


//--------------------------------------------------------------------------//
//                             The tensor product                            //
//--------------------------------------------------------------------------//

// (S tensor T)^{alpha beta} = S^{alpha} T^{beta}, of rank p + q.
//
// The component is a phase-1 product, so its upper index is the sum of the
// two operands' -- which is eq:N applied to the concatenated multi-index, and
// the theory note says so in as many words. Nothing here has to arrange that;
// it is what "upper indices add under pointwise multiplication" means.
template <typename LeftOperand, typename RightOperand>
class TensorProductNode {
 public:
  using Int = std::ptrdiff_t;
  using Left = std::remove_cvref_t<LeftOperand>;
  using Right = std::remove_cvref_t<RightOperand>;

  static constexpr Int Rank = Left::Rank + Right::Rank;
  using GridType = typename Left::GridType;

  static_assert(std::same_as<GridType, typename Right::GridType>,
                "A tensor product needs both operands on the same kind of "
                "grid");

  TensorProductNode(LeftOperand&& left, RightOperand&& right)
      : _left{std::forward<LeftOperand>(left)},
        _right{std::forward<RightOperand>(right)} {
    if (_left.Grid().Identity() != _right.Grid().Identity()) {
      throw std::invalid_argument(
          "A tensor product needs both operands on the same grid");
    }
  }

  const GridType& Grid() const { return _left.Grid(); }

  // The multi-index splits: the first p slots name the left operand's
  // component and the last q the right's.
  template <Int... Alphas>
  static constexpr auto Indices = std::array<Int, Rank>{Alphas...};

  template <Int... Alphas>
  static constexpr auto LeftIndices = TensorDetails::Head<Left::Rank>(
      Indices<Alphas...>);

  template <Int... Alphas>
  static constexpr auto RightIndices = TensorDetails::Tail<Right::Rank>(
      Indices<Alphas...>);

  template <Int... Alphas>
  static constexpr bool RepresentsFn() {
    if constexpr (sizeof...(Alphas) != static_cast<std::size_t>(Rank)) {
      return false;
    } else {
      return TensorDetails::Represents<LeftIndices<Alphas...>, Left>() and
             TensorDetails::Represents<RightIndices<Alphas...>, Right>();
    }
  }

  template <Int... Alphas>
  static constexpr bool Represents = RepresentsFn<Alphas...>();

  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const {
    return TensorDetails::ComponentOf<LeftIndices<Alphas...>>(_left) *
           TensorDetails::ComponentOf<RightIndices<Alphas...>>(_right);
  }

 private:
  OperandStorage<LeftOperand> _left;
  OperandStorage<RightOperand> _right;
};

template <typename L, typename R>
requires TensorExpr<std::remove_cvref_t<L>> && TensorExpr<std::remove_cvref_t<R>>
auto TensorProduct(L&& left, R&& right) {
  return TensorProductNode<L, R>(std::forward<L>(left), std::forward<R>(right));
}

//--------------------------------------------------------------------------//
//                                Contraction                                //
//--------------------------------------------------------------------------//

// Contract slots J and K against the metric of eq:metric,
//
//   (tr_{JK} T)^{...} = sum_{ab} g_{ab} T^{... a ... b ...}
//                     = sum_a (-1)^a T^{... a ... -a ...},
//
// which for the three values of a is -T^{-+} + T^{00} - T^{+-} in the
// contracted pair.
//
// The result lands at the right upper index by construction rather than by
// arrangement: the contracted pair contributes a + (-a) = 0 whatever a is, so
// all three terms carry the same upper index and phase 1's Equal rule admits
// their sum. A contraction that paired its slots wrongly would not compile.
template <std::ptrdiff_t J, std::ptrdiff_t K, typename Operand>
class ContractionNode {
 public:
  using Int = std::ptrdiff_t;
  using OperandType = std::remove_cvref_t<Operand>;

  static constexpr Int Rank = OperandType::Rank - 2;
  using GridType = typename OperandType::GridType;

  static_assert(J != K, "A contraction needs two different slots");
  static_assert(J >= 0 && K >= 0 && J < OperandType::Rank &&
                    K < OperandType::Rank,
                "A contracted slot must be one the tensor has");

  explicit ContractionNode(Operand&& operand)
      : _operand{std::forward<Operand>(operand)} {}

  const GridType& Grid() const { return _operand.Grid(); }

  // The operand's multi-index: the surviving slots in order, with a inserted
  // at J and -a at K.
  template <Int Alpha, Int... Alphas>
  static constexpr auto Source = TensorDetails::Insert<J, K, OperandType::Rank>(
      std::array<Int, Rank>{Alphas...}, Alpha);

  template <Int... Alphas>
  static constexpr bool RepresentsFn() {
    if constexpr (sizeof...(Alphas) != static_cast<std::size_t>(Rank)) {
      return false;
    } else {
      // Every term of the sum has to exist, since they are added. A tensor
      // with a vanishing component in the middle of a contraction is not
      // something this layer tries to be clever about.
      return TensorDetails::Represents<Source<-1, Alphas...>, OperandType>() and
             TensorDetails::Represents<Source<0, Alphas...>, OperandType>() and
             TensorDetails::Represents<Source<1, Alphas...>, OperandType>();
    }
  }

  template <Int... Alphas>
  static constexpr bool Represents = RepresentsFn<Alphas...>();

  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const {
    return -TensorDetails::ComponentOf<Source<-1, Alphas...>>(_operand) +
           TensorDetails::ComponentOf<Source<0, Alphas...>>(_operand) -
           TensorDetails::ComponentOf<Source<1, Alphas...>>(_operand);
  }

 private:
  OperandStorage<Operand> _operand;
};

template <std::ptrdiff_t J, std::ptrdiff_t K, typename T>
requires TensorExpr<std::remove_cvref_t<T>>
auto Contract(T&& tensor) {
  return ContractionNode<J, K, T>(std::forward<T>(tensor));
}

// The trace of a rank-2 tensor, which is the contraction with a name.
template <typename T>
requires TensorExpr<std::remove_cvref_t<T>> &&
         (std::remove_cvref_t<T>::Rank == 2)
auto Trace(T&& tensor) {
  return Contract<0, 1>(std::forward<T>(tensor)).template Component<>();
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_TENSOR_EXPR_GUARD_H
