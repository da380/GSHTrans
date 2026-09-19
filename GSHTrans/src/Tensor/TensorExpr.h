#ifndef GSH_TRANS_TENSOR_EXPR_GUARD_H
#define GSH_TRANS_TENSOR_EXPR_GUARD_H

#include <array>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <stdexcept>
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
// operation on that accessor. Every node below returns a *spin-weighted* node
// from Component<...>(), so the index algebra of upper indices is already
// enforced and there is no second expression system to keep consistent with the
// first. What this layer has to get right is only *which* components it asks
// for.
//
// The concept cannot require Component<...>() itself: the multi-index is a
// pack whose admissible values depend on the tensor's symmetry -- the all-zero
// component of an antisymmetric tensor does not exist -- so there is no single
// instantiation every tensor must provide. It requires the parts that are
// universal, and `Represents` is what a caller consults before asking.
// SlotSet is part of the concept because two operations read it: the tensor
// product, which requires its operands to be drawn from the same alphabet,
// and the contraction, whose sum runs over that alphabet's letters rather
// than over a fixed {-1, 0, +1}.
template <typename T>
concept TensorExpr =
    requires {
      requires std::same_as<std::remove_cv_t<decltype(T::Rank)>,
                            std::ptrdiff_t>;
      requires T::Rank >= 0;
      typename T::GridType;
      requires AngularGrid<typename T::GridType>;
      typename T::SlotSet;
      requires SlotAlphabet<typename T::SlotSet>;
      { T::template Represents<> } -> std::convertible_to<bool>;
    } and
    requires(const T& tensor) {
      { tensor.Grid() } -> std::convertible_to<const typename T::GridType&>;
    }
    // And it is a *spatial* tensor, which has to be said because a spectral one
    // answers every question above: TensorExpansion carries a rank, a grid, a
    // slot alphabet and a Represents, so without this it satisfies the concept
    // and Permute, the tensor product, Materialise and the bundle maps all
    // accept one. They would then compose the wrong Component -- a view over
    // coefficients rather than a spin-weighted node -- and the first symptom is
    // an overload of Tangential resolving to the spatial node for a spectral
    // operand, which is how this was found.
    //
    // The discriminator is the truncation degree, and it is the real difference
    // rather than a convenient one: an expansion is defined up to a degree and
    // a field is not.
    and not requires(const T& tensor) { tensor.MaxDegree(); };

namespace TensorDetails {

using Int = std::ptrdiff_t;

// Whether a node built from an operand passed as T keeps a tensor field alive
// *inside itself*: an rvalue field is moved in, and so is whatever a node
// passed by value already held.
//
// It matters because a component is a view into a field's storage. Over named
// fields a node's components name storage that outlives the node, and taking
// one from a temporary node is as safe as it is natural --
// `Transpose(t).Component<0, 1>()`. Over a field the node owns, the same line
// names storage that dies with the node at the end of the statement. The
// accessors are deleted on rvalues in exactly that case and no other.
template <typename T>
constexpr bool HoldsStorageFn() {
  if constexpr (IsTerminal<T>) {
    return !std::is_lvalue_reference_v<T>;
  } else {
    return std::remove_cvref_t<T>::HoldsStorage;
  }
}

template <typename T>
inline constexpr bool HoldsStorage = HoldsStorageFn<T>();

// Whether an image is a permutation of 0 .. N-1, which is what makes a
// relabelling of slots a tensor operation and not an arbitrary gather.
template <typename I, std::size_t N>
constexpr bool IsPermutation(const std::array<I, N>& image) {
  auto seen = std::array<bool, N>{};
  for (auto to : image) {
    if (to < 0 || static_cast<std::size_t>(to) >= N) return false;
    if (seen[static_cast<std::size_t>(to)]) return false;
    seen[static_cast<std::size_t>(to)] = true;
  }
  return true;
}

// An image as the index type the multi-index uses, so that it can be written
// with plain ints -- `std::array{1, 0}` -- as the documentation writes it.
template <typename I, std::size_t N>
constexpr auto AsIndexArray(const std::array<I, N>& image) {
  auto converted = std::array<Int, N>{};
  for (std::size_t i = 0; i < N; i++) converted[i] = static_cast<Int>(image[i]);
  return converted;
}

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
  return ComponentOf<Indices>(operand,
                              std::make_index_sequence<Indices.size()>{});
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

// The elements of the group a symmetry's generators generate, each a
// permutation with the sign a component picks up under it.
//
// The orbit table walks this group implicitly and never names its elements;
// symmetrisation has to sum over them, so it needs them listed. A breadth-
// first closure under composition, where composing images A then B gives
// C[i] = A[B[i]] -- the convention MultiIndex::Permuted fixes.
template <Int Rank>
constexpr Int FactorialBound() {
  auto n = Int{1};
  for (auto i = Int{2}; i <= Rank; i++) n *= i;
  return n;
}

template <Int Rank, TensorSymmetry<Rank> Symmetry>
constexpr auto GroupElements() {
  constexpr auto Bound = FactorialBound<Rank>();
  auto elements = std::array<SlotPermutation<Rank>, Bound>{};
  auto count = Int{1};
  elements[0] = SlotPermutation<Rank>{SymmetryDetails::Identity<Rank>(), 1};

  const auto generators = Symmetry::Generators();
  for (auto head = Int{0}; head < count; head++) {
    for (const auto& generator : generators) {
      auto image = std::array<Int, Rank>{};
      for (auto i = Int{0}; i < Rank; i++) {
        image[i] = elements[head].image[generator.image[i]];
      }
      const auto sign = elements[head].sign * generator.sign;

      auto seen = false;
      for (auto i = Int{0}; i < count; i++) {
        if (elements[i].image == image) seen = true;
      }
      if (seen || count >= Bound) continue;
      elements[count++] = SlotPermutation<Rank>{image, sign};
    }
  }
  return std::pair(elements, count);
}

}  // namespace TensorDetails

//--------------------------------------------------------------------------//
//                            Permuting the slots                            //
//--------------------------------------------------------------------------//

/// A tensor with its slots relabelled: (Permute T)^{alpha} = T^{pi(alpha)}.
///
/// No arithmetic at all -- the component accessor forwards to the operand with
/// the indices reordered -- which makes this the step that pins the machinery
/// the rest of the layer uses.
///
/// It permutes *tensor slots* and not grid points, so it does not touch the
/// pointwise-and-index-preserving invariant that the aliasing theorem
/// rests on. That is worth stating because a re-indexing view is exactly the
/// thing that would break it, and this is the closest the library comes to one.
template <auto Image, typename Operand>
class PermuteNode {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  /** @brief The operand, with references and cv-qualifiers stripped. */
  using OperandType = std::remove_cvref_t<Operand>;

  /** @brief The tensor rank. */
  static constexpr Int Rank = OperandType::Rank;
  /// The angular grid this is defined on.
  using GridType = typename OperandType::GridType;
  /// The alphabet the slots are drawn from.
  using SlotSet = typename OperandType::SlotSet;

  static_assert(Image.size() == static_cast<std::size_t>(Rank),
                "A slot permutation needs one image per tensor slot");
  static_assert(TensorDetails::IsPermutation(Image),
                "A slot permutation must send the slots 0 .. Rank-1 to "
                "themselves, each once");

  /** @brief Whether this node keeps a tensor field alive inside itself. */
  static constexpr bool HoldsStorage = TensorDetails::HoldsStorage<Operand>;

  /** @brief Wraps an operand whose slots are to be relabelled. */
  explicit PermuteNode(Operand&& operand)
      : operand_{std::forward<Operand>(operand)} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return operand_.Grid(); }

  /// The operand's multi-index for this tensor's alpha: its slot i carries
  /// alpha's slot Image[i], which is MultiIndex::Permuted. So
  ///
  ///     Permute<Image>(T)^{a_0 a_1 ...} = T^{a_Image[0] a_Image[1] ...},
  ///
  /// and with Image = {1, 2, 0} the result R has R^{abc} = T^{bca}. For a
  /// transposition, or any product of disjoint ones, this and the inverse
  /// convention are the same thing, which is how the prose here once came to
  /// state the inverse without any test noticing.
  template <Int... Alphas>
  static constexpr auto Source =
      MultiIndex<Rank, SlotSet>(std::array<Int, Rank>{Alphas...})
          .Permuted(TensorDetails::AsIndexArray(Image))
          .Slots();

  /// The two checks come before the multi-index is formed rather than beside
  /// it, for the reason IsSlotLetter gives: a `and` short-circuits evaluation
  /// but is not a promise about instantiation, and forming Source with a
  /// letter the alphabet does not have is a hard error.
  /** @brief Backs Represents; see there. */
  template <Int... Alphas>
  static constexpr bool RepresentsFn() {
    if constexpr (sizeof...(Alphas) != static_cast<std::size_t>(Rank)) {
      return false;
    } else if constexpr (!AreSlotLetters<SlotSet, Alphas...>()) {
      return false;
    } else {
      return TensorDetails::Represents<Source<Alphas...>, OperandType>();
    }
  }

  /** @brief Whether this holds a component at those slot letters. */
  template <Int... Alphas>
  static constexpr bool Represents = RepresentsFn<Alphas...>();

  /** @brief The component at those slot letters, as a spin-weighted node. */
  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const& {
    return TensorDetails::ComponentOf<Source<Alphas...>>(operand_);
  }

  /// Not offered on a temporary node that owns a tensor field, where the view
  /// returned would name storage gone by the end of the statement: see
  /// TensorDetails::HoldsStorage. Over named fields it is offered as ever.
  template <Int... Alphas>
  requires HoldsStorage
  void Component() const&& = delete;

 private:
  OperandStorage<Operand> operand_;
};

// Relabel a tensor's slots: Permute<Image>(T)^{a_0 a_1 ...} is
// T^{a_Image[0] a_Image[1] ...}. The image is a std::array and may be written
// with plain ints, so Permute<std::array{1, 0}>(t) is the transpose and
// Permute<std::array{1, 2, 0}>(t)^{abc} is t^{bca}. An image that is not a
// permutation of the slots is not an overload.
template <auto Image, typename T>
requires TensorExpr<std::remove_cvref_t<T>> &&
         (Image.size() ==
          static_cast<std::size_t>(std::remove_cvref_t<T>::Rank)) &&
         (TensorDetails::IsPermutation(Image))
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

/// (S tensor T)^{alpha beta} = S^{alpha} T^{beta}, of rank p + q.
///
/// The component is a spin-weighted product, so its upper index is the sum of
/// the two operands' -- which is eq:N applied to the concatenated multi-index,
/// and the theory note, docs/canonical-components.tex, says so in as many
/// words. Nothing here has to arrange that; it is what "upper indices add under
/// pointwise multiplication" means.
template <typename LeftOperand, typename RightOperand>
class TensorProductNode {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  using Left = std::remove_cvref_t<LeftOperand>;  ///< The left operand, bare.
  /// The right operand, bare.
  using Right = std::remove_cvref_t<RightOperand>;

  /** @brief The tensor rank. */
  static constexpr Int Rank = Left::Rank + Right::Rank;
  /// The angular grid this is defined on.
  using GridType = typename Left::GridType;
  /// The alphabet the slots are drawn from.
  using SlotSet = typename Left::SlotSet;

  static_assert(std::same_as<GridType, typename Right::GridType>,
                "A tensor product needs both operands on the same kind of "
                "grid");

  // A tangential tensor lives in a different bundle from a general one, so
  // the product of the two is not a tensor over either alphabet. Crossing
  // bundles is done by Embed at the call site and never implicitly. The
  // constraint is on the free function below, where it is a clean
  // overload-resolution failure; this is
  // the backstop for anyone building the node directly.
  static_assert(std::same_as<SlotSet, typename Right::SlotSet>,
                "A tensor product needs both operands drawn from the same "
                "slot alphabet; embed one of them first");

  /** @brief Whether this node keeps a tensor field alive inside itself. */
  static constexpr bool HoldsStorage =
      TensorDetails::HoldsStorage<LeftOperand> or
      TensorDetails::HoldsStorage<RightOperand>;

  /**
   * @brief Wraps the two operands of a tensor product.
   * @throws std::invalid_argument if they are not on the same grid.
   */
  TensorProductNode(LeftOperand&& left, RightOperand&& right)
      : left_{std::forward<LeftOperand>(left)},
        right_{std::forward<RightOperand>(right)} {
    if (left_.Grid().Identity() != right_.Grid().Identity()) {
      throw std::invalid_argument(
          "A tensor product needs both operands on the same grid");
    }
  }

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return left_.Grid(); }

  /// The multi-index splits: the first p slots name the left operand's
  /// component and the last q the right's.
  template <Int... Alphas>
  static constexpr auto Indices = std::array<Int, Rank>{Alphas...};

  /** @brief The first p slots, naming the left operand's component. */
  template <Int... Alphas>
  static constexpr auto LeftIndices =
      TensorDetails::Head<Left::Rank>(Indices<Alphas...>);

  /** @brief The last q slots, naming the right operand's component. */
  template <Int... Alphas>
  static constexpr auto RightIndices =
      TensorDetails::Tail<Right::Rank>(Indices<Alphas...>);

  /** @brief Backs Represents; see there. */
  template <Int... Alphas>
  static constexpr bool RepresentsFn() {
    if constexpr (sizeof...(Alphas) != static_cast<std::size_t>(Rank)) {
      return false;
    } else if constexpr (!AreSlotLetters<SlotSet, Alphas...>()) {
      return false;
    } else {
      return TensorDetails::Represents<LeftIndices<Alphas...>, Left>() and
             TensorDetails::Represents<RightIndices<Alphas...>, Right>();
    }
  }

  /** @brief Whether this holds a component at those slot letters. */
  template <Int... Alphas>
  static constexpr bool Represents = RepresentsFn<Alphas...>();

  /** @brief The component at those slot letters, as a spin-weighted node. */
  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const& {
    return TensorDetails::ComponentOf<LeftIndices<Alphas...>>(left_) *
           TensorDetails::ComponentOf<RightIndices<Alphas...>>(right_);
  }

  /// Not offered on a temporary node that owns a tensor field, where the view
  /// returned would name storage gone by the end of the statement: see
  /// TensorDetails::HoldsStorage. Over named fields it is offered as ever.
  template <Int... Alphas>
  requires HoldsStorage
  void Component() const&& = delete;

 private:
  OperandStorage<LeftOperand> left_;
  OperandStorage<RightOperand> right_;
};

template <typename L, typename R>
requires TensorExpr<std::remove_cvref_t<L>> &&
         TensorExpr<std::remove_cvref_t<R>> &&
         std::same_as<typename std::remove_cvref_t<L>::SlotSet,
                      typename std::remove_cvref_t<R>::SlotSet>
auto TensorProduct(L&& left, R&& right) {
  return TensorProductNode<L, R>(std::forward<L>(left), std::forward<R>(right));
}

//--------------------------------------------------------------------------//
//                                Contraction                                //
//--------------------------------------------------------------------------//

/// Contract slots J and K against the metric of eq:metric,
///
///   (tr_{JK} T)^{...} = sum_{ab} g_{ab} T^{... a ... b ...}
///                     = sum_a (-1)^a T^{... a ... -a ...},
///
/// which for an ordinary tensor's three values of a is
/// -T^{-+} + T^{00} - T^{+-} in the contracted pair.
///
/// The sum runs over the *alphabet's* letters, not over a fixed {-1, 0, +1}.
/// For a tangential tensor that leaves -T^{-+} - T^{+-}, which is the metric
/// of the sphere induced on the tangent plane -- the same expression with the
/// radial term absent because there is no radial slot to contribute one. The
/// result stays in the same bundle, at rank p - 2, and nothing here has to
/// arrange that either.
///
/// The result lands at the right upper index by construction rather than by
/// arrangement: the contracted pair contributes a + (-a) = 0 whatever a is, so
/// all three terms carry the same upper index and the Equal rule admits
/// their sum. A contraction that paired its slots wrongly would not compile.
template <std::ptrdiff_t J, std::ptrdiff_t K, typename Operand>
class ContractionNode {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  /** @brief The operand, with references and cv-qualifiers stripped. */
  using OperandType = std::remove_cvref_t<Operand>;

  /** @brief The tensor rank. */
  static constexpr Int Rank = OperandType::Rank - 2;
  /// The angular grid this is defined on.
  using GridType = typename OperandType::GridType;
  /// The alphabet the slots are drawn from.
  using SlotSet = typename OperandType::SlotSet;
  using Real = typename GridType::Real;  ///< The precision.

  /** @brief The letters the contracted slots run over. */
  static constexpr auto& Alphabet = SlotSet::Alphabet;
  /** @brief The slot letters contracted over. */
  static constexpr auto Letters = Alphabet.size();

  static_assert(J != K, "A contraction needs two different slots");
  static_assert(J >= 0 && K >= 0 && J < OperandType::Rank &&
                    K < OperandType::Rank,
                "A contracted slot must be one the tensor has");

  /** @brief Wraps an operand, two of whose slots are to be contracted. */
  explicit ContractionNode(Operand&& operand)
      : operand_{std::forward<Operand>(operand)} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return operand_.Grid(); }

  /// The operand's multi-index: the surviving slots in order, with a inserted
  /// at J and -a at K.
  template <Int Alpha, Int... Alphas>
  static constexpr auto Source = TensorDetails::Insert<J, K, OperandType::Rank>(
      std::array<Int, Rank>{Alphas...}, Alpha);

  /** @brief Whether the term of the sum at letter number @p A exists. */
  template <std::size_t A, Int... Alphas>
  static constexpr bool TermExists =
      TensorDetails::Represents<Source<Alphabet[A], Alphas...>, OperandType>();

  /** @brief Whether this node keeps a tensor field alive inside itself. */
  static constexpr bool HoldsStorage = TensorDetails::HoldsStorage<Operand>;

  /** @brief Backs Represents; see there. */
  template <Int... Alphas>
  static constexpr bool RepresentsFn() {
    if constexpr (sizeof...(Alphas) != static_cast<std::size_t>(Rank)) {
      return false;
    } else if constexpr (!AreSlotLetters<SlotSet, Alphas...>()) {
      return false;
    } else {
      // A term the operand does not represent is identically zero -- the
      // diagonal of an antisymmetric tensor, a radial component of an
      // embedded tangential one -- so it is a term to leave out of the sum,
      // and the sum exists if any term does. It was once required that every
      // term exist, which made A.v for an antisymmetric A wholly
      // unrepresented, and Materialise wrote nothing for it.
      return [&]<std::size_t... A>(std::index_sequence<A...>) {
        return (TermExists<A, Alphas...> or ...);
      }(std::make_index_sequence<Letters>{});
    }
  }

  /** @brief Whether this holds a component at those slot letters. */
  template <Int... Alphas>
  static constexpr bool Represents = RepresentsFn<Alphas...>();

  /** @brief The component at those slot letters, as a spin-weighted node. */
  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const& {
    return SumFrom<0, Alphas...>();
  }

  /// Not offered on a temporary node that owns a tensor field, where the view
  /// returned would name storage gone by the end of the statement: see
  /// TensorDetails::HoldsStorage. Over named fields it is offered as ever.
  template <Int... Alphas>
  requires HoldsStorage
  void Component() const&& = delete;

 private:
  OperandStorage<Operand> operand_;

  // One term: the metric's (-1)^a times the operand's component. Written as a
  // scalar multiplication rather than as unary minus so that one expression
  // covers both alphabets; the factor is exactly +-1, so no arithmetic is
  // added.
  template <std::size_t A, Int... Alphas>
  auto Term() const {
    return MinusOneToPower<Real>(Alphabet[A]) *
           TensorDetails::ComponentOf<Source<Alphabet[A], Alphas...>>(operand_);
  }

  template <std::size_t A, Int... Alphas>
  static constexpr bool AnyTermExistsFrom() {
    const auto exists = [&]<std::size_t... All>(std::index_sequence<All...>) {
      return std::array<bool, Letters>{TermExists<All, Alphas...>...};
    }(std::make_index_sequence<Letters>{});
    for (auto a = A; a < Letters; a++) {
      if (exists[a]) return true;
    }
    return false;
  }

  // The terms that exist, from letter number A on, summed from the right --
  // t_A + (t_B + t_C) -- which is the association a fold over the alphabet
  // gives, so where every term exists the value is what it always was.
  template <std::size_t A, Int... Alphas>
  auto SumFrom() const {
    if constexpr (!TermExists<A, Alphas...>) {
      return SumFrom<A + 1, Alphas...>();
    } else if constexpr (!AnyTermExistsFrom<A + 1, Alphas...>()) {
      return Term<A, Alphas...>();
    } else {
      return Term<A, Alphas...>() + SumFrom<A + 1, Alphas...>();
    }
  }
};

template <std::ptrdiff_t J, std::ptrdiff_t K, typename T>
requires TensorExpr<std::remove_cvref_t<T>>
auto Contract(T&& tensor) {
  return ContractionNode<J, K, T>(std::forward<T>(tensor));
}

// The trace of a rank-2 tensor, which is the contraction with a name.
//
// It returns a spin-weighted field and not a tensor, and that field is an
// expression over views into the argument's storage. So the argument has to
// outlive it, and a temporary tensor field -- or an expression that has taken
// ownership of one -- does not: the contraction would own the field, hand out
// views into it, and die on the way out of this function. That case is not an
// overload. Name the tensor first.
template <typename T>
requires TensorExpr<std::remove_cvref_t<T>> &&
         (std::remove_cvref_t<T>::Rank == 2) &&
         (!TensorDetails::HoldsStorage<T>)
auto Trace(T&& tensor) {
  const auto contracted = Contract<0, 1>(std::forward<T>(tensor));
  return contracted.template Component<>();
}

//--------------------------------------------------------------------------//
//                              Symmetrisation                               //
//--------------------------------------------------------------------------//

/// The projection onto a symmetry:
///
///   Sym(T)^{alpha} = (1/|G|) sum_{pi in G} sign(pi) T^{pi(alpha)}.
///
/// Every term carries the same upper index, since a permutation preserves the
/// slot sum, so the sum is admissible and the result is again a tensor. What it
/// is *not* is a tensor whose type records the symmetry: the value has the
/// property, and only Materialise can be asked to store it that way.
template <typename Symmetry, typename Operand>
class SymmetriseNode {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  /** @brief The operand, with references and cv-qualifiers stripped. */
  using OperandType = std::remove_cvref_t<Operand>;

  /** @brief The tensor rank. */
  static constexpr Int Rank = OperandType::Rank;
  /// The angular grid this is defined on.
  using GridType = typename OperandType::GridType;
  /// The alphabet the slots are drawn from.
  using SlotSet = typename OperandType::SlotSet;
  using Real = typename GridType::Real;  ///< The precision.

  /** @brief The elements of the symmetry group, and how many there are. */
  static constexpr auto Group = TensorDetails::GroupElements<Rank, Symmetry>();
  /** @brief The order of the symmetry group, which is what the sum divides by.
   */
  static constexpr Int GroupSize = Group.second;

  /** @brief Wraps an operand to be averaged over the symmetry group. */
  explicit SymmetriseNode(Operand&& operand)
      : operand_{std::forward<Operand>(operand)} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return operand_.Grid(); }

  /** @brief The operand's slot letters under group element @p Element. */
  template <std::size_t Element, Int... Alphas>
  static constexpr auto Source =
      MultiIndex<Rank, SlotSet>(std::array<Int, Rank>{Alphas...})
          .Permuted(Group.first[Element].image)
          .Slots();

  /** @brief Whether the term of the sum at group element @p Element exists. */
  template <std::size_t Element, Int... Alphas>
  static constexpr bool TermExists =
      TensorDetails::Represents<Source<Element, Alphas...>, OperandType>();

  /** @brief Whether this node keeps a tensor field alive inside itself. */
  static constexpr bool HoldsStorage = TensorDetails::HoldsStorage<Operand>;

  /** @brief Backs Represents; see there. */
  template <Int... Alphas>
  static constexpr bool RepresentsFn() {
    if constexpr (sizeof...(Alphas) != static_cast<std::size_t>(Rank)) {
      return false;
    } else if constexpr (!AreSlotLetters<SlotSet, Alphas...>()) {
      return false;
    } else {
      // Any and not all, as for a contraction and for the same reason: a
      // term the operand does not represent is zero, and is left out.
      return [&]<std::size_t... E>(std::index_sequence<E...>) {
        return (TermExists<E, Alphas...> or ...);
      }(std::make_index_sequence<GroupSize>{});
    }
  }

  /** @brief Whether this holds a component at those slot letters. */
  template <Int... Alphas>
  static constexpr bool Represents = RepresentsFn<Alphas...>();

  /** @brief The component at those slot letters, as a spin-weighted node. */
  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const& {
    // Divided by the order of the group whatever is left out: a term that
    // does not exist is a zero in the average, not a smaller average.
    return SumFrom<0, Alphas...>() / static_cast<Real>(GroupSize);
  }

  /// Not offered on a temporary node that owns a tensor field, where the view
  /// returned would name storage gone by the end of the statement: see
  /// TensorDetails::HoldsStorage. Over named fields it is offered as ever.
  template <Int... Alphas>
  requires HoldsStorage
  void Component() const&& = delete;

 private:
  OperandStorage<Operand> operand_;

  // One term: the operand's component under one element of the group, scaled
  // by that element's sign. The terms need not have the same type, since a
  // component may come back as a view from one element and as an expression
  // from another, and the binary node does not care.
  template <std::size_t E, Int... Alphas>
  auto Term() const {
    return static_cast<Real>(Group.first[E].sign) *
           TensorDetails::ComponentOf<Source<E, Alphas...>>(operand_);
  }

  template <std::size_t E, Int... Alphas>
  static constexpr bool AnyTermExistsFrom() {
    const auto exists = [&]<std::size_t... All>(std::index_sequence<All...>) {
      return std::array<bool, static_cast<std::size_t>(GroupSize)>{
          TermExists<All, Alphas...>...};
    }(std::make_index_sequence<GroupSize>{});
    for (auto e = E; e < static_cast<std::size_t>(GroupSize); e++) {
      if (exists[e]) return true;
    }
    return false;
  }

  // The terms that exist, summed from the right as a fold over the group
  // would sum them, so that where every term exists nothing has changed.
  template <std::size_t E, Int... Alphas>
  auto SumFrom() const {
    if constexpr (!TermExists<E, Alphas...>) {
      return SumFrom<E + 1, Alphas...>();
    } else if constexpr (!AnyTermExistsFrom<E + 1, Alphas...>()) {
      return Term<E, Alphas...>();
    } else {
      return Term<E, Alphas...>() + SumFrom<E + 1, Alphas...>();
    }
  }
};

// Project onto a symmetry. Symmetrise<Symmetric<2>>(t) is the symmetric part,
// Symmetrise<Antisymmetric<2>>(t) the antisymmetric one.
template <typename Symmetry, typename T>
requires TensorExpr<std::remove_cvref_t<T>>
auto Symmetrise(T&& tensor) {
  return SymmetriseNode<Symmetry, T>(std::forward<T>(tensor));
}

//--------------------------------------------------------------------------//
//                              Materialisation                              //
//--------------------------------------------------------------------------//

namespace TensorDetails {

// Copy one component of an expression into the same component of a field.
//
// Element by element rather than through EvaluateInto, because the target may
// be strided: a point-major field's component is not contiguous, and
// EvaluateInto writes a contiguous span. The contiguous case could take the
// faster path; it is not worth two code paths until something measures it.
template <auto Indices, typename Field, typename Expr, std::size_t... I>
void AssignComponent(Field& field, const Expr& expr,
                     std::index_sequence<I...>) {
  auto target = field.template Component<Indices[I]...>();
  const auto source = expr.template Component<Indices[I]...>();
  using Target = std::remove_cvref_t<decltype(target)>;
  using Source = std::remove_cvref_t<decltype(source)>;

  constexpr auto flat =
      MultiIndex<Field::Rank, typename Field::SlotSet>(Indices).Flat();
  constexpr auto constraint = Field::Orbits.constraint[flat];
  constexpr bool narrowing =
      std::same_as<typename Target::Value, RealValued> and
      std::same_as<typename Source::Value, ComplexValued>;

  const auto& grid = field.Grid();
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      if constexpr (narrowing) {
        // The target is a component the reality condition pins to a single
        // real number and the expression does not know that. Taking the part
        // that survives is the caller's assertion that the expression really
        // is a real tensor -- the same kind of assertion asking for a
        // symmetry is, and unprovable here for the same reason.
        if constexpr (constraint == ComponentConstraint::Imaginary) {
          target[iTheta, iPhi] = std::imag(source[iTheta, iPhi]);
        } else {
          target[iTheta, iPhi] = std::real(source[iTheta, iPhi]);
        }
      } else {
        target[iTheta, iPhi] = source[iTheta, iPhi];
      }
    }
  }
}

template <typename Field, typename Expr, std::size_t Slot>
void AssignSlot(Field& field, const Expr& expr) {
  constexpr auto flat = Field::ComponentLayout.flatOfSlot[Slot];
  constexpr auto indices =
      MultiIndex<Field::Rank, typename Field::SlotSet>::FromFlat(flat).Slots();
  if constexpr (Represents<indices, std::remove_cvref_t<Expr>>()) {
    AssignComponent<indices>(field, expr,
                             std::make_index_sequence<Field::Rank>{});
  }
}

}  // namespace TensorDetails

// Evaluate a tensor expression into a field, which is where a lazy tensor
// stops being lazy.
//
// The symmetry and the reality are the caller's to state, and default to none
// and complex. The product of two symmetric tensors is not symmetric, and
// inferring either property from an expression tree is a research problem
// rather than a design -- so asking for one here is an assertion about the
// value, honoured by storing only the components it keeps. Asking for
// RealTensor stores half as much and derives the rest, and where a component
// is pinned to one real number, that number is taken from the expression.
template <typename Symmetry = void, TensorReality Reality = ComplexTensor,
          typename Expr>
requires TensorExpr<std::remove_cvref_t<Expr>>
auto Materialise(const Expr& expr) {
  using E = std::remove_cvref_t<Expr>;
  using Chosen = std::conditional_t<std::same_as<Symmetry, void>,
                                    NoSymmetry<E::Rank>, Symmetry>;
  // The alphabet is the expression's, not a choice: materialising cannot move
  // a tensor between bundles. Layout is named explicitly only because SlotSet
  // sits behind it in the parameter list.
  using Field = TensorField<E::Rank, Chosen, Reality, typename E::GridType,
                            ComponentMajor, typename E::SlotSet>;

  auto field = Field(expr.Grid());
  [&]<std::size_t... Slot>(std::index_sequence<Slot...>) {
    (TensorDetails::AssignSlot<Field, E, Slot>(field, expr), ...);
  }(std::make_index_sequence<static_cast<std::size_t>(
        Field::StoredComponents)>{});
  return field;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_TENSOR_EXPR_GUARD_H
