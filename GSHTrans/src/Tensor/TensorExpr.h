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
// operation on that accessor. Every node below returns a *spin-weighted* node from
// Component<...>(), so the index algebra of upper indices is already enforced
// and there is no second expression system to keep consistent with the first.
// What this layer has to get right is only *which* components it asks for.
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
concept TensorExpr = requires {
  requires std::same_as<std::remove_cv_t<decltype(T::Rank)>, std::ptrdiff_t>;
  requires T::Rank >= 0;
  typename T::GridType;
  requires AngularGrid<typename T::GridType>;
  typename T::SlotSet;
  requires SlotAlphabet<typename T::SlotSet>;
  { T::template Represents<> } -> std::convertible_to<bool>;
} and requires(const T& tensor) {
  { tensor.Grid() } -> std::convertible_to<const typename T::GridType&>;
}
// And it is a *spatial* tensor, which has to be said because a spectral one
// answers every question above: TensorExpansion carries a rank, a grid, a slot
// alphabet and a Represents, so without this it satisfies the concept and
// Permute, the tensor product, Materialise and the bundle maps all accept one.
// They would then compose the wrong Component -- a view over coefficients
// rather than a spin-weighted node -- and the first symptom is an overload of
// Tangential resolving to the spatial node for a spectral operand, which is
// how this was found.
//
// The discriminator is the truncation degree, and it is the real difference
// rather than a convenient one: an expansion is defined up to a degree and a
// field is not.
and not requires(const T& tensor) { tensor.MaxDegree(); };

template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid, TensorLayout Layout,
          SlotAlphabet Slots>
struct IsTerminalTrait<
    TensorField<Rank, Symmetry, Reality, Grid, Layout, Slots>>
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

  /** @brief Wraps an operand whose slots are to be relabelled. */
  explicit PermuteNode(Operand&& operand)
      : _operand{std::forward<Operand>(operand)} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _operand.Grid(); }

  /// Slot i of this tensor is slot Image[i] of the operand, which is the same
  /// convention MultiIndex::Permuted uses.
  template <Int... Alphas>
  static constexpr auto Source =
      MultiIndex<Rank, SlotSet>(std::array<Int, Rank>{Alphas...})
          .Permuted(Image)
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

/// (S tensor T)^{alpha beta} = S^{alpha} T^{beta}, of rank p + q.
///
/// The component is a spin-weighted product, so its upper index is the sum of the
/// two operands' -- which is eq:N applied to the concatenated multi-index, and
/// the theory note, docs/canonical-components.tex, says so in as many words.
/// Nothing here has to arrange that;
/// it is what "upper indices add under pointwise multiplication" means.
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

  /**
   * @brief Wraps the two operands of a tensor product.
   * @throws std::invalid_argument if they are not on the same grid.
   */
  TensorProductNode(LeftOperand&& left, RightOperand&& right)
      : _left{std::forward<LeftOperand>(left)},
        _right{std::forward<RightOperand>(right)} {
    if (_left.Grid().Identity() != _right.Grid().Identity()) {
      throw std::invalid_argument(
          "A tensor product needs both operands on the same grid");
    }
  }

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _left.Grid(); }

  /// The multi-index splits: the first p slots name the left operand's
  /// component and the last q the right's.
  template <Int... Alphas>
  static constexpr auto Indices = std::array<Int, Rank>{Alphas...};

  /** @brief The first p slots, naming the left operand's component. */
  template <Int... Alphas>
  static constexpr auto LeftIndices = TensorDetails::Head<Left::Rank>(
      Indices<Alphas...>);

  /** @brief The last q slots, naming the right operand's component. */
  template <Int... Alphas>
  static constexpr auto RightIndices = TensorDetails::Tail<Right::Rank>(
      Indices<Alphas...>);

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
  auto Component() const {
    return TensorDetails::ComponentOf<LeftIndices<Alphas...>>(_left) *
           TensorDetails::ComponentOf<RightIndices<Alphas...>>(_right);
  }

 private:
  OperandStorage<LeftOperand> _left;
  OperandStorage<RightOperand> _right;
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
      : _operand{std::forward<Operand>(operand)} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _operand.Grid(); }

  /// The operand's multi-index: the surviving slots in order, with a inserted
  /// at J and -a at K.
  template <Int Alpha, Int... Alphas>
  static constexpr auto Source = TensorDetails::Insert<J, K, OperandType::Rank>(
      std::array<Int, Rank>{Alphas...}, Alpha);

  /** @brief Backs Represents; see there. */
  template <Int... Alphas>
  static constexpr bool RepresentsFn() {
    if constexpr (sizeof...(Alphas) != static_cast<std::size_t>(Rank)) {
      return false;
    } else if constexpr (!AreSlotLetters<SlotSet, Alphas...>()) {
      return false;
    } else {
      // Every term of the sum has to exist, since they are added. A tensor
      // with a vanishing component in the middle of a contraction is not
      // something this layer tries to be clever about.
      return [&]<std::size_t... A>(std::index_sequence<A...>) {
        return (TensorDetails::Represents<Source<Alphabet[A], Alphas...>,
                                          OperandType>() and
                ...);
      }(std::make_index_sequence<Letters>{});
    }
  }

  /** @brief Whether this holds a component at those slot letters. */
  template <Int... Alphas>
  static constexpr bool Represents = RepresentsFn<Alphas...>();

  /** @brief The component at those slot letters, as a spin-weighted node. */
  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const {
    return Sum<Alphas...>(std::make_index_sequence<Letters>{});
  }

 private:
  OperandStorage<Operand> _operand;

  // The metric's (-1)^a, folded over the alphabet. Written as a scalar
  // multiplication rather than as unary minus so that one expression covers
  // both alphabets; the factor is exactly +-1, so no arithmetic is added.
  template <Int... Alphas, std::size_t... A>
  auto Sum(std::index_sequence<A...>) const {
    return ((MinusOneToPower<Real>(Alphabet[A]) *
             TensorDetails::ComponentOf<Source<Alphabet[A], Alphas...>>(
                 _operand)) +
            ...);
  }
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
  /** @brief The order of the symmetry group, which is what the sum divides by. */
  static constexpr Int GroupSize = Group.second;

  /** @brief Wraps an operand to be averaged over the symmetry group. */
  explicit SymmetriseNode(Operand&& operand)
      : _operand{std::forward<Operand>(operand)} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _operand.Grid(); }

  /** @brief The operand's slot letters under group element @p Element. */
  template <std::size_t Element, Int... Alphas>
  static constexpr auto Source =
      MultiIndex<Rank, SlotSet>(std::array<Int, Rank>{Alphas...})
          .Permuted(Group.first[Element].image)
          .Slots();

  /** @brief Backs Represents; see there. */
  template <Int... Alphas>
  static constexpr bool RepresentsFn() {
    if constexpr (sizeof...(Alphas) != static_cast<std::size_t>(Rank)) {
      return false;
    } else if constexpr (!AreSlotLetters<SlotSet, Alphas...>()) {
      return false;
    } else {
      return [&]<std::size_t... E>(std::index_sequence<E...>) {
        return (TensorDetails::Represents<Source<E, Alphas...>,
                                          OperandType>() and
                ...);
      }(std::make_index_sequence<GroupSize>{});
    }
  }

  /** @brief Whether this holds a component at those slot letters. */
  template <Int... Alphas>
  static constexpr bool Represents = RepresentsFn<Alphas...>();

  /** @brief The component at those slot letters, as a spin-weighted node. */
  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const {
    return Sum<Alphas...>(std::make_index_sequence<GroupSize>{}) /
           static_cast<Real>(GroupSize);
  }

 private:
  OperandStorage<Operand> _operand;

  // A fold over the group. Each term is the operand's component under one
  // element, scaled by that element's sign; the terms need not have the same
  // type, since a component may come back as a view from one element and as an
  // expression from another, and the binary node does not care.
  template <Int... Alphas, std::size_t... E>
  auto Sum(std::index_sequence<E...>) const {
    return ((static_cast<Real>(Group.first[E].sign) *
             TensorDetails::ComponentOf<Source<E, Alphas...>>(_operand)) +
            ...);
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
void AssignComponent(Field& field, const Expr& expr, std::index_sequence<I...>) {
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
  using Chosen =
      std::conditional_t<std::same_as<Symmetry, void>, NoSymmetry<E::Rank>,
                         Symmetry>;
  // The alphabet is the expression's, not a choice: materialising cannot move
  // a tensor between bundles. Layout is named explicitly only because SlotSet
  // sits behind it in the parameter list.
  using Field = TensorField<E::Rank, Chosen, Reality, typename E::GridType,
                            ComponentMajor, typename E::SlotSet>;

  auto field = Field(expr.Grid());
  [&]<std::size_t... Slot>(std::index_sequence<Slot...>) {
    (TensorDetails::AssignSlot<Field, E, Slot>(field, expr), ...);
  }(std::make_index_sequence<static_cast<std::size_t>(Field::StoredComponents)>{});
  return field;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_TENSOR_EXPR_GUARD_H
