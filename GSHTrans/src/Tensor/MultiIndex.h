#ifndef GSH_TRANS_MULTI_INDEX_GUARD_H
#define GSH_TRANS_MULTI_INDEX_GUARD_H

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <stdexcept>
#include <type_traits>

#include "../Concepts.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                              Which slots exist                            //
//--------------------------------------------------------------------------//

/// The alphabet a tensor's slots are drawn from.
///
/// An ordinary canonical tensor has slots in {-1, 0, +1}. A *tangential* one --
/// surface strain and stress, the metric of the sphere, the spin-2 fields of
/// CMB and geodesy -- has no radial slot at all, so its slots run over
/// {-1, +1}. That is the T^{Omega Omega} piece of the decomposition Dahlen &
/// Tromp open Appendix C with, and it is what most geophysical surface objects
/// are.
///
/// A policy rather than a second class, for the reason every other such choice
/// here is one: it keeps a single MultiIndex, so that every algorithm written
/// against it stays written once. Orbits.h walks the group it is given over the
/// components it is told about, the symmetry policies permute slots by position
/// and never look at their contents, and storage groups by the slot sum -- none
/// of which knows or needs to know how many values a slot can take.
///
/// The alphabet is deliberately uniform across slots rather than per slot. A
/// mixed object such as D&T's T^{r Omega} would break an invariant the symmetry
/// machinery relies on -- a permutation may only exchange slots drawn from the
/// same alphabet -- and widening to a per-slot list later costs only this
/// header.
///
/// Both alphabets are closed under negation, which is what the reality
/// condition needs: it maps -1 to +1 and 0 to itself.
struct AllSlots {
  /** @brief The three canonical letters. */
  static constexpr std::array<std::ptrdiff_t, 3> Alphabet{-1, 0, 1};
};

/// The slots of a tangential tensor, which has no radial one: its components
/// number 2^Rank rather than 3^Rank, and the derivative closed on it is the
/// intrinsic one.
struct TangentialSlots {
  /** @brief The two canonical letters, the radial one being absent. */
  static constexpr std::array<std::ptrdiff_t, 2> Alphabet{-1, 1};
};

template <typename S>
concept SlotAlphabet = requires {
  { S::Alphabet } -> std::ranges::sized_range;
  requires std::same_as<std::ranges::range_value_t<decltype(S::Alphabet)>,
                        std::ptrdiff_t>;
  requires S::Alphabet.size() > 0;
};

// Whether a value is one of the alphabet's letters, and whether every letter
// of a candidate multi-index is.
//
// This has to be askable *before* a multi-index is formed. The constructor
// throws on a letter it does not recognise, and a throw in a constant
// expression is a hard error at the point of use rather than something a
// requires-expression can see: over TangentialSlots, `requires { Flat<0, 1>; }`
// is **true** and the use then fails to compile. So every accessor taking a
// component's letters as template arguments asks this first, in the same
// `if constexpr` shape the pack-size check uses. It is written down here
// rather than at one call site because every one of them needs it.
template <SlotAlphabet Slots>
constexpr bool IsSlotLetter(std::ptrdiff_t alpha) {
  for (auto letter : Slots::Alphabet) {
    if (letter == alpha) return true;
  }
  return false;
}

template <SlotAlphabet Slots, std::ptrdiff_t... Alphas>
constexpr bool AreSlotLetters() {
  return (IsSlotLetter<Slots>(Alphas) and ...);
}

// Whether a multi-index has a radial slot, in the two spellings the bundle
// maps meet it in: as a pack of letters and as a compile-time array.
//
// The radial letter is zero, so this is what tells a component of the general
// bundle from one that is also a component of the tangential bundle -- the
// question both Embed and Tangential are built on, spatially and spectrally.
// It lives here rather than in either of them because it is a fact about the
// alphabet, and because two copies of it drifting apart is exactly the kind
// of thing that would not show up as a compile error.
/** @brief Whether a multi-index given as a pack of letters has a radial
 * slot. */
template <std::ptrdiff_t... Alphas>
constexpr bool HasRadialSlot() {
  return ((Alphas == 0) or ...);
}

/** @brief Whether a multi-index given as a compile-time array has a radial
 * slot. */
template <auto Indices>
constexpr bool HasRadialSlotIn() {
  for (auto alpha : Indices) {
    if (alpha == 0) return true;
  }
  return false;
}

//--------------------------------------------------------------------------//
//                              The multi-index                              //
//--------------------------------------------------------------------------//

/// The label of one canonical component of a rank-p tensor: p slots, each
/// carrying alpha drawn from the alphabet above; see section 2 of the theory
/// note, docs/canonical-components.tex.
///
/// The distinction this type exists to keep is that the multi-index is *not*
/// the upper index. The upper index is the signed sum of the slots (eq:N), and
/// for rank >= 2 several multi-indices share one: a rank-2 tensor has three
/// components at N = 0, namely (-+), (00) and (+-). A collection labelled only
/// by N therefore does not determine a tensor, which is why the unit below is
/// called SpinField and why components are addressed by multi-index here.
template <std::ptrdiff_t _Rank, SlotAlphabet _Slots = AllSlots>
class MultiIndex {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  // Named SlotSet rather than Slots, which is taken by the accessor below.
  using SlotSet = _Slots;  ///< The alphabet the slots are drawn from.

  /** @brief The tensor rank. */
  static constexpr Int Rank = _Rank;
  static_assert(Rank >= 0, "A tensor rank cannot be negative");

  /// How many values one slot can take: the radix everything below counts in.
  static constexpr Int Base = static_cast<Int>(SlotSet::Alphabet.size());

  /// Base^Rank, the number of components -- 3^Rank for an ordinary tensor and
  /// 2^Rank for a tangential one. Formed by repeated multiplication rather than
  /// by pow, which is neither constexpr nor exact.
  static constexpr Int Size = [] {
    auto size = Int{1};
    for (auto i = Int{0}; i < Rank; i++) size *= Base;
    return size;
  }();

  constexpr MultiIndex() : _slots{} {
    // Not the all-zero index in general: zero is not in every alphabet. The
    // first letter is, and is the natural default for the same reason it is
    // the first flat index.
    for (auto& alpha : _slots) alpha = SlotSet::Alphabet[0];
  }

  /**
   * @brief The index with those slot letters.
   * @param slots One letter per slot, each drawn from the alphabet.
   * @throws std::invalid_argument if a letter is not in the alphabet.
   */
  constexpr explicit MultiIndex(std::array<Int, Rank> slots) : _slots{slots} {
    for (auto alpha : _slots) {
      if (DigitOf(alpha) < 0) {
        throw std::invalid_argument(
            "A canonical index must be drawn from the slot alphabet");
      }
    }
  }

  /// Slot 0 is the most significant digit, so that flat indices run in
  /// lexicographic order of (alpha_1 ... alpha_p) starting from the alphabet's
  /// first letter. The choice is arbitrary but has to be pinned: it is the
  /// order components are laid out in and the order any table below is
  /// generated in.
  static constexpr MultiIndex FromFlat(Int flat) {
    if (flat < 0 || flat >= Size) {
      throw std::invalid_argument("Flat component index out of range");
    }
    auto slots = std::array<Int, Rank>{};
    for (auto i = Rank - 1; i >= 0; i--) {
      slots[i] = SlotSet::Alphabet[flat % Base];
      flat /= Base;
    }
    return MultiIndex(slots);
  }

  /// The flat component index, in the lexicographic order FromFlat describes.
  constexpr Int Flat() const {
    auto flat = Int{0};
    for (auto alpha : _slots) flat = Base * flat + DigitOf(alpha);
    return flat;
  }

  /// The upper index, eq:N. This is the spin weight of the component's field
  /// and therefore a compile-time quantity wherever it is used.
  constexpr Int UpperIndex() const {
    auto n = Int{0};
    for (auto alpha : _slots) n += alpha;
    return n;
  }

  /** @brief The letter at slot @p slot. */
  constexpr Int operator[](Int slot) const { return _slots[slot]; }
  /** @brief Every letter, one per slot. */
  constexpr const std::array<Int, Rank>& Slots() const { return _slots; }

  /// The involution of the reality condition (eq:reality). Its only fixed point
  /// is the all-zero index -- which for a tangential tensor of rank >= 1 does
  /// not exist, so negation is then free of fixed points and every orbit has
  /// size two. That is what removes the pinned components, and with them the
  /// second buffer, from a tangential real tensor with no permutation
  /// symmetry.
  constexpr MultiIndex Negated() const {
    auto slots = _slots;
    for (auto& alpha : slots) alpha = -alpha;
    return MultiIndex(slots);
  }

  /// Slot i of the result is slot image[i] of this one.
  constexpr MultiIndex Permuted(const std::array<Int, Rank>& image) const {
    auto slots = std::array<Int, Rank>{};
    for (auto i = Int{0}; i < Rank; i++) slots[i] = _slots[image[i]];
    return MultiIndex(slots);
  }

  /** @brief Compares slot by slot. */
  constexpr bool operator==(const MultiIndex&) const = default;

 private:
  // The position of a letter in the alphabet, or -1 if it is not one. This is
  // the only place that knows the alphabet is a list rather than a range, and
  // a linear scan over two or three entries is not worth improving.
  static constexpr Int DigitOf(Int alpha) {
    for (auto i = Int{0}; i < Base; i++) {
      if (SlotSet::Alphabet[i] == alpha) return i;
    }
    return -1;
  }

  std::array<Int, Rank> _slots;
};

// How many components of a rank-p tensor carry upper index N: the trinomial
// coefficient, which for p = 2 gives 1, 2, 3, 2, 1 (theory note section 2).
// For a tangential tensor it is the binomial coefficient instead, and it
// vanishes unless N has the same parity as the rank. Counted rather than
// derived, since the enumeration is cheap and either closed form is one more
// thing to get wrong.
template <std::ptrdiff_t Rank, SlotAlphabet Slots = AllSlots>
constexpr auto ComponentsAtUpperIndex(std::ptrdiff_t n) {
  using Index = MultiIndex<Rank, Slots>;
  auto count = std::ptrdiff_t{0};
  for (auto flat = std::ptrdiff_t{0}; flat < Index::Size; flat++) {
    if (Index::FromFlat(flat).UpperIndex() == n) count++;
  }
  return count;
}

//--------------------------------------------------------------------------//
//                            Permutation symmetry                           //
//--------------------------------------------------------------------------//

/// A permutation of a tensor's slots, together with the sign the component
/// picks up under it: +1 for a symmetric generator, -1 for an antisymmetric
/// one.
template <std::ptrdiff_t Rank>
struct SlotPermutation {
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.

  std::array<Int, Rank> image;  ///< Slot i of the result is slot image[i].
  Int sign;  ///< The sign the component picks up, +1 or -1.

  /** @brief Compares componentwise. */
  constexpr bool operator==(const SlotPermutation&) const = default;
};

namespace SymmetryDetails {

using Int = std::ptrdiff_t;

template <Int Rank>
constexpr auto Identity() {
  auto image = std::array<Int, Rank>{};
  for (auto i = Int{0}; i < Rank; i++) image[i] = i;
  return image;
}

template <Int Rank>
constexpr auto Transposition(Int a, Int b, Int sign) {
  auto image = Identity<Rank>();
  image[a] = b;
  image[b] = a;
  return SlotPermutation<Rank>{image, sign};
}

}  // namespace SymmetryDetails

/// No relation between components: every one of the 3^Rank is independent.
template <std::ptrdiff_t Rank>
struct NoSymmetry {
  /** @brief Generators of the symmetry group. */
  static constexpr auto Generators() {
    return std::array<SlotPermutation<Rank>, 0>{};
  }
};

/// Symmetric or antisymmetric in every pair of slots. The adjacent
/// transpositions generate the whole symmetric group, so two generators
/// suffice at any rank and the orbit search finds the rest.
template <std::ptrdiff_t Rank, std::ptrdiff_t Sign>
struct FullPermutationSymmetry {
  static_assert(Sign == 1 || Sign == -1);

  /** @brief Generators of the symmetry group. */
  static constexpr auto Generators() {
    constexpr auto count = Rank > 1 ? Rank - 1 : 0;
    auto generators = std::array<SlotPermutation<Rank>, count>{};
    for (auto i = std::ptrdiff_t{0}; i < count; i++) {
      generators[i] = SymmetryDetails::Transposition<Rank>(i, i + 1, Sign);
    }
    return generators;
  }
};

template <std::ptrdiff_t Rank>
using Symmetric = FullPermutationSymmetry<Rank, 1>;

template <std::ptrdiff_t Rank>
using Antisymmetric = FullPermutationSymmetry<Rank, -1>;

/// Anything else, given by its generators. The case this exists for is the
/// elastic tensor, c_{ijkl} = c_{jikl} = c_{ijlk} = c_{klij}, which is three
/// generators and is why rank 4 is exposed at all.
template <std::ptrdiff_t Rank, SlotPermutation<Rank>... Gs>
struct GeneratedBy {
  /** @brief Generators of the symmetry group. */
  static constexpr auto Generators() {
    return std::array<SlotPermutation<Rank>, sizeof...(Gs)>{Gs...};
  }
};

/// The symmetry of an elastic tensor: symmetric within each pair of slots, and
/// under exchange of the pairs.
struct ElasticSymmetry {
  /** @brief The tensor rank this symmetry is defined for. */
  static constexpr std::ptrdiff_t Rank = 4;

  /** @brief Generators of the symmetry group. */
  static constexpr auto Generators() {
    return std::array<SlotPermutation<4>, 3>{
        SymmetryDetails::Transposition<4>(0, 1, 1),
        SymmetryDetails::Transposition<4>(2, 3, 1),
        SlotPermutation<4>{std::array<std::ptrdiff_t, 4>{2, 3, 0, 1}, 1}};
  }
};

template <typename S, std::ptrdiff_t Rank>
concept TensorSymmetry = requires {
  { S::Generators() } -> std::ranges::sized_range;
  requires std::same_as<
      std::ranges::range_value_t<decltype(S::Generators())>,
      SlotPermutation<Rank>>;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_MULTI_INDEX_GUARD_H
