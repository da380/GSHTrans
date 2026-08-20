#ifndef GSH_TRANS_MULTI_INDEX_GUARD_H
#define GSH_TRANS_MULTI_INDEX_GUARD_H

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <stdexcept>

#include "../Concepts.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                              The multi-index                              //
//--------------------------------------------------------------------------//

// The label of one canonical component of a rank-p tensor: p slots, each
// carrying alpha in {-1, 0, +1} (theory note section 2).
//
// The distinction this type exists to keep is that the multi-index is *not*
// the upper index. The upper index is the signed sum of the slots (eq:N), and
// for rank >= 2 several multi-indices share one: a rank-2 tensor has three
// components at N = 0, namely (-+), (00) and (+-). A collection labelled only
// by N therefore does not determine a tensor, which is why phase 1's unit is
// called SpinField and why components are addressed by multi-index here.
template <std::ptrdiff_t _Rank>
class MultiIndex {
 public:
  using Int = std::ptrdiff_t;

  static constexpr Int Rank = _Rank;
  static_assert(Rank >= 0, "A tensor rank cannot be negative");

  // 3^Rank, the number of components. Formed by repeated multiplication
  // rather than by pow, which is neither constexpr nor exact.
  static constexpr Int Size = [] {
    auto size = Int{1};
    for (auto i = Int{0}; i < Rank; i++) size *= 3;
    return size;
  }();

  constexpr MultiIndex() : _slots{} {}

  constexpr explicit MultiIndex(std::array<Int, Rank> slots) : _slots{slots} {
    for (auto alpha : _slots) {
      if (alpha < -1 || alpha > 1) {
        throw std::invalid_argument(
            "A canonical index must be -1, 0 or +1");
      }
    }
  }

  // Slot 0 is the most significant digit, so that flat indices run in
  // lexicographic order of (alpha_1 ... alpha_p) starting from all -1. The
  // choice is arbitrary but has to be pinned: it is the order components are
  // laid out in and the order any table below is generated in.
  static constexpr MultiIndex FromFlat(Int flat) {
    if (flat < 0 || flat >= Size) {
      throw std::invalid_argument("Flat component index out of range");
    }
    auto slots = std::array<Int, Rank>{};
    for (auto i = Rank - 1; i >= 0; i--) {
      slots[i] = flat % 3 - 1;
      flat /= 3;
    }
    return MultiIndex(slots);
  }

  constexpr Int Flat() const {
    auto flat = Int{0};
    for (auto alpha : _slots) flat = 3 * flat + (alpha + 1);
    return flat;
  }

  // The upper index, eq:N. This is the spin weight of the component's field
  // and therefore a compile-time quantity wherever it is used.
  constexpr Int UpperIndex() const {
    auto n = Int{0};
    for (auto alpha : _slots) n += alpha;
    return n;
  }

  constexpr Int operator[](Int slot) const { return _slots[slot]; }
  constexpr const std::array<Int, Rank>& Slots() const { return _slots; }

  // The involution of the reality condition (eq:reality). Its only fixed
  // point is the all-zero index.
  constexpr MultiIndex Negated() const {
    auto slots = _slots;
    for (auto& alpha : slots) alpha = -alpha;
    return MultiIndex(slots);
  }

  // Slot i of the result is slot image[i] of this one.
  constexpr MultiIndex Permuted(const std::array<Int, Rank>& image) const {
    auto slots = std::array<Int, Rank>{};
    for (auto i = Int{0}; i < Rank; i++) slots[i] = _slots[image[i]];
    return MultiIndex(slots);
  }

  constexpr bool operator==(const MultiIndex&) const = default;

 private:
  std::array<Int, Rank> _slots;
};

// How many components of a rank-p tensor carry upper index N: the trinomial
// coefficient, which for p = 2 gives 1, 2, 3, 2, 1 (theory note section 2).
// Counted rather than derived, since the enumeration is cheap and the closed
// form is one more thing to get wrong.
template <std::ptrdiff_t Rank>
constexpr auto ComponentsAtUpperIndex(std::ptrdiff_t n) {
  auto count = std::ptrdiff_t{0};
  for (auto flat = std::ptrdiff_t{0}; flat < MultiIndex<Rank>::Size; flat++) {
    if (MultiIndex<Rank>::FromFlat(flat).UpperIndex() == n) count++;
  }
  return count;
}

//--------------------------------------------------------------------------//
//                            Permutation symmetry                           //
//--------------------------------------------------------------------------//

// A permutation of a tensor's slots, together with the sign the component
// picks up under it: +1 for a symmetric generator, -1 for an antisymmetric
// one.
template <std::ptrdiff_t Rank>
struct SlotPermutation {
  using Int = std::ptrdiff_t;

  std::array<Int, Rank> image;
  Int sign;

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

// No relation between components: every one of the 3^Rank is independent.
template <std::ptrdiff_t Rank>
struct NoSymmetry {
  static constexpr auto Generators() {
    return std::array<SlotPermutation<Rank>, 0>{};
  }
};

// Symmetric or antisymmetric in every pair of slots. The adjacent
// transpositions generate the whole symmetric group, so two generators
// suffice at any rank and the orbit search finds the rest.
template <std::ptrdiff_t Rank, std::ptrdiff_t Sign>
struct FullPermutationSymmetry {
  static_assert(Sign == 1 || Sign == -1);

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

// Anything else, given by its generators. The case this exists for is the
// elastic tensor, c_{ijkl} = c_{jikl} = c_{ijlk} = c_{klij}, which is three
// generators and is why rank 4 is exposed at all.
template <std::ptrdiff_t Rank, SlotPermutation<Rank>... Gs>
struct GeneratedBy {
  static constexpr auto Generators() {
    return std::array<SlotPermutation<Rank>, sizeof...(Gs)>{Gs...};
  }
};

// The symmetry of an elastic tensor: symmetric within each pair of slots, and
// under exchange of the pairs.
struct ElasticSymmetry {
  static constexpr std::ptrdiff_t Rank = 4;

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
