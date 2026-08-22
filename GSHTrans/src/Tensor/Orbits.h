#ifndef GSH_TRANS_ORBITS_GUARD_H
#define GSH_TRANS_ORBITS_GUARD_H

#include <array>
#include <cstddef>

#include "../Utility.h"
#include "MultiIndex.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                        Which components must be stored                    //
//--------------------------------------------------------------------------//

// A tensor's components are not independent. Two things relate them, and this
// file treats both with one algorithm because they compose:
//
//   permutation  T^{pi(alpha)} = s T^{alpha},  s = +1 symmetric, -1 anti
//   negation     T^{-alpha}    = (-1)^N conj(T^{alpha})        (eq:reality)
//
// Together they generate a group acting on the 3^Rank multi-indices, and the
// components that must be stored are one representative of each orbit. The
// permutation half applies to every tensor; the negation half applies only to
// a real one, and is what phase 4 of the field-algebra plan turns on. The
// algorithm does not care which generators it is given, which is the point:
// phase 4 adds negation to the set and nothing else changes.
//
// The theory note's own tables are the check. Under negation alone a real
// rank-2 tensor stores 5 of its 9 components and a rank-4 one stores 41 of 81;
// a symmetric real rank-2 tensor has four orbits carrying six reals per point,
// which is the number of independent entries of a real symmetric 3x3 matrix.
//
// The same algorithm serves a tangential tensor, whose slots run over {-1, +1}
// rather than {-1, 0, +1}, and it is told nothing about the difference beyond
// which multi-indices exist. Two answers fall out that are worth knowing in
// advance, because they are the check that this is a generalisation and not a
// special case: with no permutation symmetry, negation has no fixed point --
// the all-zero index does not exist -- so every orbit has size two and there
// are *no pinned components at all*, which removes the second buffer phase 4
// had to introduce. Under a symmetry there can still be one: for a symmetric
// tangential rank-2 tensor, negation maps (-+) to (+-) and the symmetry maps
// it back, so that component is pinned real and the tensor is one complex
// number plus one real, three reals a point -- a real symmetric 2x2 matrix.

// What an orbit forces on its own representative.
//
// Zero arises without reality: a component of an antisymmetric tensor that a
// permutation maps to itself with sign -1 satisfies T = -T. Real and Imaginary
// arise only once negation is in the group, where a component related to its
// own conjugate is pinned to one or the other -- the self-paired components of
// the theory note, which always sit at N = 0 because permutation preserves the
// slot sum while negation reverses it.
enum class ComponentConstraint { None, Zero, Real, Imaginary };

template <std::ptrdiff_t _Rank, SlotAlphabet _Slots = AllSlots>
struct OrbitTable {
  using Int = std::ptrdiff_t;
  using SlotSet = _Slots;

  static constexpr Int Rank = _Rank;
  static constexpr Int Size = MultiIndex<Rank, SlotSet>::Size;

  // For each of the 3^Rank components: the flat index of the component
  // actually stored for it, the sign relating the two, and whether the
  // relation conjugates. So
  //
  //   T^{flat} = sign[flat] * (conjugate[flat] ? conj : id)(T^{representative[flat]}).
  std::array<Int, Size> representative{};
  std::array<Int, Size> sign{};
  std::array<bool, Size> conjugate{};

  // Where a component lives in the buffer, or -1 if it is not stored -- either
  // because another member of its orbit is, or because the orbit is
  // identically zero.
  std::array<Int, Size> slot{};

  // The stored components in slot order, and how many there are.
  std::array<Int, Size> stored{};
  Int storedCount{};

  // Per component, the constraint its orbit imposes. Constant along an orbit.
  std::array<ComponentConstraint, Size> constraint{};

  constexpr Int UpperIndexOf(Int flat) const {
    return MultiIndex<Rank, SlotSet>::FromFlat(flat).UpperIndex();
  }
};

// Build the table for a rank, a permutation symmetry, and a choice of whether
// the reality condition is one of the generators.
//
// A breadth-first walk of each orbit. Reaching a component twice by different
// routes is not a failure but the interesting case: the two routes give two
// expressions for the same component in terms of the representative, and
// equating them is what produces the constraint.
//
//   same conjugation, opposite signs   ->  T = -T,        so the orbit is zero
//   opposite conjugation, signs equal  ->  conj(T) = T,   so it is real
//   opposite conjugation, signs differ ->  conj(T) = -T,  so it is imaginary
//
// The sign a negation contributes is (-1)^N, and N belongs to the component
// being transformed. That is well defined along any route because permutation
// preserves N, negation reverses it, and (-1)^N = (-1)^{-N}.
template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          bool IncludeNegation, SlotAlphabet Slots = AllSlots>
constexpr auto MakeOrbitTable() {
  using Int = std::ptrdiff_t;
  using Index = MultiIndex<Rank, Slots>;
  constexpr auto Size = Index::Size;

  auto table = OrbitTable<Rank, Slots>{};
  const auto generators = Symmetry::Generators();

  for (auto& r : table.representative) r = -1;
  for (auto& s : table.slot) s = -1;

  for (auto start = Int{0}; start < Size; start++) {
    if (table.representative[start] != -1) continue;

    auto orbitConstraint = ComponentConstraint::None;

    // The orbit, as a queue of (component, sign, conjugated) triples. It can
    // never hold more than the whole index set.
    auto queue = std::array<Int, Size>{};
    auto head = Int{0};
    auto tail = Int{0};

    table.representative[start] = start;
    table.sign[start] = 1;
    table.conjugate[start] = false;
    queue[tail++] = start;

    // What a route says about a component, given what is already recorded.
    auto Reach = [&](Int to, Int sign, bool conjugated) {
      if (table.representative[to] == -1) {
        table.representative[to] = start;
        table.sign[to] = sign;
        table.conjugate[to] = conjugated;
        queue[tail++] = to;
        return;
      }
      if (table.conjugate[to] == conjugated) {
        if (table.sign[to] != sign) orbitConstraint = ComponentConstraint::Zero;
      } else if (orbitConstraint != ComponentConstraint::Zero) {
        orbitConstraint = table.sign[to] * sign > 0
                              ? ComponentConstraint::Real
                              : ComponentConstraint::Imaginary;
      }
    };

    while (head != tail) {
      const auto from = queue[head++];
      const auto index = Index::FromFlat(from);
      const auto fromSign = table.sign[from];
      const auto fromConjugated = table.conjugate[from];

      for (const auto& generator : generators) {
        Reach(index.Permuted(generator.image).Flat(),
              fromSign * generator.sign, fromConjugated);
      }

      if constexpr (IncludeNegation) {
        Reach(index.Negated().Flat(),
              fromSign * MinusOneToPower(index.UpperIndex()), !fromConjugated);
      }
    }

    for (auto i = Int{0}; i < tail; i++) {
      table.constraint[queue[i]] = orbitConstraint;
    }

    if (orbitConstraint != ComponentConstraint::Zero) {
      table.slot[start] = table.storedCount;
      table.stored[table.storedCount++] = start;
    }
  }

  return table;
}

// The stored-component table for a tensor: permutation symmetry always, the
// reality condition only for a real tensor.
//
// This is the switch phase 2 left for phase 4, and turning it on was the whole
// of the change to this file. A real tensor's components are related by
// T^{-alpha} = (-1)^N conj(T^alpha), so the negation joins the generating set
// and the orbits get larger and fewer; a complex tensor has no such relation
// and keeps every component of every permutation orbit.
struct RealTensor {
  static constexpr bool ReducesOnReality = true;
};

struct ComplexTensor {
  static constexpr bool ReducesOnReality = false;
};

template <typename R>
concept TensorReality =
    std::same_as<R, RealTensor> or std::same_as<R, ComplexTensor>;

template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, SlotAlphabet Slots = AllSlots>
inline constexpr auto TensorOrbits =
    MakeOrbitTable<Rank, Symmetry, Reality::ReducesOnReality, Slots>();

}  // namespace GSHTrans

#endif  // GSH_TRANS_ORBITS_GUARD_H
