#ifndef GSH_TRANS_ORBITS_GUARD_H
#define GSH_TRANS_ORBITS_GUARD_H

#include <array>
#include <complex>
#include <cstddef>
#include <utility>

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
// a real one, and is what the reality reduction turns on. The
// algorithm does not care which generators it is given, which is the point:
// reality adds negation to the set and nothing else changes.
//
// The tables in the theory note, docs/canonical-components.tex, are the
// check. Under negation alone a real
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
// are *no pinned components at all*, which removes the second buffer reality
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

/**
 * @brief How one component is obtained from the representative stored for it.
 * @details A value rather than three loose constants so that it can be handed
 * whole, as a template argument, to the two functions that interpret it --
 * DerivedComponent and DerivedCoefficient -- and to nothing else.
 */
struct OrbitRelation {
  std::ptrdiff_t sign;  ///< The sign relating the two.
  bool conjugate;       ///< Whether the relation conjugates.
  /** @brief What the orbit pins its members to. */
  ComponentConstraint constraint;
};

/**
 * @brief For every component of a tensor: which component is actually stored
 * for it, how the two are related, and what its orbit pins it to.
 *
 * @tparam Rank_ The tensor rank.
 * @tparam Slots_ The alphabet the slots are drawn from.
 */
template <std::ptrdiff_t Rank_, SlotAlphabet Slots_ = AllSlots>
struct OrbitTable {
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  using SlotSet = Slots_;      ///< The alphabet the slots are drawn from.

  /** @brief The tensor rank. */
  static constexpr Int Rank = Rank_;
  /** @brief How many components there are, stored or not. */
  static constexpr Int Size = MultiIndex<Rank, SlotSet>::Size;

  /// For each of the 3^Rank components: the flat index of the component
  /// actually stored for it, the sign relating the two, and whether the
  /// relation conjugates. So
  ///
  ///   T^{flat} = sign[flat] * (conjugate[flat] ? conj :
  ///   id)(T^{representative[flat]}).
  std::array<Int, Size> representative{};  ///< The component actually stored.
  std::array<Int, Size> sign{};            ///< The sign relating the two.
  std::array<bool, Size> conjugate{};      ///< Whether the relation conjugates.

  /// Where a component lives in the buffer, or -1 if it is not stored -- either
  /// because another member of its orbit is, or because the orbit is
  /// identically zero.
  std::array<Int, Size> slot{};

  /// The stored components in slot order, and how many there are.
  std::array<Int, Size> stored{};
  Int storedCount{};  ///< How many components are stored.

  /// Per component, the constraint its orbit imposes. Constant along an orbit.
  std::array<ComponentConstraint, Size> constraint{};

  /** @brief The upper index of the component with that flat index. */
  constexpr Int UpperIndexOf(Int flat) const {
    return MultiIndex<Rank, SlotSet>::FromFlat(flat).UpperIndex();
  }

  /** @brief How the component with that flat index is obtained. */
  constexpr OrbitRelation RelationOf(Int flat) const {
    return {sign[flat], conjugate[flat], constraint[flat]};
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
        Reach(index.Permuted(generator.image).Flat(), fromSign * generator.sign,
              fromConjugated);
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

//--------------------------------------------------------------------------//
//                    A derived component from its representative            //
//--------------------------------------------------------------------------//

// What the table records is the relation
//
//   T = sign * (conjugate ? conj : id)(R),
//
// with R the representative's *value*. What is stored for R depends on what
// its orbit pins it to, and the three cases differ in what conjugation does:
//
//   None        R is the stored complex field,    conj(R) is its conjugate
//   Real        R = r, the stored real field,     conj(R) = R
//   Imaginary   R = i r,                          conj(R) = -R
//
// So conjugation is a sign change for an imaginary orbit and *nothing at all*
// for a real one. That distinction is the whole content of these two
// functions, and it is why there are two functions and not four accessors:
// written out separately on the flat and the layered type, and again on their
// expansions, the four copies came to disagree, each wrong for a different
// case -- a sign lost on the real orbits of an elastic tensor, which no rank-2
// test could see because no pinned member there is reached by conjugating.
// tests/TestTensorOrbitValues.cpp checks the defining relations on every
// component of every type, which is a check on this and not a copy of it.

/**
 * @brief The field a component takes, from the view of what is stored for it.
 *
 * @details The view is taken **by value** and moved into whatever is returned.
 * A spin-weighted node holds an lvalue terminal by reference, which is right
 * at a call site and wrong in an accessor whose view is a local: the node
 * would outlive it. Taking the view by value makes that mistake impossible to
 * write here rather than something each caller must remember.
 *
 * @tparam relation How the component is obtained from its representative.
 * @tparam Real The precision.
 * @param view A view of the representative's storage: complex at the
 * representative's upper index if the orbit is unpinned, real at upper index
 * zero if it is pinned.
 */
template <OrbitRelation relation, typename Real, typename View>
auto DerivedComponent(View view) {
  constexpr auto scale = static_cast<Real>(relation.sign);
  if constexpr (relation.constraint == ComponentConstraint::None) {
    if constexpr (relation.conjugate) {
      return scale * conj(std::move(view));
    } else if constexpr (relation.sign == 1) {
      return view;
    } else {
      return -std::move(view);
    }
  } else if constexpr (relation.constraint == ComponentConstraint::Real) {
    return scale * std::move(view);
  } else {
    constexpr auto turn = relation.conjugate ? -scale : scale;
    return std::complex<Real>{0, turn} * std::move(view);
  }
}

/**
 * @brief The coefficient a component takes at order m, from the coefficients
 * of what is stored for it.
 *
 * @details The spectral form of the same relation. Conjugating a field of
 * upper index N sends its coefficients to
 * @f$(-1)^{m+N}\,\overline{f_{l,-m}}@f$, and the factor i of an imaginary
 * orbit is conjugated with the field it multiplies.
 *
 * @tparam relation How the component is obtained from its representative.
 * @tparam Real The precision.
 * @param m The order wanted.
 * @param storedUpperIndex The upper index of the stored field: the
 * representative's, or zero for a pinned orbit.
 * @param stored Returns the stored field's coefficient at a given order, at
 * the degree the caller has fixed. It is asked for order m or for order -m.
 */
template <OrbitRelation relation, typename Real, typename Stored>
std::complex<Real> DerivedCoefficient(std::ptrdiff_t m,
                                      std::ptrdiff_t storedUpperIndex,
                                      Stored&& stored) {
  using Complex = std::complex<Real>;
  constexpr auto scale = static_cast<Real>(relation.sign);
  constexpr auto turn = relation.constraint == ComponentConstraint::Imaginary
                            ? Complex{0, 1}
                            : Complex{1, 0};
  if constexpr (relation.conjugate) {
    return scale * std::conj(turn) *
           static_cast<Real>(MinusOneToPower(m + storedUpperIndex)) *
           std::conj(Complex{stored(-m)});
  } else {
    return scale * turn * Complex{stored(m)};
  }
}

/// The stored-component table for a tensor: permutation symmetry always, the
/// reality condition only for a real tensor.
///
/// This is the switch the reality reduction turns on. A real tensor's
/// components are related by
/// T^{-alpha} = (-1)^N conj(T^alpha), so the negation joins the generating set
/// and the orbits get larger and fewer; a complex tensor has no such relation
/// and keeps every component of every permutation orbit.
struct RealTensor {
  /** @brief The negation joins the generating set. */
  static constexpr bool ReducesOnReality = true;
};

/// A tensor with no reality relation between its components.
struct ComplexTensor {
  /** @brief Every component of every permutation orbit is kept. */
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
