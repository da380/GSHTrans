#include <gtest/gtest.h>

#include <GSHTrans/All>

#include <array>
#include <cstddef>

namespace {

using namespace GSHTrans;
using Int = std::ptrdiff_t;

// The reals per grid point a stored set costs: two for each stored component,
// less one for each that a constraint pins to a single real number.
//
// This is the invariant worth testing above any individual count, because it
// is what the abstraction is *for*. It must come out equal to the number of
// independent real degrees of freedom of the tensor, and it does so for every
// rank and symmetry below without any of them being special-cased.
template <Int Rank, typename Symmetry, bool Reality>
constexpr auto RealsPerPoint() {
  constexpr auto table = MakeOrbitTable<Rank, Symmetry, Reality>();
  auto reals = Int{0};
  for (auto i = Int{0}; i < table.Size; i++) {
    if (table.slot[i] < 0) continue;
    reals += table.constraint[i] == ComponentConstraint::None ? 2 : 1;
  }
  return reals;
}

template <Int Rank, typename Symmetry, bool Reality>
constexpr auto StoredCount() {
  return MakeOrbitTable<Rank, Symmetry, Reality>().storedCount;
}

}  // namespace

//--------------------------------------------------------------------------//
//                               The multi-index                             //
//--------------------------------------------------------------------------//

TEST(MultiIndex, FlatIndexRoundTrips) {
  constexpr auto check = []<Int Rank>() {
    for (auto flat = Int{0}; flat < MultiIndex<Rank>::Size; flat++) {
      if (MultiIndex<Rank>::FromFlat(flat).Flat() != flat) return false;
    }
    return true;
  };
  static_assert(check.template operator()<0>());
  static_assert(check.template operator()<1>());
  static_assert(check.template operator()<2>());
  static_assert(check.template operator()<3>());
  static_assert(check.template operator()<4>());
  SUCCEED();
}

TEST(MultiIndex, SizeIsThreeToTheRank) {
  static_assert(MultiIndex<0>::Size == 1);
  static_assert(MultiIndex<1>::Size == 3);
  static_assert(MultiIndex<2>::Size == 9);
  static_assert(MultiIndex<4>::Size == 81);
  SUCCEED();
}

// Theory note table 1: the upper index of each component of a rank-2 tensor.
// The point of the table is that the multi-index and the upper index are
// different things -- three distinct components share N = 0.
TEST(MultiIndex, UpperIndexIsTheSignedSum) {
  constexpr auto Index = [](Int a, Int b) {
    return MultiIndex<2>(std::array<Int, 2>{a, b});
  };

  static_assert(Index(-1, -1).UpperIndex() == -2);
  static_assert(Index(-1, 0).UpperIndex() == -1);
  static_assert(Index(0, -1).UpperIndex() == -1);
  static_assert(Index(-1, 1).UpperIndex() == 0);
  static_assert(Index(0, 0).UpperIndex() == 0);
  static_assert(Index(1, -1).UpperIndex() == 0);
  static_assert(Index(0, 1).UpperIndex() == 1);
  static_assert(Index(1, 0).UpperIndex() == 1);
  static_assert(Index(1, 1).UpperIndex() == 2);

  // The count at each upper index is the trinomial coefficient: 1, 2, 3, 2, 1.
  static_assert(ComponentsAtUpperIndex<2>(-2) == 1);
  static_assert(ComponentsAtUpperIndex<2>(-1) == 2);
  static_assert(ComponentsAtUpperIndex<2>(0) == 3);
  static_assert(ComponentsAtUpperIndex<2>(1) == 2);
  static_assert(ComponentsAtUpperIndex<2>(2) == 1);

  // At rank 1 they coincide, which is the coincidence that misleads.
  static_assert(MultiIndex<1>(std::array<Int, 1>{-1}).UpperIndex() == -1);
  static_assert(MultiIndex<1>(std::array<Int, 1>{1}).UpperIndex() == 1);
  SUCCEED();
}

TEST(MultiIndex, NegationIsAnInvolutionWithOneFixedPoint) {
  constexpr auto fixedPoints = []<Int Rank>() {
    auto count = Int{0};
    for (auto flat = Int{0}; flat < MultiIndex<Rank>::Size; flat++) {
      const auto index = MultiIndex<Rank>::FromFlat(flat);
      if (index.Negated().Negated() != index) return Int{-1};
      if (index.Negated() == index) count++;
    }
    return count;
  };
  static_assert(fixedPoints.template operator()<1>() == 1);
  static_assert(fixedPoints.template operator()<2>() == 1);
  static_assert(fixedPoints.template operator()<4>() == 1);
  SUCCEED();
}

//--------------------------------------------------------------------------//
//                                  Orbits                                   //
//--------------------------------------------------------------------------//

// Theory note table 2, the storage table for a real tensor under negation
// alone. These are the numbers the reduction is claimed to achieve.
TEST(Orbits, ReproducesTheStorageTableOfTheTheoryNote) {
  static_assert((StoredCount<0, NoSymmetry<0>, true>()) == 1);
  static_assert((StoredCount<1, NoSymmetry<1>, true>()) == 2);
  static_assert((StoredCount<2, NoSymmetry<2>, true>()) == 5);
  static_assert((StoredCount<4, NoSymmetry<4>, true>()) == 41);

  // Exactly one of them is real-valued in each case: the all-zero component,
  // which is the only fixed point of the negation.
  constexpr auto realCount = []<Int Rank>() {
    constexpr auto table = MakeOrbitTable<Rank, NoSymmetry<Rank>, true>();
    auto count = Int{0};
    for (auto i = Int{0}; i < table.Size; i++) {
      if (table.slot[i] >= 0 &&
          table.constraint[i] == ComponentConstraint::Real) {
        count++;
      }
    }
    return count;
  };
  static_assert(realCount.template operator()<1>() == 1);
  static_assert(realCount.template operator()<2>() == 1);
  static_assert(realCount.template operator()<4>() == 1);

  // And the storage in reals is 3^p, the real degrees of freedom of a real
  // rank-p tensor. "The reduction loses nothing."
  static_assert((RealsPerPoint<1, NoSymmetry<1>, true>()) == 3);
  static_assert((RealsPerPoint<2, NoSymmetry<2>, true>()) == 9);
  static_assert((RealsPerPoint<4, NoSymmetry<4>, true>()) == 81);
  SUCCEED();
}

// The theory note's worked check, which it offers as evidence that the orbit
// construction is the right abstraction. Four orbits, six reals, and the
// membership is what makes it a check rather than a coincidence of counts.
TEST(Orbits, ReproducesTheWorkedSymmetricRankTwoExample) {
  constexpr auto table = MakeOrbitTable<2, Symmetric<2>, true>();
  constexpr auto Flat = [](Int a, Int b) {
    return MultiIndex<2>(std::array<Int, 2>{a, b}).Flat();
  };

  static_assert(table.storedCount == 4);
  static_assert((RealsPerPoint<2, Symmetric<2>, true>()) == 6);

  // {(00)} alone, and real.
  static_assert(table.constraint[Flat(0, 0)] == ComponentConstraint::Real);

  // {(-+), (+-)}: one orbit, self-paired, and therefore also real.
  static_assert(table.representative[Flat(-1, 1)] ==
                table.representative[Flat(1, -1)]);
  static_assert(table.constraint[Flat(-1, 1)] == ComponentConstraint::Real);

  // {(--), (++)}: one orbit, complex.
  static_assert(table.representative[Flat(-1, -1)] ==
                table.representative[Flat(1, 1)]);
  static_assert(table.constraint[Flat(-1, -1)] == ComponentConstraint::None);

  // {(-0), (0-), (+0), (0+)}: all four together, complex.
  static_assert(table.representative[Flat(-1, 0)] ==
                table.representative[Flat(0, -1)]);
  static_assert(table.representative[Flat(-1, 0)] ==
                table.representative[Flat(1, 0)]);
  static_assert(table.representative[Flat(-1, 0)] ==
                table.representative[Flat(0, 1)]);
  static_assert(table.constraint[Flat(-1, 0)] == ComponentConstraint::None);

  // Self-paired components always sit at N = 0: permutation preserves the
  // slot sum and negation reverses it, so a component fixed by the pair
  // satisfies N = -N.
  constexpr auto selfPairedAtZero = [](const auto& t) {
    for (auto i = Int{0}; i < t.Size; i++) {
      if (t.constraint[i] == ComponentConstraint::Real ||
          t.constraint[i] == ComponentConstraint::Imaginary) {
        if (MultiIndex<2>::FromFlat(i).UpperIndex() != 0) return false;
      }
    }
    return true;
  };
  static_assert(selfPairedAtZero(table));
  SUCCEED();
}

// Antisymmetry produces components that vanish identically, which is the one
// kind of constraint that arises without the reality condition. Worth pinning
// because it is the case a design assuming "every component is stored or
// derived from a stored one" gets wrong.
TEST(Orbits, AntisymmetryAnnihilatesTheDiagonal) {
  constexpr auto table = MakeOrbitTable<2, Antisymmetric<2>, false>();
  constexpr auto Flat = [](Int a, Int b) {
    return MultiIndex<2>(std::array<Int, 2>{a, b}).Flat();
  };

  static_assert(table.constraint[Flat(-1, -1)] == ComponentConstraint::Zero);
  static_assert(table.constraint[Flat(0, 0)] == ComponentConstraint::Zero);
  static_assert(table.constraint[Flat(1, 1)] == ComponentConstraint::Zero);
  static_assert(table.slot[Flat(0, 0)] == -1);

  // Three independent components, which is the count for an antisymmetric
  // 3x3 matrix, and three reals once it is real.
  static_assert(table.storedCount == 3);
  static_assert((RealsPerPoint<2, Antisymmetric<2>, true>()) == 3);
  SUCCEED();
}

// The case rank 4 is exposed for. The theory note does not mention it, which
// makes it the strongest check here: nothing in the machinery knows about
// elasticity, and 21 is the answer everyone already knows.
TEST(Orbits, ElasticSymmetryGivesTwentyOneIndependentComponents) {
  static_assert((StoredCount<4, ElasticSymmetry, false>()) == 21);
  static_assert((RealsPerPoint<4, ElasticSymmetry, true>()) == 21);
  SUCCEED();
}

// Whatever the group, the table has to describe a consistent set of
// relations: every component points at a stored representative, every
// representative points at itself, and no component outside a vanishing orbit
// is left without storage.
TEST(Orbits, TheTableIsInternallyConsistent) {
  constexpr auto consistent = []<Int Rank, typename Symmetry, bool Reality>() {
    constexpr auto table = MakeOrbitTable<Rank, Symmetry, Reality>();
    for (auto i = Int{0}; i < table.Size; i++) {
      const auto rep = table.representative[i];
      if (rep < 0 || rep >= table.Size) return false;
      if (table.representative[rep] != rep) return false;
      if (table.constraint[rep] != table.constraint[i]) return false;
      if (table.sign[i] != 1 && table.sign[i] != -1) return false;

      const auto vanishes = table.constraint[i] == ComponentConstraint::Zero;
      if (vanishes != (table.slot[rep] < 0)) return false;
      if (!vanishes && table.stored[table.slot[rep]] != rep) return false;
      if (i != rep && table.slot[i] >= 0) return false;

      // A relation that conjugates can only appear when reality is in the
      // generating set.
      if (table.conjugate[i] && !Reality) return false;
    }
    return true;
  };

  static_assert(consistent.template operator()<2, NoSymmetry<2>, false>());
  static_assert(consistent.template operator()<2, NoSymmetry<2>, true>());
  static_assert(consistent.template operator()<2, Symmetric<2>, true>());
  static_assert(consistent.template operator()<2, Antisymmetric<2>, true>());
  static_assert(consistent.template operator()<3, Symmetric<3>, true>());
  static_assert(consistent.template operator()<4, NoSymmetry<4>, true>());
  static_assert(consistent.template operator()<4, ElasticSymmetry, true>());
  SUCCEED();
}

// A real tensor reduces on both relations, a complex one on permutation
// symmetry alone. This is the switch phase 2 left for phase 4.
TEST(Orbits, RealityReducesOnlyARealTensor) {
  static_assert(RealTensor::ReducesOnReality);
  static_assert(!ComplexTensor::ReducesOnReality);

  static_assert((TensorOrbits<2, NoSymmetry<2>, ComplexTensor>.storedCount) ==
                9);
  static_assert((TensorOrbits<2, NoSymmetry<2>, RealTensor>.storedCount) == 5);

  static_assert((TensorOrbits<2, Symmetric<2>, ComplexTensor>.storedCount) ==
                6);
  static_assert((TensorOrbits<2, Symmetric<2>, RealTensor>.storedCount) == 4);

  static_assert((TensorOrbits<4, NoSymmetry<4>, RealTensor>.storedCount) == 41);
  static_assert((TensorOrbits<4, ElasticSymmetry, RealTensor>.storedCount) ==
                13);
  SUCCEED();
}
