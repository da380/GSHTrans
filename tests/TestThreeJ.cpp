#include <gtest/gtest.h>

#include <GSHTrans/Core.hpp>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include "RacahReference.hpp"

// The test family for 3j.hpp.
//
// `3j.hpp` had no coverage at all although it is public API, and the reason
// this file can be written without a reference implementation is the
// completeness relation:
//
//     sum over the (m1, m3) plane of the squared symbols  =  1
//
// exactly, for any triple satisfying the triangle rule. That is one loop over
// a table that has just been built, it needs nothing to compare against, and
// it is what the rest of this file is built on.
//
// Completeness cannot see everything -- Schulten-Gordon normalises each row
// by it, so it holds by construction -- and the tests further down say what
// each of the other checks can and cannot see. The stretched region, where
// the recursion once ran in its unstable direction, is covered at the end of
// the file against Racah's closed form, which is exact there.

namespace {

using namespace GSHTrans;

// Sum of squares over the whole (m1, m3) plane at fixed degrees.
template <typename T>
T Completeness(int l1, int l2, int l3) {
  const auto table = Wigner3jMatrix<T>(l1, l2, l3);
  auto sum = T{0};
  for (auto m1 : table.M1Axis()) {
    for (auto m3 : table.M3Axis()) {
      const auto value = table(m1, m3);
      sum += value * value;
    }
  }
  return sum;
}

//--------------------------------------------------------------------------//
//                     Where the identity holds                              //
//--------------------------------------------------------------------------//

TEST(ThreeJ, CompletenessHoldsForFatTriangles) {
  for (auto l : {1, 2, 3, 5, 8, 13, 20, 32}) {
    EXPECT_NEAR(Completeness<double>(l, l, l), 1.0, 1e-12)
        << "(l,l,l), l = " << l;
  }
}

TEST(ThreeJ, CompletenessHoldsAwayFromStretched) {
  // A spread of shapes, none of them close to l3 = l1 + l2.
  const auto triples = std::vector<std::array<int, 3>>{
      {4, 4, 2},  {6, 4, 4},    {10, 7, 5},   {12, 12, 6},
      {16, 9, 9}, {20, 20, 10}, {24, 16, 12}, {30, 20, 16}};
  for (const auto& t : triples) {
    EXPECT_NEAR(Completeness<double>(t[0], t[1], t[2]), 1.0, 1e-12)
        << "(" << t[0] << "," << t[1] << "," << t[2] << ")";
  }
}

// The degenerate edge of the triangle rule, and where the present recursion
// starts to lose it. Measured, the departure from one is:
//
//     l = 16 : 1e-16     l = 25 : 1.8e-08
//     l = 20 : 4.4e-14   l = 28 : 1.2e-05
//                        l = 30 : 4.1e-03
//
// so the failure is **exponential from about l = 20** rather than a cliff at
// 30 rather than a cliff at 30. Two assertions
// therefore, at two tolerances, so that the shape of the decay is pinned and
// not just its ends.
TEST(ThreeJ, CompletenessHoldsAtTheStretchedEdgeForModestDegrees) {
  for (auto l : {1, 2, 4, 8, 12, 16, 20}) {
    EXPECT_NEAR(Completeness<double>(l, l, 2 * l), 1.0, 1e-12)
        << "(l, l, 2l), l = " << l;
  }
  for (auto l : {22, 25}) {
    EXPECT_NEAR(Completeness<double>(l, l, 2 * l), 1.0, 1e-6)
        << "(l, l, 2l), l = " << l << " is degrading but still usable";
  }
}

//--------------------------------------------------------------------------//
//                     Where the identity fails, and it must                 //
//--------------------------------------------------------------------------//

// The old boundary is gone. Schulten-Gordon recurses inward from both
// forbidden ends and matches in the middle, so the region that destroyed the
// one-directional scheme -- and the band that neither it nor Racah's closed
// form could reach -- is answered like anywhere else.
TEST(ThreeJ, AnswersEverythingTheOldSchemesCouldNot) {
  // Stretched: fatal to the old recursion past l = 30.
  for (auto l : {35, 40, 50, 64, 100, 128}) {
    ASSERT_NO_THROW(Wigner3jMatrix<double>(l, l, 2 * l)) << "l = " << l;
    EXPECT_NEAR(Completeness<double>(l, l, 2 * l), 1.0, 1e-12) << "l = " << l;
  }
  // The band neither classical method reached.
  for (const auto& t : std::vector<std::array<int, 3>>{{80, 80, 120},
                                                       {90, 90, 135},
                                                       {100, 100, 150},
                                                       {128, 128, 192},
                                                       {128, 128, 160},
                                                       {160, 160, 160},
                                                       {200, 200, 200}}) {
    ASSERT_NO_THROW(Wigner3jMatrix<double>(t[0], t[1], t[2]))
        << "(" << t[0] << "," << t[1] << "," << t[2] << ")";
    EXPECT_NEAR(Completeness<double>(t[0], t[1], t[2]), 1.0, 1e-12)
        << "(" << t[0] << "," << t[1] << "," << t[2] << ")";
  }
}

// The runtime check is the recurrence residual, not completeness.
// Completeness cannot fail here -- every row is normalised by that identity --
// so it is asserted above as a property and relied on nowhere.
//
// This asserts the other half: that nothing good is refused, over the whole
// range the suite exercises and a non-triangle besides.
TEST(ThreeJ, RefusesNothingItShouldAnswer) {
  for (auto l : {0, 1, 2, 3, 5, 8, 13, 20, 32, 64, 128}) {
    EXPECT_NO_THROW(Wigner3jMatrix<double>(l, l, l)) << "(l,l,l), l = " << l;
  }
  for (auto l : {1, 2, 4, 8, 16, 32, 64}) {
    EXPECT_NO_THROW(Wigner3jMatrix<double>(l, l, 2 * l));
    EXPECT_NO_THROW(Wigner3jMatrix<double>(l, l, 3 * l / 2 + 1));
  }
  for (const auto& t : std::vector<std::array<int, 3>>{
           {128, 128, 250}, {100, 150, 200}, {60, 60, 90}, {70, 70, 105}}) {
    EXPECT_NO_THROW(Wigner3jMatrix<double>(t[0], t[1], t[2]));
  }
  // A non-triangle has an identically zero table, and that is not a failure.
  EXPECT_NO_THROW(Wigner3jMatrix<double>(1, 1, 5));
}

//--------------------------------------------------------------------------//
//                              Closed forms                                 //
//--------------------------------------------------------------------------//

// (j j 0; m -m 0) = (-1)^{j-m} / sqrt(2j+1). Measured exact -- zero absolute
// difference -- so the tolerance discriminates a wrong phase or a wrong index
// and nothing else.
TEST(ThreeJ, MatchesTheClosedFormWithAZeroDegree) {
  for (auto j : {0, 1, 2, 5, 13, 30, 64}) {
    const auto table = Wigner3jMatrix<double>(j, j, 0);
    for (auto m : table.M1Axis()) {
      const auto want =
          ((j - m) % 2 == 0 ? 1.0 : -1.0) / std::sqrt(2.0 * j + 1);
      EXPECT_NEAR(table(m, 0), want, 1e-15) << "j = " << j << ", m = " << m;
    }
  }
}

// The fully stretched symbol, which is Racah's one-term case and collapses to
// something very simple:
//
//     (l1 l2 l1+l2 ; l1 l2 -(l1+l2))  =  1 / sqrt(2 l3 + 1)
//
// Every factorial in the general formula cancels at that corner, which is
// what makes this a strong test: there is no arithmetic in the expected value
// to be wrong in the same way the code is. It is also the case Racah's
// path must reproduce, since the sum there has exactly one term.
TEST(ThreeJ, MatchesTheClosedFormAtTheStretchedCorner) {
  for (auto l1 : {1, 2, 3, 5, 8, 16}) {
    for (auto l2 : {1, 2, 4, 7, 11}) {
      const auto l3 = l1 + l2;
      const auto want = 1 / std::sqrt(2.0 * l3 + 1);
      const auto table = Wigner3jMatrix<double>(l1, l2, l3);
      EXPECT_NEAR(table(l1, l2, -l3), want, 1e-13)
          << "l1 = " << l1 << ", l2 = " << l2;
    }
  }
}

//--------------------------------------------------------------------------//
//                          Selection rules                                  //
//--------------------------------------------------------------------------//

TEST(ThreeJ, VanishesWhenTheSelectionRulesAreBroken) {
  const auto table = Wigner3jMatrix<double>(3, 4, 5);

  // Orders must sum to zero.
  EXPECT_EQ(table(1, 1, 1), 0.0);
  // An order outside its degree.
  EXPECT_EQ(table(7, 0), 0.0);
  EXPECT_EQ(table(0, 9), 0.0);
  // And a triple that is not a triangle has no symbols at all.
  EXPECT_FALSE(SatisfiesTriangle(1, 1, 5));
  EXPECT_EQ(Wigner3jSymbol<double>(1, 1, 5, 0, 0, 0), 0.0);
}

// The (l1 l2 l3; 0 0 0) family vanishes unless the degrees sum to an even
// number. This is the case Gaunt integrals lean on hardest, so it is worth
// pinning separately from the general rule.
TEST(ThreeJ, VanishesAtAllZeroOrdersWhenTheDegreeSumIsOdd) {
  for (auto l1 = 0; l1 <= 6; ++l1) {
    for (auto l2 = 0; l2 <= 6; ++l2) {
      for (auto l3 = std::abs(l1 - l2); l3 <= l1 + l2; ++l3) {
        const auto value = Wigner3jSymbol<double>(l1, l2, l3, 0, 0, 0);
        if ((l1 + l2 + l3) % 2 != 0) {
          EXPECT_NEAR(value, 0.0, 1e-15)
              << "(" << l1 << "," << l2 << "," << l3 << ")";
        }
      }
    }
  }
}

//--------------------------------------------------------------------------//
//                              Symmetries                                   //
//--------------------------------------------------------------------------//

// Reversing every order multiplies the symbol by (-1)^{l1+l2+l3}.
TEST(ThreeJ, ObeysTheOrderReversalSymmetry) {
  for (const auto& t : std::vector<std::array<int, 3>>{
           {3, 4, 5}, {6, 6, 6}, {7, 5, 4}, {8, 3, 9}}) {
    const auto table = Wigner3jMatrix<double>(t[0], t[1], t[2]);
    const auto phase = (t[0] + t[1] + t[2]) % 2 == 0 ? 1.0 : -1.0;
    for (auto m1 : table.M1Axis()) {
      for (auto m3 : table.M3Axis()) {
        EXPECT_NEAR(table(-m1, -m3), phase * table(m1, m3), 1e-14)
            << "(" << t[0] << "," << t[1] << "," << t[2] << ") at m1 = " << m1
            << ", m3 = " << m3;
      }
    }
  }
}

// An even permutation of the columns leaves the symbol unchanged. Taken as
// (l1 l2 l3) -> (l2 l3 l1), with the orders carried with them.
TEST(ThreeJ, IsInvariantUnderAnEvenColumnPermutation) {
  const auto a = Wigner3jMatrix<double>(4, 5, 6);
  const auto b = Wigner3jMatrix<double>(5, 6, 4);

  for (auto m1 : a.M1Axis()) {
    for (auto m3 : a.M3Axis()) {
      const auto m2 = -(m1 + m3);
      if (std::abs(m2) > 5) continue;
      EXPECT_NEAR(b(m2, m1), a(m1, m3), 1e-14)
          << "m1 = " << m1 << ", m3 = " << m3;
    }
  }
}

//--------------------------------------------------------------------------//
//                 The stack, and the convenience entry point                //
//--------------------------------------------------------------------------//

TEST(ThreeJ, TheStackAgreesWithItsMatrices) {
  const auto stack = Wigner3jStack<double>(4, 6);
  for (auto l2 : stack.L2Axis()) {
    const auto direct = Wigner3jMatrix<double>(4, l2, 6);
    for (auto m1 : direct.M1Axis()) {
      for (auto m3 : direct.M3Axis()) {
        EXPECT_EQ(stack(l2, m1, m3), direct(m1, m3))
            << "l2 = " << l2 << ", m1 = " << m1 << ", m3 = " << m3;
      }
    }
  }
}

TEST(ThreeJ, TheSingleSymbolEntryPointAgreesWithTheTable) {
  const auto table = Wigner3jMatrix<double>(3, 5, 4);
  for (auto m1 : table.M1Axis()) {
    for (auto m3 : table.M3Axis()) {
      const auto m2 = -(m1 + m3);
      if (std::abs(m2) > 5) continue;
      EXPECT_EQ(Wigner3jSymbol<double>(3, 5, 4, m1, m2, m3), table(m1, m3))
          << "m1 = " << m1 << ", m3 = " << m3;
    }
  }
}

//--------------------------------------------------------------------------//
//        the structural checks, which carry the weight now            //
//--------------------------------------------------------------------------//
//
// With one implementation there is no second one to compare against, so the
// suite rests on properties a single implementation cannot satisfy by
// accident. The sharpest is column-permutation invariance: the recursion runs
// over m2 at fixed m1, so permuting the columns makes it run along entirely
// different lines through different data -- and it is independent of the
// per-row normalisation, which is what makes it stronger here than
// completeness.
//
// This is what caught the one real bug in the implementation. The phase was
// recovered from the last stored value, and at near-stretched triples of high
// degree the row's dynamic range reaches 1e201, so the rescaling flushed that
// element to zero and the sign with it. Whole rows came out negated with
// every magnitude correct to 1e-16. Completeness is a sum of squares and saw
// nothing; the recurrence is homogeneous and saw nothing; only this saw it.
TEST(ThreeJ, IsInvariantUnderBothCyclicPermutations) {
  for (const auto& t : std::vector<std::array<int, 3>>{{3, 4, 5},
                                                       {12, 12, 12},
                                                       {40, 40, 80},
                                                       {80, 80, 120},
                                                       {100, 100, 150},
                                                       {128, 128, 192},
                                                       {160, 160, 160},
                                                       {100, 150, 200},
                                                       {128, 128, 250}}) {
    const auto a = Wigner3jMatrix<double>(t[0], t[1], t[2]);
    const auto b = Wigner3jMatrix<double>(t[1], t[2], t[0]);
    const auto c = Wigner3jMatrix<double>(t[2], t[0], t[1]);

    auto worst = 0.0;
    for (auto m1 : a.M1Axis()) {
      for (auto m3 : a.M3Axis()) {
        const auto m2 = -(m1 + m3);
        if (std::abs(m2) > t[1]) continue;
        const auto value = a(m1, m3);
        // (l1 l2 l3; m1 m2 m3) = (l2 l3 l1; m2 m3 m1) = (l3 l1 l2; m3 m1 m2)
        worst = std::max(worst, std::abs(b(m2, m1) - value));
        worst = std::max(worst, std::abs(c(m3, m2) - value));
      }
    }
    EXPECT_LT(worst, 1e-13)
        << "(" << t[0] << "," << t[1] << "," << t[2] << ") worst " << worst;
  }
}

// Exact rational values at small degree, written out as literals. These pin
// the convention and the phase rather than the accuracy -- there is no
// arithmetic here that could be wrong in the same way the code is.
TEST(ThreeJ, MatchesExactValuesAtSmallDegree) {
  const auto third = 1.0 / 3;
  const auto t111 = Wigner3jMatrix<double>(1, 1, 1);
  // (1 1 1; 1 -1 0) = 1/sqrt(6), (1 1 1; 0 0 0) = 0.
  EXPECT_NEAR(t111(1, 0), 1 / std::sqrt(6.0), 1e-15);
  EXPECT_NEAR(t111(0, 0), 0.0, 1e-15);

  const auto t112 = Wigner3jMatrix<double>(1, 1, 2);
  // (1 1 2; 0 0 0) = sqrt(2/15), (1 1 2; 1 1 -2) = 1/sqrt(5).
  EXPECT_NEAR(t112(0, 0), std::sqrt(2.0 / 15), 1e-15);
  EXPECT_NEAR(t112(1, -2), 1 / std::sqrt(5.0), 1e-15);

  const auto t222 = Wigner3jMatrix<double>(2, 2, 2);
  // (2 2 2; 0 0 0) = -sqrt(2/35).
  EXPECT_NEAR(t222(0, 0), -std::sqrt(2.0 / 35), 1e-15);

  const auto t110 = Wigner3jMatrix<double>(1, 1, 0);
  // (1 1 0; m -m 0) = (-1)^(1-m)/sqrt(3).
  EXPECT_NEAR(t110(1, 0), third * std::sqrt(3.0), 1e-15);
  EXPECT_NEAR(t110(0, 0), -third * std::sqrt(3.0), 1e-15);
}

// The recurrence residual is the runtime check, so it must actually fire on a
// table that does not satisfy the recurrence. Perturbing one interior value
// is the cheapest way to be sure the check is not vacuous -- done through the
// detail function, since a Wigner3jMatrix that has been built is by
// construction one that passed.
TEST(ThreeJ, TheRecurrenceCheckIsNotVacuous) {
  constexpr auto l1 = 8, l2 = 8, l3 = 8, m1 = 0;
  auto row = std::vector<double>(2 * l2 + 2);
  const auto n = ThreeJDetails::SchultenGordonRow<double>(
      l1, l2, l3, m1, std::span<double>(row));
  ASSERT_GT(n, 4);

  EXPECT_TRUE(ThreeJDetails::RowSatisfiesRecurrence<double>(
      l1, l2, l3, m1, std::span<const double>(row.data(), n), n));

  auto broken = row;
  broken[static_cast<std::size_t>(n / 2)] *= 1.5;
  EXPECT_FALSE(ThreeJDetails::RowSatisfiesRecurrence<double>(
      l1, l2, l3, m1, std::span<const double>(broken.data(), n), n));
}

// An independent formula, which is the one check the algorithm cannot make
// about itself. Completeness holds by construction here -- Schulten-Gordon
// normalises each row by it -- and the recurrence residual is homogeneous, so
// it is blind to an overall scale or sign. Racah sees all of that.
//
// Restricted to entries whose Racah sum is short, because that is where the
// oracle is trustworthy: at l3 = l1 + l2 the sum has one term and cannot
// cancel at all, and it degrades as the sum lengthens. Comparing against a
// long Racah sum would measure Racah rather than the library.
TEST(ThreeJ, AgreesWithRacahWhereRacahIsExact) {
  for (const auto& t : std::vector<std::array<int, 3>>{{3, 4, 5},
                                                       {12, 12, 12},
                                                       {40, 40, 80},
                                                       {64, 64, 128},
                                                       {80, 80, 120},
                                                       {100, 100, 150},
                                                       {128, 128, 192},
                                                       {128, 128, 250},
                                                       {160, 160, 160},
                                                       {200, 200, 200}}) {
    const auto table = Wigner3jMatrix<double>(t[0], t[1], t[2]);
    auto worst = 0.0;
    auto checked = 0;
    for (auto m1 : table.M1Axis()) {
      for (auto m3 : table.M3Axis()) {
        const auto m2 = -(m1 + m3);
        if (std::abs(m2) > t[1]) continue;
        if (RacahSumLength(t[0], t[1], t[2], m1, m2) > 6) continue;
        worst = std::max(
            worst, std::abs(table(m1, m3) -
                            RacahSymbol<double>(t[0], t[1], t[2], m1, m2, m3)));
        ++checked;
      }
    }
    EXPECT_GT(checked, 0) << "(" << t[0] << "," << t[1] << "," << t[2]
                          << ") gave the oracle nothing to check";

    // The tolerance is a bound on *Racah's* accuracy, not the library's.
    // Racah exponentiates a logarithm of size O(l), so it loses bits in
    // proportion to the degree -- measured, 1.8e-11 at (128,128,192) and
    // 6.6e-11 at (200,200,200), against 1e-16 for the library's own
    // cyclic-permutation check on the same triples. So this test is much
    // blunter than the structural ones and is written to say so.
    const auto tolerance = 1e-12 * (t[0] + t[1] + t[2]);
    EXPECT_LT(worst, tolerance)
        << "(" << t[0] << "," << t[1] << "," << t[2] << ") worst " << worst;
  }
}

// And the oracle is worth having only if it is right where it is used, so its
// short-sum values are themselves pinned against the two closed forms.
TEST(ThreeJ, TheRacahOracleIsExactWhereItIsTrusted) {
  for (auto l1 : {1, 5, 16, 64, 128}) {
    for (auto l2 : {1, 7, 32}) {
      const auto l3 = l1 + l2;
      ASSERT_EQ(RacahSumLength(l1, l2, l3, l1, l2), 1);
      EXPECT_NEAR(RacahSymbol<double>(l1, l2, l3, l1, l2, -l3),
                  1 / std::sqrt(2.0 * l3 + 1), 1e-13)
          << "l1 = " << l1 << ", l2 = " << l2;
    }
  }
  for (auto j : {1, 5, 30}) {
    for (auto m = -j; m <= j; ++m) {
      EXPECT_NEAR(RacahSymbol<double>(j, j, 0, m, -m, 0),
                  ((j - m) % 2 == 0 ? 1.0 : -1.0) / std::sqrt(2.0 * j + 1),
                  1e-13)
          << "j = " << j << ", m = " << m;
    }
  }
}

// The coupling layout is a convention: a mirror in m1 and a phase. Applying
// the swap twice must give the table back, which is what says the two
// conventions cannot drift apart.
TEST(ThreeJ, TheCouplingLayoutIsAConventionAndIsItsOwnInverse) {
  constexpr auto l1 = 4, l2 = 5, l3 = 6;
  const auto plain = Wigner3jMatrix<double>(l1, l2, l3);

  auto coupling =
      std::vector<double>(static_cast<std::size_t>(2 * l1 + 1) * (2 * l3 + 1));
  FillCouplingMatrix(l1, l2, l3, coupling);

  // c(m, mp) = (-1)^m (l1 l2 l3; -m, m-mp, mp)
  const auto columns = 2 * l3 + 1;
  for (auto m = -l1; m <= l1; ++m) {
    for (auto mp = -l3; mp <= l3; ++mp) {
      const auto phase = (m % 2 == 0) ? 1.0 : -1.0;
      const auto want = phase * plain(-m, mp);
      const auto got =
          coupling[static_cast<std::size_t>(m + l1) * columns + (mp + l3)];
      EXPECT_NEAR(got, want, 1e-15) << "m = " << m << ", mp = " << mp;
      EXPECT_NEAR(plain.CouplingElement(m, mp), want, 1e-15);
    }
  }

  auto twice = coupling;
  ThreeJDetails::SwapCouplingConvention<double>(l1, l3,
                                                std::span<double>(twice));
  for (std::size_t i = 0; i < twice.size(); ++i) {
    EXPECT_NEAR(twice[i], plain.Data()[i], 1e-15) << "entry " << i;
  }
}

//--------------------------------------------------------------------------//
//                      Stretched triangles at high degree                   //
//--------------------------------------------------------------------------//

// A row of a near-stretched table is a single hump with no classically
// allowed region, and the rule that decides where the two halves of the
// recursion meet has to find the top of it. It once did not -- it watched the
// recurrence coefficient, which is a proxy for the values and a poor one
// here -- so one half ran downhill, which costs digits, and beyond a certain
// disparity the join overflowed and the table threw. More than half of all
// (l, 2l, l) tables up to l = 1000 could not be built, and in single precision
// tables came back wrong by a third of their own scale without complaint.
//
// The oracle is Racah's closed form, which is exact to rounding where its sum
// is short. A short sum *is* the stretched region, reached by whichever cyclic
// permutation of the columns puts the stretched degree third.
//
// GSHTRANS_TEST_THOROUGH=1 runs the full sweeps the fix was accepted on, which
// take a few seconds optimised and rather longer under a sanitiser.

namespace Stretched {

inline bool Thorough() {
  const auto* value = std::getenv("GSHTRANS_TEST_THOROUGH");
  return value != nullptr && std::string(value) != "0";
}

inline std::vector<int> Degrees() {
  if (Thorough()) {
    auto all = std::vector<int>{};
    for (auto l = 2; l <= 1000; l += 13) all.push_back(l);
    return all;
  }
  // Small, the first failures, the band where they were densest, and one well
  // beyond it.
  return {2, 33, 130, 260, 264, 275, 302, 520};
}

// The worst disagreement with Racah over the entries it can be trusted on, as
// a fraction of the table's largest entry.
template <typename T>
double WorstAgainstRacah(int l1, int l2, int l3) {
  const auto table = Wigner3jMatrix<T>(l1, l2, l3);
  auto scale = 0.0;
  for (auto v : table.Data()) scale = std::max<double>(scale, std::abs(v));

  auto worst = 0.0;
  for (auto m1 = -l1; m1 <= l1; m1 += std::max(1, l1 / 12)) {
    for (auto m3 = -l3; m3 <= l3; m3 += std::max(1, l3 / 12)) {
      const auto m2 = -m1 - m3;
      if (std::abs(m2) > l2) continue;
      // An even permutation of the columns leaves the symbol unchanged, so
      // use whichever arrangement Racah finds shortest.
      const std::array<std::array<int, 6>, 3> cyclic{
          {{l1, l2, l3, m1, m2, m3},
           {l2, l3, l1, m2, m3, m1},
           {l3, l1, l2, m3, m1, m2}}};
      for (const auto& a : cyclic) {
        // Two terms at most: at three the oracle's own cancellation shows.
        if (RacahSumLength(a[0], a[1], a[2], a[3], a[4]) > 2) continue;
        const auto exact =
            RacahSymbol<long double>(a[0], a[1], a[2], a[3], a[4], a[5]);
        worst = std::max<double>(
            worst, std::abs(static_cast<long double>(table(m1, m3)) - exact));
        break;
      }
    }
  }
  return worst / scale;
}

template <typename A, typename B>
double WorstBetween(int l1, int l2, int l3) {
  const auto a = Wigner3jMatrix<A>(l1, l2, l3);
  const auto b = Wigner3jMatrix<B>(l1, l2, l3);
  auto scale = 0.0;
  auto worst = 0.0;
  for (std::size_t i = 0; i < a.Data().size(); ++i) {
    scale = std::max<double>(scale, std::abs(b.Data()[i]));
    worst = std::max<double>(worst,
                             std::abs(static_cast<long double>(a.Data()[i]) -
                                      static_cast<long double>(b.Data()[i])));
  }
  return worst / scale;
}

inline std::vector<std::array<int, 3>> RandomTriples(int count, int lTop) {
  auto gen = std::mt19937(20260919u);
  auto triples = std::vector<std::array<int, 3>>{};
  for (auto i = 0; i < count; ++i) {
    const auto l1 = static_cast<int>(gen() % static_cast<unsigned>(lTop + 1));
    const auto l3 = static_cast<int>(gen() % static_cast<unsigned>(lTop + 1));
    const auto lowest = std::abs(l1 - l3);
    const auto l2 =
        lowest +
        static_cast<int>(gen() % static_cast<unsigned>(l1 + l3 - lowest + 1));
    triples.push_back({l1, l2, l3});
  }
  return triples;
}

}  // namespace Stretched

TEST(ThreeJ, StretchedTablesCanBeBuiltAndAgreeWithRacah) {
  for (auto l : Stretched::Degrees()) {
    for (const auto& t : std::vector<std::array<int, 3>>{{l, 2 * l, l},
                                                         {l, 2 * l - 1, l},
                                                         {l, 2 * l - 2, l},
                                                         {l, l, 2 * l}}) {
      auto worst = 0.0;
      ASSERT_NO_THROW(
          worst = Stretched::WorstAgainstRacah<double>(t[0], t[1], t[2]))
          << "(" << t[0] << "," << t[1] << "," << t[2] << ")";
      // The oracle's floor, not the library's: a logarithm of size O(l)
      // exponentiated in double precision.
      EXPECT_LT(worst, 1e-11)
          << "(" << t[0] << "," << t[1] << "," << t[2] << ")";
    }
  }
}

TEST(ThreeJ, StretchedInTheMiddleDegreeCanBeBuilt) {
  // l2 = l1 + l3 with the outer two different, which is what Wigner3jStack
  // reaches by default at the top of its range.
  const auto step = Stretched::Thorough() ? 25 : 125;
  for (auto l1 = 150; l1 <= 400; l1 += step) {
    for (auto l3 = 150; l3 <= 400; l3 += step) {
      auto worst = 0.0;
      ASSERT_NO_THROW(worst =
                          Stretched::WorstAgainstRacah<double>(l1, l1 + l3, l3))
          << "(" << l1 << "," << l1 + l3 << "," << l3 << ")";
      EXPECT_LT(worst, 1e-11)
          << "(" << l1 << "," << l1 + l3 << "," << l3 << ")";
    }
  }
  EXPECT_NO_THROW((Wigner3jStack<double>(40, 30)));
}

TEST(ThreeJ, DoublePrecisionHoldsAwayFromTheStretchedRegion) {
  // Against long double, which is the same algorithm and so no oracle for a
  // wrong answer, but a good one for lost digits. The bound is set by flat
  // rows -- l1 = 0 has 2 l + 1 equal entries -- where any three-term
  // recurrence accumulates about n^2 epsilon.
  const auto triples = Stretched::RandomTriples(
      Stretched::Thorough() ? 300 : 20, Stretched::Thorough() ? 500 : 350);
  for (const auto& t : triples) {
    auto worst = 0.0;
    ASSERT_NO_THROW(worst = (Stretched::WorstBetween<double, long double>(
                        t[0], t[1], t[2])))
        << "(" << t[0] << "," << t[1] << "," << t[2] << ")";
    EXPECT_LT(worst, 1e-10) << "(" << t[0] << "," << t[1] << "," << t[2] << ")";
  }
}

TEST(ThreeJ, SinglePrecisionIsRightAndNotMerelyReturned) {
  // Single precision used to come back wrong by up to six tenths of the
  // table's own scale, having passed the residual check, and threw from
  // l = 33 on stretched tables. The first three named here are from the
  // review. The flat rows are here for a different reason: even a correct
  // recursion loses n^2 epsilon along a row, which in single precision left
  // two digits at l = 450, so the rows are computed in double and a
  // single-precision table is the double one, rounded.
  auto triples = Stretched::RandomTriples(Stretched::Thorough() ? 300 : 20,
                                          Stretched::Thorough() ? 500 : 350);
  triples.push_back({74, 466, 392});
  triples.push_back({203, 480, 277});
  triples.push_back({40, 70, 30});
  triples.push_back({170, 340, 170});
  triples.push_back({0, 448, 448});
  triples.push_back({2, 300, 301});
  for (const auto& t : triples) {
    auto worst = 0.0;
    ASSERT_NO_THROW(
        worst = (Stretched::WorstBetween<float, double>(t[0], t[1], t[2])))
        << "(" << t[0] << "," << t[1] << "," << t[2] << ")";
    // Rounding to single precision and nothing more.
    EXPECT_LT(worst, 1e-6) << "(" << t[0] << "," << t[1] << "," << t[2] << ")";
  }
}

}  // namespace
