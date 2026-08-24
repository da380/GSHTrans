#include <gtest/gtest.h>

#include <GSHTrans/Core>

#include <cmath>
#include <array>
#include <cstddef>
#include <vector>

// T1 of docs/3j-plan.md: the test family, before any algorithm changes.
//
// `3j.h` had no coverage at all although it is public API, and the reason
// this file can be written without a reference implementation is the
// completeness relation:
//
//     sum over the (m1, m3) plane of the squared symbols  =  1
//
// exactly, for any triple satisfying the triangle rule. That is one loop over
// a table that has just been built, it needs nothing to compare against, and
// it is what the rest of the plan is built on.
//
// Half of what is asserted here is where the identity **fails**, which is as
// important as where it holds: the existing recursion is run in its unstable
// direction near stretched triangles, and pinning the boundary is what turns
// an unknown into a known one. T4 will invert those assertions when Racah's
// closed form covers the region.

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
    EXPECT_NEAR(Completeness<double>(l, l, l), 1.0, 1e-12) << "(l,l,l), l = " << l;
  }
}

TEST(ThreeJ, CompletenessHoldsAwayFromStretched) {
  // A spread of shapes, none of them close to l3 = l1 + l2.
  const auto triples = std::vector<std::array<int, 3>>{
      {4, 4, 2},  {6, 4, 4},  {10, 7, 5},  {12, 12, 6},
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
// 30, which sharpens what 3j-plan.md section 1 records. Two assertions
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

// [J1]: the boundary of section 1 of 3j-plan.md, and it is now a *refusal*
// rather than a wrong answer. Before the self-check these triples returned
// finite, plausibly-shaped tables holding numbers of order 1e112.
//
// These are the assertions T4 inverts: when Racah's closed form covers the
// stretched region they become successes, and the commit that changes them is
// the visible record of what was fixed.
TEST(ThreeJ, RefusesStretchedTrianglesOfHighDegree) {
  for (auto l : {35, 40, 50}) {
    EXPECT_THROW(Wigner3jMatrix<double>(l, l, 2 * l), std::runtime_error)
        << "(l, l, 2l) at l = " << l
        << ": the present recursion is badly wrong here and must say so. If "
           "this now fails, the algorithm has been fixed and the test should "
           "be inverted (3j-plan.md T4)";
  }
}

TEST(ThreeJ, RefusesNearStretchedTrianglesOfHighDegree) {
  for (auto l : {100, 128}) {
    EXPECT_THROW(Wigner3jMatrix<double>(l, l, 3 * l / 2), std::runtime_error)
        << "(l, l, 3l/2) at l = " << l;
  }
}

// And the check does not fire on anything good, which is the half that would
// otherwise make it a nuisance. Every triple the tests above assert is
// accurate must construct without complaint.
TEST(ThreeJ, AcceptsEveryTripleItShould) {
  for (auto l : {1, 2, 3, 5, 8, 13, 20, 32, 64}) {
    EXPECT_NO_THROW(Wigner3jMatrix<double>(l, l, l)) << "(l,l,l), l = " << l;
  }
  for (auto l : {1, 2, 4, 8, 12, 16, 20}) {
    EXPECT_NO_THROW(Wigner3jMatrix<double>(l, l, 2 * l))
        << "(l, l, 2l), l = " << l;
  }
  for (auto l : {10, 20, 30, 40}) {
    EXPECT_NO_THROW(Wigner3jMatrix<double>(l, l, 3 * l / 2))
        << "(l, l, 3l/2), l = " << l;
  }
  // A triple that is not a triangle has an identically zero table, so its
  // completeness sum is zero and correctly so: the check must skip it rather
  // than refuse it.
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
// to be wrong in the same way the code is. It is also the case T4's Racah
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

}  // namespace
