#ifndef RACAH_REFERENCE_GUARD
#define RACAH_REFERENCE_GUARD

// Racah's closed form for a Wigner 3-j symbol, as a test oracle.
//
// This is a *second implementation*, and that is the whole of why it is here.
// The library computes its symbols by Schulten-Gordon recursion, and the
// checks that algorithm can make about itself are limited: the completeness
// relation is what it normalises by, so it holds by construction, and the
// recurrence residual is homogeneous, so it cannot see an overall scale or
// sign. An independent formula can see all of that.
//
// It lives in the test tree rather than the library because it is not an
// alternative production path -- it was measured against the recursion and
// loses badly over most of the triangle space (docs/3j-plan.md T2). What it
// is good at is the region where its sum is *short*:
//
//     length = min(l1+l2-l3, l1-m1, l2+m2) - max(0, l2-l3-m1, l1-l3+m2) + 1
//
// At l3 = l1 + l2 that is one term, so there is no cancellation at all and
// the result is exact to rounding. As the sum lengthens the alternating terms
// cancel and it degrades, catastrophically for fat triangles at high degree.
//
// **So a test must restrict itself to short sums**, and RacahSumLength is
// provided so that it can. Comparing against a long Racah sum measures Racah,
// not the library.
//
// The formula is evaluated in log space -- the factorials are formed as
// lgamma and the prefactor is folded into each term -- so nothing overflows
// however large the degrees. The residual drift from exponentiating a
// logarithm grows with degree, and is about 1e-12 at l = 500.

#include <GSHTrans/Core>
#include <algorithm>
#include <cmath>
#include <concepts>

// The number of terms in Racah's alternating sum for one symbol. A test
// should skip entries where this is large: the oracle is untrustworthy there
// and the library is not.
inline int RacahSumLength(int l1, int l2, int l3, int m1, int m2) {
  const auto lowest = std::max({0, l2 - l3 - m1, l1 - l3 + m2});
  const auto highest = std::min({l1 + l2 - l3, l1 - m1, l2 + m2});
  return highest - lowest + 1;
}

template <std::floating_point Real>
Real RacahSymbol(int l1, int l2, int l3, int m1, int m2, int m3) {
  using namespace GSHTrans;

  if (m1 + m2 + m3 != 0) return Real{0};
  if (std::abs(m1) > l1 or std::abs(m2) > l2 or std::abs(m3) > l3) {
    return Real{0};
  }
  if (not SatisfiesTriangle(l1, l2, l3)) return Real{0};

  const auto logFactorial = [](int n) {
    return static_cast<Real>(std::lgamma(static_cast<double>(n) + 1.0));
  };

  const auto logDelta =
      logFactorial(l1 + l2 - l3) + logFactorial(l1 - l2 + l3) +
      logFactorial(-l1 + l2 + l3) - logFactorial(l1 + l2 + l3 + 1);
  const auto logNumerator = logFactorial(l1 + m1) + logFactorial(l1 - m1) +
                            logFactorial(l2 + m2) + logFactorial(l2 - m2) +
                            logFactorial(l3 + m3) + logFactorial(l3 - m3);
  const auto logPrefactor = (logDelta + logNumerator) / 2;

  const auto lowest = std::max({0, l2 - l3 - m1, l1 - l3 + m2});
  const auto highest = std::min({l1 + l2 - l3, l1 - m1, l2 + m2});

  auto sum = Real{0};
  for (auto k = lowest; k <= highest; ++k) {
    const auto logDenominator =
        logFactorial(k) + logFactorial(l1 + l2 - l3 - k) +
        logFactorial(l1 - m1 - k) + logFactorial(l2 + m2 - k) +
        logFactorial(l3 - l2 + m1 + k) + logFactorial(l3 - l1 - m2 + k);
    const auto term = std::exp(logPrefactor - logDenominator);
    sum += (k % 2 == 0) ? term : -term;
  }
  return ((l1 - l2 - m3) % 2 == 0) ? sum : -sum;
}

#endif  // RACAH_REFERENCE_GUARD
