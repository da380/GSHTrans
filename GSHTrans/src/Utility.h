#ifndef GSH_TRANS_UTILITY_GUARD_H
#define GSH_TRANS_UTILITY_GUARD_H

/**
 * @file Utility.h
 * @brief Small numerical helpers used across the library.
 */

#include <cstddef>
#include <initializer_list>

namespace GSHTrans {

/**
 * @brief Returns @f$(-1)^m@f$.
 * @tparam Scalar The type of the result.
 * @param m The exponent; only its parity is used.
 */
template <typename Scalar = std::ptrdiff_t>
constexpr auto MinusOneToPower(std::ptrdiff_t m) -> Scalar {
  return m % 2 ? -1 : 1;
}

/**
 * @brief Whether FFTW has hard-coded codelets for a transform of length @p n.
 *
 * @details The lengths it does are @f$n = 2^a 3^b 5^c 7^d 11^e 13^f@f$ with
 * @f$e + f@f$ either zero or one. Any other length falls back to a general
 * algorithm, which for a large prime factor is much slower — and the grid's
 * longitude count runs into exactly that, since the smallest length resolving
 * orders @f$|m| \le 256@f$ is 513, whose next even neighbour 514 is
 * @f$2 \times 257@f$ with 257 prime.
 *
 * @param n The length to test.
 * @return True if @p n has only the small prime factors above.
 */
constexpr auto IsFastFFTSize(std::ptrdiff_t n) -> bool {
  if (n < 1) return false;
  for (auto p : {2, 3, 5, 7}) {
    while (n % p == 0) n /= p;
  }
  auto large = 0;
  for (auto p : {11, 13}) {
    if (n % p == 0) {
      n /= p;
      large++;
    }
  }
  return n == 1 && large <= 1;
}

/**
 * @brief The smallest length at least @p n for which IsFastFFTSize holds.
 * @param n The lower bound.
 */
constexpr auto FastFFTSize(std::ptrdiff_t n) -> std::ptrdiff_t {
  if (n < 1) return 1;
  while (!IsFastFFTSize(n)) n++;
  return n;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_UTILITY_GUARD_H
