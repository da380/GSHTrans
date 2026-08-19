#ifndef GSH_TRANS_UTILITY_GUARD_H
#define GSH_TRANS_UTILITY_GUARD_H

#include <cstddef>
#include <initializer_list>

namespace GSHTrans {

// Returns (-1)^m for integer m.
template <typename Scalar = std::ptrdiff_t>
constexpr auto MinusOneToPower(std::ptrdiff_t m) -> Scalar {
  return m % 2 ? -1 : 1;
}

// True if n is a length for which FFTW has hard-coded codelets, namely
// n = 2^a 3^b 5^c 7^d 11^e 13^f with e + f either zero or one. Other lengths
// fall back to a general algorithm, which for a large prime factor is much
// slower: 514 = 2 * 257 is the case that matters here.
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

// The smallest such length that is at least n.
constexpr auto FastFFTSize(std::ptrdiff_t n) -> std::ptrdiff_t {
  if (n < 1) return 1;
  while (!IsFastFFTSize(n)) n++;
  return n;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_UTILITY_GUARD_H