#ifndef GSH_TRANS_UTILITY_GUARD_H
#define GSH_TRANS_UTILITY_GUARD_H

namespace GSHTrans {

// Returns (-1)^m for integer m.
template <typename Scalar = std::ptrdiff_t>
constexpr auto MinusOneToPower(std::ptrdiff_t m) -> Scalar {
  return m % 2 ? -1 : 1;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_UTILITY_GUARD_H