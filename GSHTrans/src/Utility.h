#ifndef GSH_TRANS_UTILITY_GUARD_H
#define GSH_TRANS_UTILITY_GUARD_H

#include <limits>

namespace GSHTrans {

// Returns (-1)^m for integer m.
template <std::integral Int>
constexpr auto MinusOneToPower(Int m) {
  return m % 2 ? -1 : 1;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_UTILITY_GUARD_H