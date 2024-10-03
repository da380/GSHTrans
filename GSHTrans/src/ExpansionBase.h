#ifndef GSH_TRANS_EXPANSION_BASE_GUARD_H
#define GSH_TRANS_EXPANSION_BASE_GUARD_H

#include <algorithm>
#include <cassert>
#include <cmath>

#include "Traits.h"

namespace GSHTrans {

template <typename _Derived>
class ExpansionBase {
 public:
  using Int = std::ptrdiff_t;

  auto& Grid() const { return Derived().Grid(); }
  auto MaxDegree() const { return Derived().Grid().MaxDegree(); }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_EXPANSION_BASE_GUARD_H