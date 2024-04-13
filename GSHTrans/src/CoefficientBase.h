#ifndef GSH_TRANS_COEFFICIENT_BASE_GUARD_H
#define GSH_TRANS_COEFFICIENT_BASE_GUARD_H

#include <algorithm>
#include <cassert>
#include <cmath>

namespace GHSTrans {

template <typename _Derived>
class CoefficientBase {
  using Int = std::ptrdiff_t;

 public:
  auto GetGrid() const { return Derived().GetGrid(); }

  auto MaxDegree() const { return GetGrid().MaxDegree(); }

  void CheckDegree(Int l, Int n) const {
    assert(l >= std::abs(n) && l <= MaxDegree());
  }

  void CheckUpperIndex(Int n) const {
    assert(std::ranges::contains(GetGrid().UpperIndices(), n));
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GHSTrans

#endif  // GSH_TRANS_COEFFICIENT_BASE_GUARD_H