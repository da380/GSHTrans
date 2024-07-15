#ifndef GSH_TRANS_COEFFICIENT_BASE_GUARD_H
#define GSH_TRANS_COEFFICIENT_BASE_GUARD_H

#include <algorithm>
#include <cassert>
#include <cmath>

#include "Traits.h"

namespace GSHTrans {

template <typename Derived>
class CoefficientBase {
 public:
  using Int = std::ptrdiff_t;

  auto GetGrid() const { return GetDerived().GetGrid(); }
  auto MaxDegree() const { return GetDerived().GetGrid().MaxDegree(); }

 private:
  auto& GetDerived() const { return static_cast<const Derived&>(*this); }
  auto& GetDerived() { return static_cast<Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_COEFFICIENT_BASE_GUARD_H