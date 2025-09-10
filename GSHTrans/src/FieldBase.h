#ifndef GSH_TRANS_FIELD_BASE_GUARD_H
#define GSH_TRANS_FIELD_BASE_GUARD_H

#include <cassert>
#include <ranges>

#include "Traits.h"

namespace GSHTrans {

template <typename _Derived>
class FieldBase {
 public:
  using Int = std::ptrdiff_t;

  auto& Grid() const { return Derived().Grid(); }

  auto NumberOfCoLatitudes() const { return Grid().NumberOfCoLatitudes(); }
  auto CoLatitudes() const { return Grid().CoLatitudes(); }
  auto CoLatitudeIndices() const { return Grid().CoLatitudeIndices(); }

  auto NumberOfLongitudes() const { return Grid().NumberOfLongitudes(); }
  auto Longitudes() const { return Grid().Longitudes(); }
  auto LongitudeIndices() const { return Grid().LongitudeIndices(); }

  constexpr void CheckPointIndices(Int iTheta, Int iPhi) const {
    assert(iTheta >= 0 && iTheta <= NumberOfCoLatitudes());
    assert(iPhi >= 0 && iPhi <= NumberOfLongitudes());
  }

  auto Points() const { return Grid().Points(); }
  auto PointIndices() const { return Grid().PointIndices(); }
  auto EnumeratePointIndices() const {
    return std::ranges::views::enumerate(PointIndices());
  }

  auto Weights() const { return Grid().Weights(); }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELD_BASE_GUARD_H