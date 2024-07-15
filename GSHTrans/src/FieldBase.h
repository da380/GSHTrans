#ifndef GSH_TRANS_FIELD_BASE_GUARD_H
#define GSH_TRANS_FIELD_BASE_GUARD_H

#include <cassert>

#include "Traits.h"

namespace GSHTrans {

template <typename Derived>
class FieldBase {
 public:
  using Int = std::ptrdiff_t;

  auto GetGrid() const { return GetDerived().GetGrid(); }

  auto NumberOfCoLatitudes() const { return GetGrid().NumberOfCoLatitudes(); }
  auto CoLatitudes() const { return GetGrid().CoLatitudes(); }
  auto CoLatitudeIndices() const { return GetGrid().CoLatitudeIndices(); }

  auto NumberOfLongitudes() const { return GetGrid().NumberOfLongitudes(); }
  auto Longitudes() const { return GetGrid().Longitudes(); }
  auto LongitudeIndices() const { return GetGrid().LongitudeIndices(); }

  auto FieldSize() const {
    return NumberOfCoLatitudes() * NumberOfLongitudes();
  }

  constexpr void CheckPointIndices(Int iTheta, Int iPhi) const {
    assert(iTheta >= 0 && iTheta <= this->NumberOfCoLatitudes());
    assert(iPhi >= 0 && iPhi <= this->NumberOfLongitudes());
  }

  auto Points() const { return GetGrid().Points(); }
  auto PointIndices() const { return GetGrid().PointIndices(); }
  auto EnumeratePointIndices() const {
    return std::ranges::views::enumerate(PointIndices());
  }

  auto Weights() const { return GetGrid().Weights(); }

 private:
  auto& GetDerived() const { return static_cast<const Derived&>(*this); }
  auto& GetDerived() { return static_cast<Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELD_BASE_GUARD_H