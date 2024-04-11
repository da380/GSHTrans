#ifndef GSH_TRANS_FIELD_BASE_GUARD_H
#define GSH_TRANS_FIELD_BASE_GUARD_H

#include <cassert>

namespace GSHTrans {

template <typename _Derived>
class FieldBase {
  using Int = std::ptrdiff_t;

 public:
  // Methods related to the grid.
  auto GetGrid() const { return Derived().GetGrid(); }

  auto NumberOfCoLatitudes() const { return GetGrid().NumberOfCoLatitudes(); }
  auto CoLatitudes() const { return GetGrid().CoLatitudes(); }
  auto CoLatitudeIndices() const { return GetGrid().CoLatitudeIndices(); }

  auto NumberOfLongitudes() const { return GetGrid().NumberOfLongitudes(); }
  auto Longitudes() const { return GetGrid().Longitudes(); }
  auto LongitudeIndices() const { return GetGrid().LongitudeIndices(); }

  void CheckPointIndices(Int iTheta, Int iPhi) const {
    assert(iTheta >= 0 && iTheta <= this->NumberOfCoLatitudes());
    assert(iPhi >= 0 && iPhi <= this->NumberOfLongitudes());
  }

  auto Points() const { return GetGrid().Points(); }
  auto PointIndices() const { return GetGrid().PointIndices(); }

  auto Weights() const { return GetGrid().Weights(); }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELD_BASE_GUARD_H