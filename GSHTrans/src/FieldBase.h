#ifndef GSH_TRANS_FIELD_BASE_GUARD_H
#define GSH_TRANS_FIELD_BASE_GUARD_H

namespace GSHTrans {

template <typename _Derived>
class FieldBase {
  using Int = std::ptrdiff_t;

 public:
  // Methods related to the grid.
  auto GetGrid() const { return Derived().GetGrid(); }
  auto NumberOfCoLatitudes() const { return GetGrid().NumberOfCoLatitudes(); }
  auto NumberOfLongitudes() const { return GetGrid().NumberOfLongitudes(); }
  auto CoLatitudes() const { return GetGrid().CoLatitudes(); }
  auto Longitudes() const { return GetGrid().Longitudes(); }
  auto Points() const { return GetGrid().Points(); }
  auto PointIndices() const { return GetGrid().PointIndices(); }

  // Methods related to the data.
  auto Index(Int iTheta, int iPhi) const {
    return iTheta * NumberOfLongitudes() + iPhi;
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELD_BASE_GUARD_H