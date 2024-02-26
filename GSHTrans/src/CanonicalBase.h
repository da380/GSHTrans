#ifndef GSH_TRANS_CANONICAL_BASE_GUARD_H
#define GSH_TRANS_CANONICAL_BASE_GUARD_H

#include <ranges>

#include "Concepts.h"
#include "GridBase.h"

namespace GSHTrans {

template <typename Derived>
class GridInterface {
  using Int = std::ptrdiff_t;

 public:
  auto NumberOfCoLatitudes() const {
    return _Derived().Grid().NumberOfCoLatitudes();
  }
  auto NumberOfLongitudes() const {
    return _Derived().Grid().NumberOfLongitudes();
  }
  auto CoLatitudes() const { return _Derived().Grid().CoLatitudes(); }
  auto Longitudes() const { return _Derived().Grid().Longitudes(); }
  auto Points() const { return _Derived().Grid().Points(); }

  auto CoLatitudeIndices() const {
    return _Derived().Grid().CoLatitudeIndices();
  }
  auto LongitudeIndices() const { return _Derived().Grid().LongitudeIndices(); }
  auto PointIndices() const { return _Derived().Grid().PointIndices(); }

 private:
  auto& _Derived() const { return static_cast<const Derived&>(*this); }
  auto& _Derived() { return static_cast<Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_BASE_GUARD_H