#ifndef GSH_TRANS_FIELDS_BASE_GUARD_H
#define GSH_TRANS_FIELDS_BASE_GUARD_H

#include <FFTWpp/All>
#include <algorithm>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>
#include <numeric>
#include <ranges>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

// Base class for canonical components.
template <typename Derived>
class CanonicalComponentBase {
  using Int = std::ptrdiff_t;

 public:
  auto& GridReference() const { return GetDerived()._GridReference(); }

  auto MaxDegree() const { return GridReference().MaxDegree(); }
  auto MaxUpperIndex() const { return GridReference().MaxUpperIndex(); }

  auto CoLatitudes() const { return GridReference().CoLatitudes(); }
  auto NumberOfCoLatitudes() const {
    return GridReference().NumberOfCoLatitudes();
  }
  auto CoLatitudeIndices() const { return GridReference().CoLatitudeIndices(); }

  auto LongitudeSpacing() const { return GridReference().LongitudeSpacing(); }
  auto LongitudeIndices() const { return GridReference().LongitudeIndices(); }
  auto Longitudes() const { return GridReference().Lognitudes(); }
  auto NumberOfLongitudes() const {
    return GridReference().NumberOfLongitudes();
  }

  auto View() { return GetDerived()._View(); }

  auto begin() { return View().begin(); }
  auto end() { return View().end(); }

  auto size() { return View().size(); }

  auto operator[](Int i) { return View()[i]; }

  auto operator()(Int iTheta, Int iPhi) {
    auto i = NumberOfLongitudes() * iTheta + iPhi;
    return operator[](i);
  }

  auto Integrate() { return GridReference().Integrate(View()); }

 private:
  auto& GetDerived() { return static_cast<Derived&>(*this); }
  const auto& GetDerived() const { return static_cast<const Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELDS_BASE_GUARD_H
