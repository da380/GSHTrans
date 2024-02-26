#ifndef GSH_TRANS_CANONICAL_BASE_GUARD_H
#define GSH_TRANS_CANONICAL_BASE_GUARD_H

#include <ranges>

#include "Concepts.h"
#include "Grid.h"

namespace GSHTrans {

template <typename Derived>
class CanonicalBase
    : public std::ranges::view_interface<CanonicalBase<Derived>> {
  using Int = std::ptrdiff_t;

 public:
  // Data access functions.
  auto View() const { return _Derived()._View(); }
  auto View() { return _Derived()._View(); }

  auto begin() { return _Derived()._begin(); }
  auto end() { return _Derived()._end(); }

  // Return grid information,
  auto GridPointer() const { return _Derived()._grid; }

  auto NumberOfCoLatitudes() const {
    return GridPointer()->NumberOfCoLatitudes();
  }
  auto NumberOfLongitudes() const {
    return GridPointer()->NumberOfLongitudes();
  }
  auto CoLatitudes() const { return GridPointer()->CoLatitudes(); }
  auto Longitudes() const { return GridPointer()->Longitudes(); }
  auto Points() const { return GridPointer()->Points(); }

 private:
  auto& _Derived() const { return static_cast<const Derived&>(*this); }
  auto& _Derived() { return static_cast<Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_BASE_GUARD_H