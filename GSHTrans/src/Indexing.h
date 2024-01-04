#ifndef GSH_TRANS_INDEXING_GUARD_H
#define GSH_TRANS_INDEXING_GUARD_H

#include <ranges>

#include "Concepts.h"

namespace GSHTrans {

// Return view to the degrees for a given maximum value.
auto Degrees(int lMax) { return std::ranges::views::iota(0, lMax + 1); }

// View to the orders for a given degree
template <TransformType Type = C2C>
auto Orders(int l) {
  if constexpr (std::same_as<Type, C2C>) {
    return std::ranges::views::iota(-l, l + 1);
  } else {
    return std::ranges::views::iota(0, l + 1);
  }
}

// Return view to the upper indices for given nMax.
template <TransformType Type = C2C>
auto UpperIndices(int nMax) {
  if constexpr (std::same_as<Type, C2C>) {
    return std::ranges::views::iota(-nMax, nMax + 1);
  } else {
    return std::ranges::views::iota(0, nMax + 1);
  }
}

// Range for the spherical harmonic indices given lMax.
template <TransformType Type = C2C>
class SphericalHarmonicIndices {
 public:
  SphericalHarmonicIndices() = delete;
  explicit SphericalHarmonicIndices(int lMax) : _lMax{lMax}, _l{0}, _m{0} {}

  const auto& begin() const { return *this; }
  const auto& end() const { return *this; }

  auto operator!=(const SphericalHarmonicIndices&) const { return _l <= _lMax; }

  void operator++() {
    if (_m == _l) {
      _l++;
      if constexpr (std::same_as<Type, C2C>) {
        _m = -_l;
      } else {
        _m = 0;
      }
    } else {
      _m++;
    }
  }

  auto operator*() const { return std::pair(_l, _m); }

 private:
  int _lMax;
  int _l;
  int _m;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_INDEXING_GUARD_H
