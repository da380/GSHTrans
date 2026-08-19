#ifndef GSH_TRANS_GRID_GUARD_H
#define GSH_TRANS_GRID_GUARD_H

#include <concepts>
#include <cstddef>
#include <ranges>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

template <typename _Derived>
class GridBase {
  using Int = std::ptrdiff_t;

 public:
  auto MinUpperIndex() const {
    using NRange = _Derived::NRange;
    if constexpr (std::same_as<NRange, All>) {
      return -Derived().MaxUpperIndex();
    }
    if constexpr (std::same_as<NRange, NonNegative>) {
      return Int{0};
    }
    if constexpr (std::same_as<NRange, Single>) {
      return Derived().MaxUpperIndex();
    }
  }

  auto UpperIndices() const {
    return std::ranges::views::iota(MinUpperIndex(),
                                    Derived().MaxUpperIndex() + 1);
  }

  auto NumberOfCoLatitudes() const { return Derived().CoLatitudes().size(); }
  auto CoLatitudeIndices() const {
    return std::ranges::views::iota(Int{0},
                                    static_cast<Int>(NumberOfCoLatitudes()));
  }

  auto NumberOfLongitudes() const { return Derived().Longitudes().size(); }
  auto LongitudeIndices() const {
    return std::ranges::views::iota(Int{0},
                                    static_cast<Int>(NumberOfLongitudes()));
  }

  auto Points() const {
    return std::ranges::views::cartesian_product(Derived().CoLatitudes(),
                                                 Derived().Longitudes());
  }

  auto PointIndices() const {
    return std::ranges::views::cartesian_product(Derived().CoLatitudeIndices(),
                                                 Derived().LongitudeIndices());
  }

  auto CoLatitudeWeights() const { return Derived().CoLatitudeWeights(); }
  auto LongitudeWeights() const { return Derived().LongitudeWeights(); }

  auto Weights() const {
    return std::ranges::views::cartesian_product(CoLatitudeWeights(),
                                                 LongitudeWeights()) |
           std::ranges::views::transform(
               [](auto pair) { return std::get<0>(pair) * std::get<1>(pair); });
  }

  // Const: a value-semantic grid should be as usable through a const handle
  // as through a mutable one, and this reads nothing but the point set.
  template <typename Function>
  auto ProjectFunction(Function f) const {
    return Points() | std::ranges::views::transform([f](auto pair) {
             auto [theta, phi] = pair;
             return f(theta, phi);
           });
  }

  auto FieldSize() const {
    return NumberOfCoLatitudes() * NumberOfLongitudes();
  }

  // Number of coefficients of a complex-valued field of degree lMax at upper
  // index n. This is the full (all orders) storage.
  auto CoefficientSize(Int lMax, Int n) const {
    return GSHIndices<All>(lMax, lMax, n).Size();
  }

  auto CoefficientSize(Int n) const {
    return CoefficientSize(Derived().MaxDegree(), n);
  }

  // Number of stored coefficients of a real-valued field of degree lMax,
  // which uses the reduced m >= 0 storage. There is no upper-index argument
  // because real-valued fields exist only at upper index zero (core-plan.md
  // step A); the reduced storage is a statement about n = 0 alone.
  auto RealCoefficientSize(Int lMax) const {
    return GSHIndices<NonNegative>(lMax, lMax, 0).Size();
  }

  auto RealCoefficientSize() const {
    return RealCoefficientSize(Derived().MaxDegree());
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GRID_GUARD_H