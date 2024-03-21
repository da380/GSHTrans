#ifndef GSH_TRANS_GRID_GUARD_H
#define GSH_TRANS_GRID_GUARD_H

#include <concepts>
#include <ranges>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

template <typename Derived>
class GridBase {
  using Int = std::ptrdiff_t;

 public:
  auto MinUpperIndex() const {
    using NRange = Derived::NRange_type;
    if constexpr (std::same_as<NRange, All>) {
      return -_Derived().MaxUpperIndex();
    }
    if constexpr (std::same_as<NRange, NonNegative>) {
      return Int{0};
    }
    if constexpr (std::same_as<NRange, Single>) {
      return _Derived().MaxUpperIndex();
    }
  }

  auto UpperIndices() const {
    return std::ranges::views::iota(MinUpperIndex(),
                                    _Derived().MaxUpperIndex() + 1);
  }

  auto NumberOfCoLatitudes() const { return _Derived().CoLatitudes().size(); }
  auto CoLatitudeIndices() const {
    return std::ranges::views::iota(std::size_t{0}, NumberOfCoLatitudes());
  }

  auto NumberOfLongitudes() const { return _Derived().Longitudes().size(); }
  auto LongitudeIndices() const {
    return std::ranges::views::iota(std::size_t{0}, NumberOfLongitudes());
  }

  auto Points() const {
    return std::ranges::views::cartesian_product(_Derived().CoLatitudes(),
                                                 _Derived().Longitudes());
  }

  auto PointIndices() const {
    return std::ranges::views::cartesian_product(_Derived().CoLatitudeIndices(),
                                                 _Derived().LongitudeIndices());
  }

  auto Weights() const {
    return std::ranges::views::cartesian_product(
               _Derived().CoLatitudeWeights(), _Derived().LongitudeWeights()) |
           std::ranges::views::transform(
               [](auto pair) { return std::get<0>(pair) * std::get<1>(pair); });
  }

  template <typename Function>
  auto InterpolateFunction(Function f) {
    return Points() | std::ranges::views::transform([f](auto pair) {
             auto [theta, phi] = pair;
             return f(theta, phi);
           });
  }

  auto ComponentSize() const {
    return NumberOfCoLatitudes() * NumberOfLongitudes();
  }

  auto RealCoefficientSize(Int lMax, Int n) const {
    return GSHIndices<NonNegative>(lMax, lMax, n).size();
  }

  auto ComplexCoefficientSize(Int lMax, Int n) const {
    return GSHIndices<All>(lMax, lMax, n).size();
  }

  auto RealCoefficientSize(Int n) const {
    auto lMax = _Derived().MaxDegree();
    return GSHIndices<NonNegative>(lMax, lMax, n).size();
  }

  auto ComplexCoefficientSize(Int n) const {
    auto lMax = _Derived().MaxDegree();
    return GSHIndices<All>(lMax, lMax, n).size();
  }

  template <ComplexFloatingPointRange Range,
            typename Distribution =
                decltype(std::normal_distribution<
                         RemoveComplex<std::ranges::range_value_t<Range>>>())>
  void RandomComplexCoefficient(
      Int lMax, Int n, Range& range,
      Distribution dist = std::normal_distribution<
          RemoveComplex<std::ranges::range_value_t<Range>>>()) const {
    using Complex = std::ranges::range_value_t<Range>;
    using Real = RemoveComplex<Complex>;
    assert(range.size() == ComplexCoefficientSize(lMax, n));
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::ranges::generate(range, [&gen, &dist]() {
      return Complex{dist(gen), dist(gen)};
    });

    if (lMax == _Derived().MaxDegree()) {
      auto i = GSHIndices<All>(lMax, lMax, n).Index(lMax, lMax);
      range[i] = 0;
    }
  }

  template <ComplexFloatingPointRange Range,
            typename Distribution =
                decltype(std::normal_distribution<
                         RemoveComplex<std::ranges::range_value_t<Range>>>())>
  void RandomRealCoefficient(
      Int lMax, Int n, Range& range,
      Distribution dist = std::normal_distribution<
          RemoveComplex<std::ranges::range_value_t<Range>>>()) const {
    using Complex = std::ranges::range_value_t<Range>;
    using Real = RemoveComplex<Complex>;
    assert(range.size() == RealCoefficientSize(lMax, n));
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::ranges::generate(range, [&gen, &dist]() {
      return Complex{dist(gen), dist(gen)};
    });

    auto indices = GSHIndices<NonNegative>(lMax, lMax, n);
    for (auto l : indices.Degrees()) {
      auto i = indices.Index(l, 0);
      range[i].imag(0);
    }
    if (lMax == _Derived().MaxDegree()) {
      auto i = indices.Index(lMax, lMax);
      range[i].imag(0);
    }
  }

 private:
  auto& _Derived() const { return static_cast<const Derived&>(*this); }
  auto& _Derived() { return static_cast<Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GRID_GUARD_H