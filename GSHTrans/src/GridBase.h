#ifndef GSH_TRANS_GRID_GUARD_H
#define GSH_TRANS_GRID_GUARD_H

#include <concepts>
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

  template <typename Function>
  auto ProjectFunction(Function f) {
    return Points() | std::ranges::views::transform([f](auto pair) {
             auto [theta, phi] = pair;
             return f(theta, phi);
           });
  }

  auto FieldSize() const {
    return NumberOfCoLatitudes() * NumberOfLongitudes();
  }

  auto RealCoefficientSize(Int lMax, Int n) const {
    return GSHIndices<NonNegative>(lMax, lMax, n).Size();
  }

  auto ComplexCoefficientSize(Int lMax, Int n) const {
    return GSHIndices<All>(lMax, lMax, n).Size();
  }

  auto RealCoefficientSize(Int n) const {
    auto lMax = Derived().MaxDegree();
    return GSHIndices<NonNegative>(lMax, lMax, n).Size();
  }

  auto ComplexCoefficientSize(Int n) const {
    auto lMax = Derived().MaxDegree();
    return GSHIndices<All>(lMax, lMax, n).Size();
  }

  auto CoefficientSize(Int lMax, Int n) const {
    return GSHIndices<All>(lMax, lMax, n).Size();
  }

  auto CoefficientSize(Int n) const {
    auto lMax = Derived().MaxDegree();
    return GSHIndices<All>(lMax, lMax, n).Size();
  }

  auto CoefficientSizeNonNegative(Int lMax, Int n) const {
    return GSHIndices<NonNegative>(lMax, lMax, n).Size();
  }

  auto CoefficientSizeNonNegative(Int n) const {
    auto lMax = Derived().MaxDegree();
    return GSHIndices<NonNegative>(lMax, lMax, n).Size();
  }

  template <std::ranges::range Range,
            typename Distribution =
                decltype(std::normal_distribution<
                         RemoveComplex<std::ranges::range_value_t<Range>>>())>
  requires requires() {
    requires ComplexFloatingPoint<std::ranges::range_value_t<Range>>;
    requires std::ranges::output_range<Range,
                                       std::ranges::range_value_t<Range>>;
  }
  void RandomComplexCoefficient(
      Int lMax, Int n, Range& range,
      Distribution dist = std::normal_distribution<
          RemoveComplex<std::ranges::range_value_t<Range>>>()) const {
    using Complex = std::ranges::range_value_t<Range>;
    using Real = RemoveComplex<Complex>;
    assert(range.size() == ComplexCoefficientSize(lMax, n));
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::ranges::generate(
        range, [&gen, &dist]() { return Complex{dist(gen), dist(gen)}; });

    if (lMax == Derived().MaxDegree()) {
      auto i = GSHIndices<All>(lMax, lMax, n).Index(lMax, lMax);
      range[i] = 0;
    }
  }

  template <std::ranges::range Range,
            typename Distribution =
                decltype(std::normal_distribution<
                         RemoveComplex<std::ranges::range_value_t<Range>>>())>
  requires requires() {
    requires ComplexFloatingPoint<std::ranges::range_value_t<Range>>;
    requires std::ranges::output_range<Range,
                                       std::ranges::range_value_t<Range>>;
  }
  void RandomRealCoefficient(
      Int lMax, Int n, Range& range,
      Distribution dist = std::normal_distribution<
          RemoveComplex<std::ranges::range_value_t<Range>>>()) const {
    using Complex = std::ranges::range_value_t<Range>;
    using Real = RemoveComplex<Complex>;
    assert(range.size() == RealCoefficientSize(lMax, n));
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::ranges::generate(
        range, [&gen, &dist]() { return Complex{dist(gen), dist(gen)}; });

    auto indices = GSHIndices<NonNegative>(lMax, lMax, n);
    for (auto l : indices.Degrees()) {
      auto i = indices.Index(l, 0);
      range[i].imag(0);
    }
    if (lMax == Derived().MaxDegree()) {
      auto i = indices.Index(lMax, lMax);
      range[i].imag(0);
    }
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GRID_GUARD_H