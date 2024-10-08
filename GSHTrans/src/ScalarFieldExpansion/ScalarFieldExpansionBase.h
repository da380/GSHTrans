#ifndef GSH_TRANS_SCALAR_FIELD_EXPANSION_BASE_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_EXPANSION_BASE_GUARD_H

#include <concepts>
#include <format>
#include <iostream>

#include "../Concepts.h"
#include "../ExpansionBase.h"
#include "../GridBase.h"
#include "../Indexing.h"

namespace GSHTrans {

template <typename _Derived>
class ScalarFieldExpansionBase
    : public ExpansionBase<ScalarFieldExpansionBase<_Derived>> {
 public:
  using Int = typename Internal::Traits<_Derived>::Int;
  using Value = typename Internal::Traits<_Derived>::Value;
  using Real = typename Internal::Traits<_Derived>::Real;
  using Complex = typename Internal::Traits<_Derived>::Complex;
  using Scalar = typename Internal::Traits<_Derived>::Scalar;
  using Writeable = typename Internal::Traits<_Derived>::Writeable;

  // Return the grid.
  auto& Grid() const { return Derived().Grid(); }

  // Read access to the data.
  auto operator[](Int l, Int m) const { return Derived()[l, m]; }

  // Return view to the data.
  auto View() const {
    return Indices() | std::ranges::views::transform([this](auto index) {
             auto [l, m] = index;
             return operator[](l, m);
           });
  }

  // Return degrees.
  constexpr auto Degrees() const { _GSHIndices().Degrees(); }

  // Return spherical harmonic indices.
  constexpr auto Indices() const { return _GSHIndices().Indices(); }

  // Return index for degree and order.
  constexpr auto Index(Int l, Int m) const { return _GSHIndices().Index(l, m); }

  // Return size of Expansion vector.
  constexpr auto ExpansionSize() const { return _GSHIndices().size(); }

  // Write access to the data.
  auto operator[](Int l, Int m)
  requires Writeable::value
  {
    return Derived()[l, m];
  }

  // Assign values from another Expansion.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Complex, Complex>
  auto& operator=(const OtherDerived& other) {
    for (auto [l, m] : Indices()) {
      operator[](l, m) = other[l, m];
    }
    return Derived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Complex, Complex>
  auto& operator=(OtherDerived&& other) {
    Derived() = other;
    return Derived();
  }

  // Compound multiply assignment with scalar.
  auto& operator*=(Complex s)
  requires Writeable::value && std::same_as<Value, ComplexValued>
  {
    for (auto [l, m] : Indices()) {
      operator[](l, m) *= s;
    }
    return Derived();
  }

  auto& operator*=(Real s)
  requires Writeable::value
  {
    for (auto [l, m] : Indices()) {
      operator[](l, m) *= s;
    }
    return Derived();
  }

  // Compound divide assignment with scalar.
  auto& operator/=(Complex s)
  requires Writeable::value && std::same_as<Value, ComplexValued>
  {
    for (auto [l, m] : Indices()) {
      operator[](l, m) /= s;
    }
    return Derived();
  }

  auto& operator/=(Real s)
  requires Writeable::value
  {
    for (auto [l, m] : Indices()) {
      operator[](l, m) /= s;
    }
    return Derived();
  }

  // Write values to ostream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const ScalarFieldExpansionBase<_Derived>& u) {
    for (auto [index, val] : std::ranges::views::zip(u.Indices(), u.View())) {
      auto [l, m] = index;
      // os << std::format("{:+.8e}  {:+.8e}  ({:+.8e},  {:+.8e})\n", l, m,
      //                           val.real(), val.imag());
      os << std::format("({:4}, {:+4}) ({:+.8e}, {:+.8e})\n", l, m, val.real(),
                        val.imag());
    }
    return os;
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }

  // Return GSHIndices.
  constexpr auto _GSHIndices() const {
    return GSHIndices<
        std::conditional_t<std::same_as<Value, RealValued>, NonNegative, All>>(
        this->MaxDegree(), this->MaxDegree(), 0);
  }

  // Return GSHSubIndices for given degree.
  constexpr auto _GSHSubIndices(Int l) const { return _GSHIndices().Index(l); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_EXPANSION_BASE_GUARD_H