#ifndef GSH_TRANS_SCALAR_FIELD_COEFFICIENT_BASE_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_COEFFICIENT_BASE_GUARD_H

#include <concepts>
#include <iostream>

#include "../CoefficientBase.h"
#include "../Concepts.h"
#include "../GridBase.h"
#include "../Indexing.h"

namespace GSHTrans {

template <typename Derived>
class ScalarFieldCoefficientBase
    : public CoefficientBase<ScalarFieldCoefficientBase<Derived>> {
 public:
  using Int = typename Internal::Traits<Derived>::Int;
  using Grid = typename Internal::Traits<Derived>::Grid;
  using Value = typename Internal::Traits<Derived>::Value;
  using Real = typename Internal::Traits<Derived>::Real;
  using Complex = typename Internal::Traits<Derived>::Complex;
  using Writeable = typename Internal::Traits<Derived>::Writeable;

  // Return the grid.
  auto GetGrid() const { return GetDerived().GetGrid(); }

  // Read access to the data.
  auto operator[](Int l, Int m) const { return GetDerived()[l, m]; }

  // Return degrees.
  constexpr auto Degrees() const { _GSHIndices().Degrees(); }

  // Return spherical harmonic indices.
  constexpr auto Indices() const { return _GSHIndices().Indices(); }

  // Return index for degree and order.
  constexpr auto Index(Int l, Int m) const { return _GSHIndices().Index(l, m); }

  // Return size of coefficient vector.
  constexpr auto CoefficientSize() const { return _GSHIndices().size(); }

  // Write access to the data.
  auto operator[](Int l, Int m)
  requires Writeable::value
  {
    return GetDerived()[l, m];
  }

  // Assign values from another coefficient.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Complex, Complex>
  auto& operator=(const OtherDerived& other) {
    for (auto [l, m] : Indices()) {
      operator[](l, m) = other[l, m];
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Complex, Complex>
  auto& operator=(OtherDerived&& other) {
    GetDerived() = other;
    return GetDerived();
  }

  // Compound multiply assigment with scalar.
  auto& operator*=(Complex s)
  requires Writeable::value && std::same_as<Value, ComplexValued>
  {
    for (auto [l, m] : Indices()) {
      operator[](l, m) *= s;
    }
    return GetDerived();
  }

  auto& operator*=(Real s)
  requires Writeable::value
  {
    for (auto [l, m] : Indices()) {
      operator[](l, m) *= s;
    }
    return GetDerived();
  }

  // Compound divide assigment with scalar.
  auto& operator/=(Complex s)
  requires Writeable::value && std::same_as<Value, ComplexValued>
  {
    for (auto [l, m] : Indices()) {
      operator[](l, m) /= s;
    }
    return GetDerived();
  }

  auto& operator/=(Real s)
  requires Writeable::value
  {
    for (auto [l, m] : Indices()) {
      operator[](l, m) /= s;
    }
    return GetDerived();
  }

  // Write values to ostream.
  friend std::ostream& operator<<(
      std::ostream& os, const ScalarFieldCoefficientBase<Derived>& u) {
    for (auto [l, m] : u.Indices()) {
      os << u[l, m];
      if (l < u.MaxDegree() || m < u.MaxDegree()) os << std::endl;
    }
    return os;
  }

 private:
  auto& GetDerived() const { return static_cast<const Derived&>(*this); }
  auto& GetDerived() { return static_cast<Derived&>(*this); }

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

#endif  // GSH_TRANS_SCALAR_FIELD_COEFFICIENT_BASE_GUARD_H