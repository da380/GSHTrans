#ifndef GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H

#include <concepts>
#include <format>
#include <iomanip>
#include <iostream>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"

namespace GSHTrans {

template <typename _Derived>
class ScalarFieldBase : public FieldBase<ScalarFieldBase<_Derived>> {
 public:
  using Int = typename Internal::Traits<_Derived>::Int;
  using Value = typename Internal::Traits<_Derived>::Value;
  using Real = typename Internal::Traits<_Derived>::Real;
  using Complex = typename Internal::Traits<_Derived>::Complex;
  using Scalar = typename Internal::Traits<_Derived>::Scalar;
  using Writeable = typename Internal::Traits<_Derived>::Writeable;

  // Return the grid.
  auto& Grid() const { return Derived().Grid(); }

  // Read access to data.
  auto operator[](Int iTheta, Int iPhi) const {
    return Derived()[iTheta, iPhi];
  }

  // Write access to data.
  auto& operator[](Int iTheta, Int iPhi)
  requires Writeable::value
  {
    return Derived()[iTheta, iPhi];
  }

  // Return a view to the values.
  auto View() const {
    return Grid().PointIndices() |
           std::ranges::views::transform([this](auto pair) {
             auto [iTheta, iPhi] = pair;
             return this->operator[](iTheta, iPhi);
           });
  }

  // Assign values from another field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator=(const ScalarFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) = other[iTheta, iPhi];
    }
    return Derived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator=(ScalarFieldBase<OtherDerived>&& other) {
    Derived() = other;
    return Derived();
  }

  // Compound plus assignment with other field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator+=(const ScalarFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) += other[iTheta, iPhi];
    }
    return Derived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator+=(ScalarFieldBase<OtherDerived>&& other) {
    Derived() += other;
    return Derived();
  }

  // Compound plus assignment with scalar.
  auto& operator+=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) += s;
    }
    return Derived();
  }

  // Compound minus assignment with other field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator-=(const ScalarFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) -= other(iTheta, iPhi);
    }
    return Derived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator-=(ScalarFieldBase<OtherDerived>&& other) {
    Derived() -= other;
    return Derived();
  }

  // Compound minus assignment with scalar.
  auto& operator-=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) -= s;
    }
    return Derived();
  }

  // Compound multiply assignment with other field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator*=(const ScalarFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) *= other(iTheta, iPhi);
    }
    return Derived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator*=(ScalarFieldBase<OtherDerived>&& other) {
    Derived() *= other;
    return Derived();
  }

  // Compound multiply assignment with scalar.
  auto& operator*=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) *= s;
    }
    return Derived();
  }

  // Compound divide assignment with other field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator/=(const ScalarFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) /= other(iTheta, iPhi);
    }
    return Derived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator/=(ScalarFieldBase<OtherDerived>&& other) {
    Derived() /= other;
    return Derived();
  }

  // Compound divide assignment with scalar.
  auto& operator/=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) /= s;
    }
    return Derived();
  }

  // Write values to ostream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const ScalarFieldBase<_Derived>& u) {
    for (auto [lat, lon, val] :
         std::ranges::views::zip(u.CoLatitudes(), u.Longitudes(), u.View())) {
      if constexpr (std::same_as<typename _Derived::Value, ComplexValued>) {
        os << std::format("{:+.8e}  {:+.8e}  ({:+.8e},  {:+.8e})\n", lat, lon,
                          val.real(), val.imag());
      } else {
        os << std::format("{:+.8e}  {:+8e}  {:+.8e}\n", lat, lon, val);
      }
    }
    return os;
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H