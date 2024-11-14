#ifndef GSH_TRANS_CANONICAL_COMPONENT_FIELD_BASE_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_FIELD_BASE_GUARD_H

#include <concepts>
#include <format>
#include <iostream>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"

namespace GSHTrans {

template <std::ptrdiff_t _N, typename _Derived>
class CanonicalComponentFieldBase
    : public FieldBase<CanonicalComponentFieldBase<_N, _Derived>> {
 public:
  using Int = typename Internal::Traits<_Derived>::Int;
  using Value = typename Internal::Traits<_Derived>::Value;
  using Real = typename Internal::Traits<_Derived>::Real;
  using Complex = typename Internal::Traits<_Derived>::Complex;
  using Scalar = typename Internal::Traits<_Derived>::Scalar;
  using Writeable = typename Internal::Traits<_Derived>::Writeable;

  // Return the upper index.
  constexpr auto UpperIndex() const { return _N; }

  // Return the grid.
  auto& Grid() const { return Derived().Grid(); }

  // Return the data size.
  auto Size() const {
    return this->NumberOfLongitudes() * this->NumberOfCoLatitudes();
  }

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

  // Assign values from another field.
  template <typename __Derived>
  requires Writeable::value &&
           std::convertible_to<typename __Derived::Scalar, Scalar>
  auto& operator=(const CanonicalComponentFieldBase<_N, __Derived>& other) {
    assert(other.Size() == Size());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) = other[iTheta, iPhi];
    }
    return Derived();
  }

  template <typename __Derived>
  requires Writeable::value &&
           std::convertible_to<typename __Derived::Scalar, Scalar>
  auto& operator=(CanonicalComponentFieldBase<_N, __Derived>&& other) {
    Derived() = other;
    return Derived();
  }

  // Compound plus assignment with other field.
  template <typename __Derived>
  requires Writeable::value &&
           std::convertible_to<typename __Derived::Scalar, Scalar>
  auto& operator+=(const CanonicalComponentFieldBase<_N, __Derived>& other) {
    assert(other.Size() == Size());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) += other[iTheta, iPhi];
    }
    return Derived();
  }

  template <typename __Derived>
  requires Writeable::value &&
           std::convertible_to<typename __Derived::Scalar, Scalar>
  auto& operator+=(CanonicalComponentFieldBase<_N, __Derived>&& other) {
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
  template <typename __Derived>
  requires Writeable::value &&
           std::convertible_to<typename __Derived::Scalar, Scalar>
  auto& operator-=(const CanonicalComponentFieldBase<_N, __Derived>& other) {
    assert(other.Size() == Size());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) -= other(iTheta, iPhi);
    }
    return Derived();
  }

  template <typename __Derived>
  requires Writeable::value &&
           std::convertible_to<typename __Derived::Scalar, Scalar>
  auto& operator-=(CanonicalComponentFieldBase<_N, __Derived>&& other) {
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
  template <typename __Derived>
  requires Writeable::value &&
           std::convertible_to<typename __Derived::Scalar, Scalar>
  auto& operator*=(const CanonicalComponentFieldBase<_N, __Derived>& other) {
    assert(other.Size() == Size());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) *= other(iTheta, iPhi);
    }
    return Derived();
  }

  template <typename __Derived>
  requires Writeable::value &&
           std::convertible_to<typename __Derived::Scalar, Scalar>
  auto& operator*=(CanonicalComponentFieldBase<_N, __Derived>&& other) {
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
  template <typename __Derived>
  requires Writeable::value &&
           std::convertible_to<typename __Derived::Scalar, Scalar>
  auto& operator/=(const CanonicalComponentFieldBase<_N, __Derived>& other) {
    assert(other.Size() == Size());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) /= other(iTheta, iPhi);
    }
    return Derived();
  }

  template <typename __Derived>
  requires Writeable::value &&
           std::convertible_to<typename __Derived::Scalar, Scalar>
  auto& operator/=(CanonicalComponentFieldBase<_N, __Derived>&& other) {
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
  friend std::ostream& operator<<(
      std::ostream& os, const CanonicalComponentFieldBase<_N, _Derived>& u) {
    auto values =
        u.PointIndices() | std::ranges::views::transform([&u](auto pair) {
          auto [iTheta, iPhi] = pair;
          return u[iTheta, iPhi];
        });
    auto range =
        std::ranges::views::zip(u.CoLatitudes(), u.Longitudes(), values);
    auto n = range.size();
    auto printValues = [&os](auto lat, auto lon, auto val, auto newLine) {
      if constexpr (std::same_as<typename _Derived::Value, ComplexValued>) {
        os << std::format("{:+.8e}  {:+.8e}  ({:+.8e},  {:+.8e})", lat, lon,
                          val.real(), val.imag());
      } else {
        os << std::format("{:+.8e}  {:+8e}  {:+.8e}", lat, lon, val);
      }
      if (newLine) os << '\n';
    };
    for (auto [lat, lon, val] : range | std::ranges::views::take(n - 1)) {
      printValues(lat, lon, val, true);
    }
    for (auto [lat, lon, val] : range | std::ranges::views::drop(n - 1)) {
      printValues(lat, lon, val, false);
    }
    return os;
  }

 private:
  constexpr auto& Derived() const {
    return static_cast<const _Derived&>(*this);
  }
  constexpr auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif