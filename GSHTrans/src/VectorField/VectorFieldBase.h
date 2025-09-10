#ifndef GSH_TRANS_VECTOR_FIELD_BASE_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_BASE_GUARD_H

#include <concepts>
#include <iostream>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../ScalarField/ScalarFieldBase.h"
#include "CanonicalVector.h"

namespace GSHTrans {

template <typename Derived>
class VectorFieldBase : public FieldBase<VectorFieldBase<Derived>>,
                        public CanonicalVectorBase {
 public:
  using Int = typename Internal::Traits<Derived>::Int;
  using Grid = typename Internal::Traits<Derived>::Grid;
  using Value = typename Internal::Traits<Derived>::Value;
  using Real = typename Internal::Traits<Derived>::Real;
  using Complex = typename Internal::Traits<Derived>::Complex;
  using Scalar = typename Internal::Traits<Derived>::Scalar;
  using Writeable = typename Internal::Traits<Derived>::Writeable;

  // Return the grid.
  auto GetGrid() const { return GetDerived().GetGrid(); }

  // Read access to data.
  auto operator[](Int alpha, Int iTheta, Int iPhi) const {
    return GetDerived()[alpha, iTheta, iPhi];
  }

  // Write access to data.
  auto& operator[](Int alpha, Int iTheta, Int iPhi)
  requires Writeable::value
  {
    return GetDerived()[alpha, iTheta, iPhi];
  }

  //  Read-write access to component.
  auto operator[](Int alpha)
  requires Writeable::value
  {
    return GetDerived()[alpha];
  }

  // Read access to component.
  auto operator[](Int alpha) const { return GetDerived()[alpha]; }

  // Assign values from another field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator=(const VectorFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator[](alpha, iTheta, iPhi) = other[alpha, iTheta, iPhi];
      }
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator=(VectorFieldBase<OtherDerived>&& other) {
    *this = other;
    return GetDerived();
  }

  // Compound plus assigment with other field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator+=(const VectorFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator[](alpha, iTheta, iPhi) += other[alpha, iTheta, iPhi];
      }
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator+=(VectorFieldBase<OtherDerived>&& other) {
    GetDerived() += other;
    return GetDerived();
  }

  // Compound multiply by scalar field
  template <typename OtherDerived>
  requires Writeable::value && std::same_as<typename OtherDerived::Value, Value>
  auto& operator*=(const ScalarFieldBase<OtherDerived>& other) {
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator[](alpha, iTheta, iPhi) *= other[iTheta, iPhi];
      }
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value && std::same_as<typename OtherDerived::Value, Value>
  auto& operator*=(ScalarFieldBase<OtherDerived>&& other) {
    GetDerived() *= other;
    return GetDerived();
  }

  // Compound multiply by scalar
  auto& operator*=(Scalar s)
  requires Writeable::value
  {
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator[](alpha, iTheta, iPhi) *= s;
      }
    }
    return GetDerived();
  }

  // Compound divide by scalar field
  template <typename OtherDerived>
  requires Writeable::value && std::same_as<typename OtherDerived::Value, Value>
  auto& operator/=(const ScalarFieldBase<OtherDerived>& other) {
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator[](alpha, iTheta, iPhi) /= other[iTheta, iPhi];
      }
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value && std::same_as<typename OtherDerived::Value, Value>
  auto& operator/=(ScalarFieldBase<OtherDerived>&& other) {
    GetDerived() /= other;
    return GetDerived();
  }

  // Compound divide by scalar
  auto& operator/=(Scalar s)
  requires Writeable::value
  {
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator[](alpha, iTheta, iPhi) /= s;
      }
    }
    return GetDerived();
  }

  // Write values to ostream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const VectorFieldBase<Derived>& u) {
    auto indices = u.PointIndices();
    auto first = indices | std::ranges::views::reverse |
                 std::ranges::views::drop(1) | std::ranges::views::reverse;
    auto last =
        indices | std::ranges::views::reverse | std::ranges::views::take(1);
    for (auto [iTheta, iPhi] : first)
      os << u[-1, iTheta, iPhi] << " " << u[0, iTheta, iPhi] << " "
         << u[1, iTheta, iPhi] << std::endl;
    for (auto [iTheta, iPhi] : last)
      os << u[-1, iTheta, iPhi] << " " << u[0, iTheta, iPhi] << " "
         << u[1, iTheta, iPhi];
    return os;
  }

 private:
  auto& GetDerived() const { return static_cast<const Derived&>(*this); }
  auto& GetDerived() { return static_cast<Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_Vector_FIELD_BASE_GUARD_H