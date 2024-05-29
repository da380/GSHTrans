#ifndef GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"

namespace GSHTrans {

namespace Internal {}  // namespace Internal

template <typename Derived>
class ScalarFieldBase : public FieldBase<ScalarFieldBase<Derived>> {
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
  auto operator()(Int iTheta, Int iPhi) const {
    return GetDerived().operator()(iTheta, iPhi);
  }

  // Write access to data.
  auto& operator()(Int iTheta, Int iPhi)
  requires Writeable::value
  {
    return GetDerived().operator()(iTheta, iPhi);
  }

  // Assign values from another field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator=(const ScalarFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      this->operator()(iTheta, iPhi) = other(iTheta, iPhi);
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator=(ScalarFieldBase<OtherDerived>&& other) {
    *this = other;
    return GetDerived();
  }

  // Compound plus assigment with other field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator+=(const ScalarFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      this->operator()(iTheta, iPhi) += other(iTheta, iPhi);
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator+=(ScalarFieldBase<OtherDerived>&& other) {
    *this += other;
    return GetDerived();
  }

  // Compound plus assigment with scalar.
  auto& operator+=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      this->operator()(iTheta, iPhi) += s;
    }
    return GetDerived();
  }

  // Compound minus assigment with other field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator-=(const ScalarFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      this->operator()(iTheta, iPhi) -= other(iTheta, iPhi);
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator-=(ScalarFieldBase<OtherDerived>&& other) {
    *this -= other;
    return GetDerived();
  }

  // Compound minus assigment with scalar.
  auto& operator-=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      this->operator()(iTheta, iPhi) -= s;
    }
    return GetDerived();
  }

  // Compound multiply assigment with other field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator*=(const ScalarFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      this->operator()(iTheta, iPhi) *= other(iTheta, iPhi);
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator*=(ScalarFieldBase<OtherDerived>&& other) {
    *this *= other;
    return GetDerived();
  }

  // Compound multiply assigment with scalar.
  auto& operator*=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      this->operator()(iTheta, iPhi) *= s;
    }
    return GetDerived();
  }

  // Compound divide assigment with other field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator/=(const ScalarFieldBase<OtherDerived>& other) {
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      this->operator()(iTheta, iPhi) /= other(iTheta, iPhi);
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator/=(ScalarFieldBase<OtherDerived>&& other) {
    *this /= other;
    return GetDerived();
  }

  // Compound multiply assigment with scalar.
  auto& operator/=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      this->operator()(iTheta, iPhi) /= s;
    }
    return GetDerived();
  }

  // Print the values.
  void Print() const {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      std::cout << operator()(iTheta, iPhi) << std::endl;
    }
  }

 private:
  auto& GetDerived() const { return static_cast<const Derived&>(*this); }
  auto& GetDerived() { return static_cast<Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H