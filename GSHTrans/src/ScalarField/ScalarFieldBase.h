#ifndef GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H

#include <concepts>
#include <iostream>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"

namespace GSHTrans {

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
  auto operator[](Int iTheta, Int iPhi) const {
    return GetDerived()[iTheta, iPhi];
  }

  // Write access to data.
  auto& operator[](Int iTheta, Int iPhi)
  requires Writeable::value
  {
    return GetDerived()[iTheta, iPhi];
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
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator=(ScalarFieldBase<OtherDerived>&& other) {
    GetDerived() = other;
    return GetDerived();
  }

  // Compound plus assigment with other field.
  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator+=(const ScalarFieldBase<OtherDerived>& other) {
    std::cout << "Hello!" << std::endl;
    assert(other.FieldSize() == this->FieldSize());
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) += other[iTheta, iPhi];
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator+=(ScalarFieldBase<OtherDerived>&& other) {
    GetDerived() += other;
    return GetDerived();
  }

  // Compound plus assigment with scalar.
  auto& operator+=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) += s;
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
      operator[](iTheta, iPhi) -= other(iTheta, iPhi);
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator-=(ScalarFieldBase<OtherDerived>&& other) {
    GetDerived() -= other;
    return GetDerived();
  }

  // Compound minus assigment with scalar.
  auto& operator-=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) -= s;
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
      operator[](iTheta, iPhi) *= other(iTheta, iPhi);
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator*=(ScalarFieldBase<OtherDerived>&& other) {
    GetDerived() *= other;
    return GetDerived();
  }

  // Compound multiply assigment with scalar.
  auto& operator*=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) *= s;
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
      operator[](iTheta, iPhi) /= other(iTheta, iPhi);
    }
    return GetDerived();
  }

  template <typename OtherDerived>
  requires Writeable::value &&
           std::convertible_to<typename OtherDerived::Scalar, Scalar>
  auto& operator/=(ScalarFieldBase<OtherDerived>&& other) {
    GetDerived() /= other;
    return GetDerived();
  }

  // Compound divide assigment with scalar.
  auto& operator/=(Scalar s)
  requires Writeable::value
  {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) /= s;
    }
    return GetDerived();
  }

  // Write values to ostream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const ScalarFieldBase<Derived>& u) {
    auto indices = u.PointIndices();
    auto first = indices | std::ranges::views::reverse |
                 std::ranges::views::drop(1) | std::ranges::views::reverse;
    auto last =
        indices | std::ranges::views::reverse | std::ranges::views::take(1);
    for (auto [iTheta, iPhi] : first) os << u[iTheta, iPhi] << std::endl;
    for (auto [iTheta, iPhi] : last) os << u[iTheta, iPhi];
    return os;
  }

 private:
  auto& GetDerived() const { return static_cast<const Derived&>(*this); }
  auto& GetDerived() { return static_cast<Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H