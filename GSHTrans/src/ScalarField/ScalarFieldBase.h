#ifndef GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H

#include <concepts>
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
  using GridType = typename Internal::Traits<_Derived>::GridType;
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
    std::cout << "Hello!" << std::endl;
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
    os << std::fixed;
    os << std::setprecision(8);

    auto first = u.View() | std::ranges::views::take(u.FieldSize() - 1);
    auto last = u.View() | std::ranges::views::drop(u.FieldSize() - 1);

    for (auto value : first) os << value << std::endl;
    for (auto value : last) os << value;

    /*
    for (auto [iTheta, iPhi] : first)
      os << theta[iTheta] << " " << phi[iPhi] << " " << u[iTheta, iPhi]
         << std::endl;
    for (auto [iTheta, iPhi] : last)
      os << theta[iTheta] << " " << phi[iPhi] << " " << u[iTheta, iPhi];

  */

    return os;
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H