#ifndef GSH_TRANS_CANONICAL_COMPONENT_FIELD_OVERLOADS_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_FIELD_OVERLOADS_GUARD_H

#include <complex>
#include <concepts>
#include <functional>
#include <type_traits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../Utility.h"
#include "CanonicalComponentFieldBase.h"
#include "CanonicalComponentFieldBinary.h"
#include "CanonicalComponentFieldUnary.h"

namespace GSHTrans {

// Overloads for unary minus.
template <typename Derived>
auto operator-(const CanonicalComponentFieldBase<Derived>& u) {
  return CanonicalComponentFieldUnary(u, [](auto x) { return -x; });
}

template <typename Derived>
auto operator-(CanonicalComponentFieldBase<Derived>&& u) {
  return -u;
}

// Overloads for complex conjugation.
template <typename Derived>
auto conj(const CanonicalComponentFieldBase<Derived>& u) {
  return CanonicalComponentFieldConj(u);
}

template <typename Derived>
auto conj(CanonicalComponentFieldBase<Derived>&& u) {
  return conj(u);
}

// Overloads for taking the real part.
template <typename Derived>
auto real(const CanonicalComponentFieldBase<Derived>& u) {
  return CanonicalComponentFieldReal(u);
}

template <typename Derived>
auto real(CanonicalComponentFieldBase<Derived>&& u) {
  return real(u);
}

// Overloads for taking the imaginary part.
template <typename Derived>
auto imag(const CanonicalComponentFieldBase<Derived>& u) {
  return CanonicalComponentFieldImag(u);
}

template <typename Derived>
auto imag(CanonicalComponentFieldBase<Derived>&& u) {
  return imag(u);
}

// Overloads for scalar addition.
template <typename Derived>
auto operator+(const CanonicalComponentFieldBase<Derived>& u,
               typename Derived::Scalar s) {
  return CanonicalComponentFieldUnaryWithScalar(u, std::plus<>(), s);
}

template <typename Derived>
auto operator+(CanonicalComponentFieldBase<Derived>&& u,
               typename Derived::Scalar s) {
  return u + s;
}

template <typename Derived>
auto operator+(typename Derived::Scalar s,
               const CanonicalComponentFieldBase<Derived>& u) {
  return u + s;
}

template <typename Derived>
auto operator+(typename Derived::Scalar s,
               CanonicalComponentFieldBase<Derived>&& u) {
  return u + s;
}

// Overloads for scalar multiplication.
template <typename Derived>
auto operator*(const CanonicalComponentFieldBase<Derived>& u,
               typename Derived::Scalar s) {
  return CanonicalComponentFieldUnaryWithScalar(u, std::multiplies<>(), s);
}

template <typename Derived>
auto operator*(CanonicalComponentFieldBase<Derived>&& u,
               typename Derived::Scalar s) {
  return u * s;
}

template <typename Derived>
auto operator*(typename Derived::Scalar s,
               const CanonicalComponentFieldBase<Derived>& u) {
  return u * s;
}

template <typename Derived>
auto operator*(typename Derived::Scalar s,
               CanonicalComponentFieldBase<Derived>&& u) {
  return u * s;
}

// Overloads for scalar subtraction.
template <typename Derived>
auto operator-(const CanonicalComponentFieldBase<Derived>& u,
               typename Derived::Scalar s) {
  return CanonicalComponentFieldUnaryWithScalar(u, std::minus<>(), s);
}

template <typename Derived>
auto operator-(CanonicalComponentFieldBase<Derived>&& u,
               typename Derived::Scalar s) {
  return u - s;
}

// Overloads for scalar division.
template <typename Derived>
auto operator/(const CanonicalComponentFieldBase<Derived>& u,
               typename Derived::Scalar s) {
  return CanonicalComponentFieldUnaryWithScalar(u, std::divides<>(), s);
}

template <typename Derived>
auto operator/(CanonicalComponentFieldBase<Derived>&& u,
               typename Derived::Scalar s) {
  return u / s;
}

// Overloads for addition.
template <typename Derived0, typename Derived1>
auto operator+(const CanonicalComponentFieldBase<Derived0>& u0,
               const CanonicalComponentFieldBase<Derived1>& u1) {
  assert(u0.UpperIndex() == u1.UpperIndex());
  return CanonicalComponentFieldBinary(u0, u1, std::plus<>(), u0.UpperIndex());
}

template <typename Derived0, typename Derived1>
auto operator+(const CanonicalComponentFieldBase<Derived0>& u0,
               CanonicalComponentFieldBase<Derived1>&& u1) {
  return u0 + u1;
}

template <typename Derived0, typename Derived1>
auto operator+(CanonicalComponentFieldBase<Derived0>&& u0,
               const CanonicalComponentFieldBase<Derived1>& u1) {
  return u0 + u1;
}

template <typename Derived0, typename Derived1>
auto operator+(CanonicalComponentFieldBase<Derived0>&& u0,
               CanonicalComponentFieldBase<Derived1>&& u1) {
  return u0 + u1;
}

// Overloads for subtraction.
template <typename Derived0, typename Derived1>
auto operator-(const CanonicalComponentFieldBase<Derived0>& u0,
               const CanonicalComponentFieldBase<Derived1>& u1) {
  assert(u0.UpperIndex() == u1.UpperIndex());
  return CanonicalComponentFieldBinary(u0, u1, std::minus<>(), u0.UpperIndex());
}

template <typename Derived0, typename Derived1>
auto operator-(const CanonicalComponentFieldBase<Derived0>& u0,
               CanonicalComponentFieldBase<Derived1>&& u1) {
  return u0 - u1;
}

template <typename Derived0, typename Derived1>
auto operator-(CanonicalComponentFieldBase<Derived0>&& u0,
               const CanonicalComponentFieldBase<Derived1>& u1) {
  return u0 - u1;
}

template <typename Derived0, typename Derived1>
auto operator-(CanonicalComponentFieldBase<Derived0>&& u0,
               CanonicalComponentFieldBase<Derived1>&& u1) {
  return u0 - u1;
}

// Overloads for multiplication.
template <typename Derived0, typename Derived1>
auto operator*(const CanonicalComponentFieldBase<Derived0>& u0,
               const CanonicalComponentFieldBase<Derived1>& u1) {
  return CanonicalComponentFieldBinary(u0, u1, std::multiplies<>(),
                                       u0.UpperIndex() + u1.UpperIndex());
}

template <typename Derived0, typename Derived1>
auto operator*(const CanonicalComponentFieldBase<Derived0>& u0,
               CanonicalComponentFieldBase<Derived1>&& u1) {
  return u0 * u1;
}

template <typename Derived0, typename Derived1>
auto operator*(CanonicalComponentFieldBase<Derived0>&& u0,
               const CanonicalComponentFieldBase<Derived1>& u1) {
  return u0 * u1;
}

template <typename Derived0, typename Derived1>
auto operator*(CanonicalComponentFieldBase<Derived0>&& u0,
               CanonicalComponentFieldBase<Derived1>&& u1) {
  return u0 * u1;
}

// Overloads for division.
template <typename Derived0, typename Derived1>
auto operator/(const CanonicalComponentFieldBase<Derived0>& u0,
               const CanonicalComponentFieldBase<Derived1>& u1) {
  assert(u1.UpperIndex() == 0);
  return CanonicalComponentFieldBinary(u0, u1, std::divides<>(),
                                       u0.UpperIndex());
}

template <typename Derived0, typename Derived1>
auto operator/(const CanonicalComponentFieldBase<Derived0>& u0,
               CanonicalComponentFieldBase<Derived1>&& u1) {
  return u0 / u1;
}

template <typename Derived0, typename Derived1>
auto operator/(CanonicalComponentFieldBase<Derived0>&& u0,
               const CanonicalComponentFieldBase<Derived1>& u1) {
  return u0 / u1;
}

template <typename Derived0, typename Derived1>
auto operator/(CanonicalComponentFieldBase<Derived0>&& u0,
               CanonicalComponentFieldBase<Derived1>&& u1) {
  return u0 / u1;
}

}  // namespace GSHTrans

#endif
