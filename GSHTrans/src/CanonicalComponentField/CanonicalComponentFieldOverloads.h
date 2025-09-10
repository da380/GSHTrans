#ifndef GSH_TRANS_CANONICAL_COMPONENT_FIELD_OVERLOADS_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_FIELD_OVERLOADS_GUARD_H

#include <complex>
#include <concepts>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../Utility.h"
#include "CanonicalComponentFieldBase.h"
// #include "CanonicalComponentFieldBinary.h"
#include "CanonicalComponentFieldBinary.h"
#include "CanonicalComponentFieldUnary.h"

namespace GSHTrans {

// Overloads for unary minus.
template <std::ptrdiff_t N, typename Derived>
auto operator-(const CanonicalComponentFieldBase<N, Derived>& u) {
  return CanonicalComponentFieldUnary(u, [](auto x) { return -x; });
}

template <std::ptrdiff_t N, typename Derived>
auto operator-(CanonicalComponentFieldBase<N, Derived>&& u) {
  return -u;
}

// Overloads for complex conjugation.
template <std::ptrdiff_t N, typename Derived>
auto conj(const CanonicalComponentFieldBase<N, Derived>& u) {
  return CanonicalComponentFieldConj(u);
}

template <std::ptrdiff_t N, typename Derived>
auto conj(CanonicalComponentFieldBase<N, Derived>&& u) {
  return conj(u);
}

// Overloads for taking the real part.
template <std::ptrdiff_t N, typename Derived>
auto real(const CanonicalComponentFieldBase<N, Derived>& u) {
  return CanonicalComponentFieldReal(u);
}

template <std::ptrdiff_t N, typename Derived>
auto real(CanonicalComponentFieldBase<N, Derived>&& u) {
  return real(u);
}

// Overloads for taking the imaginary part.
template <std::ptrdiff_t N, typename Derived>
auto imag(const CanonicalComponentFieldBase<N, Derived>& u) {
  return CanonicalComponentFieldImag(u);
}

template <std::ptrdiff_t N, typename Derived>
auto imag(CanonicalComponentFieldBase<N, Derived>&& u) {
  return imag(u);
}

// Overloads for scalar addition.
template <std::ptrdiff_t N, typename Derived>
requires(N == 0)
auto operator+(const CanonicalComponentFieldBase<N, Derived>& u,
               typename Derived::Scalar s) {
  return CanonicalComponentFieldUnaryWithScalar(u, std::plus<>(), s);
}

template <std::ptrdiff_t N, typename Derived>
requires(N == 0)
auto operator+(CanonicalComponentFieldBase<N, Derived>&& u,
               typename Derived::Scalar s) {
  return u + s;
}

template <std::ptrdiff_t N, typename Derived>
requires(N == 0)
auto operator+(typename Derived::Scalar s,
               const CanonicalComponentFieldBase<N, Derived>& u) {
  return u + s;
}

template <std::ptrdiff_t N, typename Derived>
requires(N == 0)
auto operator+(typename Derived::Scalar s,
               CanonicalComponentFieldBase<N, Derived>&& u) {
  return u + s;
}

// Overloads for scalar multiplication.
template <std::ptrdiff_t N, typename Derived>
auto operator*(const CanonicalComponentFieldBase<N, Derived>& u,
               typename Derived::Scalar s) {
  return CanonicalComponentFieldUnaryWithScalar(u, std::multiplies<>(), s);
}

template <std::ptrdiff_t N, typename Derived>
auto operator*(CanonicalComponentFieldBase<N, Derived>&& u,
               typename Derived::Scalar s) {
  return u * s;
}

template <std::ptrdiff_t N, typename Derived>
auto operator*(typename Derived::Scalar s,
               const CanonicalComponentFieldBase<N, Derived>& u) {
  return u * s;
}

template <std::ptrdiff_t N, typename Derived>
auto operator*(typename Derived::Scalar s,
               CanonicalComponentFieldBase<N, Derived>&& u) {
  return u * s;
}

// Overloads for scalar subtraction.
template <std::ptrdiff_t N, typename Derived>
requires(N == 0)
auto operator-(const CanonicalComponentFieldBase<N, Derived>& u,
               typename Derived::Scalar s) {
  return CanonicalComponentFieldUnaryWithScalar(u, std::minus<>(), s);
}

template <std::ptrdiff_t N, typename Derived>
requires(N == 0)
auto operator-(CanonicalComponentFieldBase<N, Derived>&& u,
               typename Derived::Scalar s) {
  return u - s;
}

// Overloads for scalar division.
template <std::ptrdiff_t N, typename Derived>
auto operator/(const CanonicalComponentFieldBase<N, Derived>& u,
               typename Derived::Scalar s) {
  return CanonicalComponentFieldUnaryWithScalar(u, std::divides<>(), s);
}

template <std::ptrdiff_t N, typename Derived>
auto operator/(CanonicalComponentFieldBase<N, Derived>&& u,
               typename Derived::Scalar s) {
  return u / s;
}

// Overloads for addition.
template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator+(const CanonicalComponentFieldBase<N0, Derived0>& u0,
               const CanonicalComponentFieldBase<N1, Derived1>& u1) {
  return CanonicalComponentFieldAdd(u0, u1);
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator+(const CanonicalComponentFieldBase<N0, Derived0>& u0,
               CanonicalComponentFieldBase<N1, Derived1>&& u1) {
  return u0 + u1;
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator+(CanonicalComponentFieldBase<N0, Derived0>&& u0,
               const CanonicalComponentFieldBase<N1, Derived1>& u1) {
  return u0 + u1;
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator+(CanonicalComponentFieldBase<N0, Derived0>&& u0,
               CanonicalComponentFieldBase<N1, Derived1>&& u1) {
  return u0 + u1;
}

// Overloads for subtraction.
template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator-(const CanonicalComponentFieldBase<N0, Derived0>& u0,
               const CanonicalComponentFieldBase<N1, Derived1>& u1) {
  return CanonicalComponentFieldSubtract(u0, u1);
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator-(const CanonicalComponentFieldBase<N0, Derived0>& u0,
               CanonicalComponentFieldBase<N1, Derived1>&& u1) {
  return u0 - u1;
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator-(CanonicalComponentFieldBase<N0, Derived0>&& u0,
               const CanonicalComponentFieldBase<N1, Derived1>& u1) {
  return u0 - u1;
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator-(CanonicalComponentFieldBase<N0, Derived0>&& u0,
               CanonicalComponentFieldBase<N1, Derived1>&& u1) {
  return u0 - u1;
}

// Overloads for multiplication.
template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator*(const CanonicalComponentFieldBase<N0, Derived0>& u0,
               const CanonicalComponentFieldBase<N1, Derived1>& u1) {
  return CanonicalComponentFieldMultiply(u0, u1);
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator*(const CanonicalComponentFieldBase<N0, Derived0>& u0,
               CanonicalComponentFieldBase<N1, Derived1>&& u1) {
  return u0 * u1;
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator*(CanonicalComponentFieldBase<N0, Derived0>&& u0,
               const CanonicalComponentFieldBase<N1, Derived1>& u1) {
  return u0 * u1;
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator*(CanonicalComponentFieldBase<N0, Derived0>&& u0,
               CanonicalComponentFieldBase<N1, Derived1>&& u1) {
  return u0 * u1;
}

// Overloads for division.
template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator/(const CanonicalComponentFieldBase<N0, Derived0>& u0,
               const CanonicalComponentFieldBase<N1, Derived1>& u1) {
  return CanonicalComponentFieldDivide(u0, u1);
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator/(const CanonicalComponentFieldBase<N0, Derived0>& u0,
               CanonicalComponentFieldBase<N1, Derived1>&& u1) {
  return u0 / u1;
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator/(CanonicalComponentFieldBase<N0, Derived0>&& u0,
               const CanonicalComponentFieldBase<N1, Derived1>& u1) {
  return u0 / u1;
}

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
auto operator/(CanonicalComponentFieldBase<N0, Derived0>&& u0,
               CanonicalComponentFieldBase<N1, Derived1>&& u1) {
  return u0 / u1;
}

}  // namespace GSHTrans

#endif
