#ifndef GSH_TRANS_SCALAR_FIELD_OVERLOADS_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_OVERLOADS_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "ScalarFieldBase.h"
#include "ScalarFieldBinary.h"
#include "ScalarFieldBinaryWithScalar.h"
#include "ScalarFieldUnary.h"

namespace GSHTrans {

// Negation.
template <typename _Derived>
auto operator-(const ScalarFieldBase<_Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return -x; });
}

template <typename _Derived>
auto operator-(ScalarFieldBase<_Derived>&& u) {
  return -u;
}

// Square root.
template <typename _Derived>
auto sqrt(const ScalarFieldBase<_Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return std::sqrt(x); });
}

template <typename _Derived>
auto sqrt(ScalarFieldBase<_Derived>&& u) {
  return sqrt(u);
}

// absolute value.
template <typename _Derived>
auto abs(const ScalarFieldBase<_Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return std::abs(x); });
}

template <typename _Derived>
auto abs(ScalarFieldBase<_Derived>&& u) {
  return abs(u);
}

// Real part.
template <typename _Derived>
requires std::same_as<typename _Derived::Value, ComplexValued>
auto real(const ScalarFieldBase<_Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return std::real(x); });
}

template <typename _Derived>
requires std::same_as<typename _Derived::Value, ComplexValued>
auto real(ScalarFieldBase<_Derived>&& u) {
  return real(u);
}

// Imaginary part.
template <typename _Derived>
requires std::same_as<typename _Derived::Value, ComplexValued>
auto imag(const ScalarFieldBase<_Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return std::imag(x); });
}

template <typename _Derived>
requires std::same_as<typename _Derived::Value, ComplexValued>
auto imag(ScalarFieldBase<_Derived>&& u) {
  return imag(u);
}

// Real to complex
template <typename _Derived>
requires std::same_as<typename _Derived::Value, RealValued>
auto complex(const ScalarFieldBase<_Derived>& u) {
  return ScalarFieldUnary(
      u, [](auto x) -> typename _Derived::Complex { return x; });
}

template <typename _Derived>
requires std::same_as<typename _Derived::Value, RealValued>
auto complex(ScalarFieldBase<_Derived>&& u) {
  return complex(u);
}

// Complex conjugation.
template <typename _Derived>
requires std::same_as<typename _Derived::Value, ComplexValued>
auto conj(const ScalarFieldBase<_Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return std::conj(x); });
}

template <typename _Derived>
requires std::same_as<typename _Derived::Value, ComplexValued>
auto conj(ScalarFieldBase<_Derived>&& u) {
  return conj(u);
}

// Scalar multiplication.
template <typename _Derived>
auto operator*(const ScalarFieldBase<_Derived>& u,
               typename _Derived::Scalar s) {
  return ScalarFieldBinaryWithScalar(u, s, std::multiplies<>());
}

template <typename _Derived>
auto operator*(typename _Derived::Scalar s,
               const ScalarFieldBase<_Derived>& u) {
  return u * s;
}

template <typename _Derived>
auto operator*(ScalarFieldBase<_Derived>&& u, typename _Derived::Scalar s) {
  return u * s;
}

template <typename _Derived>
auto operator*(typename _Derived::Scalar s, ScalarFieldBase<_Derived>&& u) {
  return u * s;
}

// Scalar addition.
template <typename _Derived>
auto operator+(const ScalarFieldBase<_Derived>& u,
               typename _Derived::Scalar s) {
  return ScalarFieldBinaryWithScalar(u, s, std::plus<>());
}

template <typename _Derived>
auto operator+(typename _Derived::Scalar s,
               const ScalarFieldBase<_Derived>& u) {
  return u + s;
}

template <typename _Derived>
auto operator+(ScalarFieldBase<_Derived>&& u, typename _Derived::Scalar s) {
  return u + s;
}

template <typename _Derived>
auto operator+(typename _Derived::Scalar s, ScalarFieldBase<_Derived>&& u) {
  return u + s;
}

// Scalar subtraction.
template <typename _Derived>
auto operator-(const ScalarFieldBase<_Derived>& u,
               typename _Derived::Scalar s) {
  return ScalarFieldBinaryWithScalar(u, s, std::minus<>());
}

template <typename _Derived>
auto operator-(ScalarFieldBase<_Derived>&& u, typename _Derived::Scalar s) {
  return u - s;
}

// Scalar division.
template <typename _Derived>
auto operator/(const ScalarFieldBase<_Derived>& u,
               typename _Derived::Scalar s) {
  return ScalarFieldBinaryWithScalar(u, s, std::divides<>());
}

template <typename _Derived>
auto operator/(ScalarFieldBase<_Derived>&& u, typename _Derived::Scalar s) {
  return u / s;
}

// Raise to scalar power.
template <typename _Derived>
auto pow(const ScalarFieldBase<_Derived>& u, typename _Derived::Scalar s) {
  return ScalarFieldBinaryWithScalar(
      u, s, [](auto x, auto y) { return std::pow(x, y); });
}

template <typename _Derived>
auto pow(ScalarFieldBase<_Derived>&& u, typename _Derived::Scalar s) {
  return pow(u, s);
}

// Addition.
template <typename _Derived1, typename _Derived2>
auto operator+(const ScalarFieldBase<_Derived1>& u1,
               const ScalarFieldBase<_Derived2>& u2) {
  return ScalarFieldBinary(u1, u2, std::plus<>());
}

template <typename _Derived1, typename _Derived2>
auto operator+(const ScalarFieldBase<_Derived1>& u1,
               ScalarFieldBase<_Derived2>&& u2) {
  return u1 + u2;
}

template <typename _Derived1, typename _Derived2>
auto operator+(ScalarFieldBase<_Derived1>&& u1,
               const ScalarFieldBase<_Derived2>& u2) {
  return u1 + u2;
}

template <typename _Derived1, typename _Derived2>
auto operator+(ScalarFieldBase<_Derived1>&& u1,
               ScalarFieldBase<_Derived2>&& u2) {
  return u1 + u2;
}

// Subtraction.
template <typename _Derived1, typename _Derived2>
auto operator-(const ScalarFieldBase<_Derived1>& u1,
               const ScalarFieldBase<_Derived2>& u2) {
  return ScalarFieldBinary(u1, u2, std::minus<>());
}

template <typename _Derived1, typename _Derived2>
auto operator-(const ScalarFieldBase<_Derived1>& u1,
               ScalarFieldBase<_Derived2>&& u2) {
  return u1 - u2;
}

template <typename _Derived1, typename _Derived2>
auto operator-(ScalarFieldBase<_Derived1>&& u1,
               const ScalarFieldBase<_Derived2>& u2) {
  return u1 - u2;
}

template <typename _Derived1, typename _Derived2>
auto operator-(ScalarFieldBase<_Derived1>&& u1,
               ScalarFieldBase<_Derived2>&& u2) {
  return u1 - u2;
}

// Multiplication.
template <typename _Derived1, typename _Derived2>
auto operator*(const ScalarFieldBase<_Derived1>& u1,
               const ScalarFieldBase<_Derived2>& u2) {
  return ScalarFieldBinary(u1, u2, std::multiplies<>());
}

template <typename _Derived1, typename _Derived2>
auto operator*(const ScalarFieldBase<_Derived1>& u1,
               ScalarFieldBase<_Derived2>&& u2) {
  return u1 * u2;
}

template <typename _Derived1, typename _Derived2>
auto operator*(ScalarFieldBase<_Derived1>&& u1,
               const ScalarFieldBase<_Derived2>& u2) {
  return u1 * u2;
}

template <typename _Derived1, typename _Derived2>
auto operator*(ScalarFieldBase<_Derived1>&& u1,
               ScalarFieldBase<_Derived2>&& u2) {
  return u1 * u2;
}

// Division.
template <typename _Derived1, typename _Derived2>
auto operator/(const ScalarFieldBase<_Derived1>& u1,
               const ScalarFieldBase<_Derived2>& u2) {
  return ScalarFieldBinary(u1, u2, std::divides<>());
}

template <typename _Derived1, typename _Derived2>
auto operator/(const ScalarFieldBase<_Derived1>& u1,
               ScalarFieldBase<_Derived2>&& u2) {
  return u1 / u2;
}

template <typename _Derived1, typename _Derived2>
auto operator/(ScalarFieldBase<_Derived1>&& u1,
               const ScalarFieldBase<_Derived2>& u2) {
  return u1 / u2;
}

template <typename _Derived1, typename _Derived2>
auto operator/(ScalarFieldBase<_Derived1>&& u1,
               ScalarFieldBase<_Derived2>&& u2) {
  return u1 / u2;
}

// Integation.
template <typename _Derived>
auto Integrate(const ScalarFieldBase<_Derived>& u) {
  using Scalar = typename _Derived::Scalar;
  auto w = u.Weights();
  auto i = 0;
  auto sum = Scalar{0};
  for (auto [iTheta, iPhi] : u.PointIndices()) {
    sum += u[iTheta, iPhi] * w[i++];
  }
  return sum;
}

template <typename _Derived>
auto Integrate(ScalarFieldBase<_Derived>&& u) {
  return Integrate(u);
}

// L2 inner product.
template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto L2InnerProduct(const ScalarFieldBase<_Derived1>& u1,
                    const ScalarFieldBase<_Derived2>& u2) {
  if constexpr (std::same_as<typename _Derived1::Value, RealValued>) {
    return Integrate(u1 * u2);
  } else {
    return Integrate(conj(u1) * u2);
  }
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto L2InnerProduct(ScalarFieldBase<_Derived1>&& u1,
                    const ScalarFieldBase<_Derived2>& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto L2InnerProduct(const ScalarFieldBase<_Derived1>& u1,
                    ScalarFieldBase<_Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto L2InnerProduct(ScalarFieldBase<_Derived1>&& u1,
                    ScalarFieldBase<_Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

// L2 norm.
template <typename _Derived>
auto L2Norm(const ScalarFieldBase<_Derived>& u) {
  return std::sqrt(std::abs(L2InnerProduct(u, u)));
}

template <typename _Derived>
auto L2Norm(ScalarFieldBase<_Derived>&& u) {
  return L2Norm(u);
}

// Equality operator (based on L2 norm).
template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto operator==(const ScalarFieldBase<_Derived1>& u1,
                const ScalarFieldBase<_Derived2>& u2) {
  assert(u1.FieldSize() == u2.FieldSize());
  auto scale = std::max(L2Norm(u1), L2Norm(u2));
  return L2Norm(u1 - u2) <
         std::numeric_limits<typename _Derived1::Real>::epsilon() * scale;
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto operator==(ScalarFieldBase<_Derived1>&& u1,
                const ScalarFieldBase<_Derived2>& u2) {
  return u1 == u2;
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto operator==(const ScalarFieldBase<_Derived1>& u1,
                ScalarFieldBase<_Derived2>&& u2) {
  return u1 == u2;
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto operator==(ScalarFieldBase<_Derived1>&& u1,
                ScalarFieldBase<_Derived2>&& u2) {
  return u1 == u2;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_EXPRESSIONS_GUARD_H