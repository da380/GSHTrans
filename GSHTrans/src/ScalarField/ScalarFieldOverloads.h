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
#include "ComplexifiedScalarField.h"
#include "ScalarFieldBase.h"
#include "ScalarFieldPointwiseBinary.h"
#include "ScalarFieldPointwiseUnary.h"
#include "ScalarFieldPointwiseUnaryWithScalar.h"

namespace GSHTrans {

// Complexification.
template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto Complexify(const ScalarFieldBase<Derived>& u) {
  return ComplexifiedScalarField(u);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto Complexify(ScalarFieldBase<Derived>&& u) {
  return Complexify(u);
}

// Negation.
template <typename Derived>
auto operator-(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldPointwiseUnary(u, [](auto x) { return -x; });
}

template <typename Derived>
auto operator-(ScalarFieldBase<Derived>&& u) {
  return -u;
}

// Square root.
template <typename Derived>
auto sqrt(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldPointwiseUnary(u, [](auto x) { return std::sqrt(x); });
}

template <typename Derived>
auto sqrt(ScalarFieldBase<Derived>&& u) {
  return sqrt(u);
}

// absolute value.
template <typename Derived>
auto abs(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldPointwiseUnary(u, [](auto x) { return std::abs(x); });
}

template <typename Derived>
auto abs(ScalarFieldBase<Derived>&& u) {
  return abs(u);
}

// Real part.
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldPointwiseUnary(u, [](auto x) { return std::real(x); });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(ScalarFieldBase<Derived>&& u) {
  return real(u);
}

// Imaginary part.
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto imag(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldPointwiseUnary(u, [](auto x) { return std::imag(x); });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto imag(ScalarFieldBase<Derived>&& u) {
  return imag(u);
}

// Complex conjugation.
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldPointwiseUnary(u, [](auto x) { return std::conj(x); });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(ScalarFieldBase<Derived>&& u) {
  return conj(u);
}

// Scalar multiplication.
template <typename Derived>
auto operator*(const ScalarFieldBase<Derived>& u, typename Derived::Scalar s) {
  return ScalarFieldPointwiseUnaryWithScalar(u, std::multiplies<>(), s);
}

template <typename Derived>
auto operator*(typename Derived::Scalar s, const ScalarFieldBase<Derived>& u) {
  return u * s;
}

template <typename Derived>
auto operator*(ScalarFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return u * s;
}

template <typename Derived>
auto operator*(typename Derived::Scalar s, ScalarFieldBase<Derived>&& u) {
  return u * s;
}

// Scalar addition.
template <typename Derived>
auto operator+(const ScalarFieldBase<Derived>& u, typename Derived::Scalar s) {
  return ScalarFieldPointwiseUnaryWithScalar(u, std::plus<>(), s);
}

template <typename Derived>
auto operator+(typename Derived::Scalar s, const ScalarFieldBase<Derived>& u) {
  return u + s;
}

template <typename Derived>
auto operator+(ScalarFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return u + s;
}

template <typename Derived>
auto operator+(typename Derived::Scalar s, ScalarFieldBase<Derived>&& u) {
  return u + s;
}

// Scalar subtraction.
template <typename Derived>
auto operator-(const ScalarFieldBase<Derived>& u, typename Derived::Scalar s) {
  return ScalarFieldPointwiseUnaryWithScalar(u, std::minus<>(), s);
}

template <typename Derived>
auto operator-(ScalarFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return u - s;
}

// Scalar division.
template <typename Derived>
auto operator/(const ScalarFieldBase<Derived>& u, typename Derived::Scalar s) {
  return ScalarFieldPointwiseUnaryWithScalar(u, std::divides<>(), s);
}

template <typename Derived>
auto operator/(ScalarFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return u / s;
}

// Raise to scalar power.
template <typename Derived>
auto pow(const ScalarFieldBase<Derived>& u, typename Derived::Scalar s) {
  return ScalarFieldPointwiseUnaryWithScalar(
      u, [](auto x, auto y) { return std::pow(x, y); }, s);
}

template <typename Derived>
auto pow(ScalarFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return pow(u, s);
}

// Addition.
template <typename Derived1, typename Derived2>
auto operator+(const ScalarFieldBase<Derived1>& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return ScalarFieldPointwiseBinary(u1, u2, std::plus<>());
}

template <typename Derived1, typename Derived2>
auto operator+(const ScalarFieldBase<Derived1>& u1,
               ScalarFieldBase<Derived2>&& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(ScalarFieldBase<Derived1>&& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(ScalarFieldBase<Derived1>&& u1, ScalarFieldBase<Derived2>&& u2) {
  return u1 + u2;
}

// Subtraction.
template <typename Derived1, typename Derived2>
auto operator-(const ScalarFieldBase<Derived1>& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return ScalarFieldPointwiseBinary(u1, u2, std::minus<>());
}

template <typename Derived1, typename Derived2>
auto operator-(const ScalarFieldBase<Derived1>& u1,
               ScalarFieldBase<Derived2>&& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(ScalarFieldBase<Derived1>&& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(ScalarFieldBase<Derived1>&& u1, ScalarFieldBase<Derived2>&& u2) {
  return u1 - u2;
}

// Multiplication.
template <typename Derived1, typename Derived2>
auto operator*(const ScalarFieldBase<Derived1>& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return ScalarFieldPointwiseBinary(u1, u2, std::multiplies<>());
}

template <typename Derived1, typename Derived2>
auto operator*(const ScalarFieldBase<Derived1>& u1,
               ScalarFieldBase<Derived2>&& u2) {
  return u1 * u2;
}

template <typename Derived1, typename Derived2>
auto operator*(ScalarFieldBase<Derived1>&& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return u1 * u2;
}

template <typename Derived1, typename Derived2>
auto operator*(ScalarFieldBase<Derived1>&& u1, ScalarFieldBase<Derived2>&& u2) {
  return u1 * u2;
}

// Division.
template <typename Derived1, typename Derived2>
auto operator/(const ScalarFieldBase<Derived1>& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return ScalarFieldPointwiseBinary(u1, u2, std::divides<>());
}

template <typename Derived1, typename Derived2>
auto operator/(const ScalarFieldBase<Derived1>& u1,
               ScalarFieldBase<Derived2>&& u2) {
  return u1 / u2;
}

template <typename Derived1, typename Derived2>
auto operator/(ScalarFieldBase<Derived1>&& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return u1 / u2;
}

template <typename Derived1, typename Derived2>
auto operator/(ScalarFieldBase<Derived1>&& u1, ScalarFieldBase<Derived2>&& u2) {
  return u1 / u2;
}

// Integation.
template <typename Derived>
auto Integrate(const ScalarFieldBase<Derived>& u) {
  using Scalar = typename Derived::Scalar;
  auto w = u.Weights();
  auto i = 0;
  auto sum = Scalar{0};
  for (auto [iTheta, iPhi] : u.PointIndices()) {
    sum += u(iTheta, iPhi) * w[i++];
  }
  return sum;
}

template <typename Derived>
auto Integrate(ScalarFieldBase<Derived>&& u) {
  return Integrate(u);
}

// L2 inner product.
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(const ScalarFieldBase<Derived1>& u1,
                    const ScalarFieldBase<Derived2>& u2) {
  if constexpr (std::same_as<typename Derived1::Value, RealValued>) {
    return Integrate(u1 * u2);
  } else {
    return Integrate(conj(u1) * u2);
  }
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(ScalarFieldBase<Derived1>&& u1,
                    const ScalarFieldBase<Derived2>& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(const ScalarFieldBase<Derived1>& u1,
                    ScalarFieldBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(ScalarFieldBase<Derived1>&& u1,
                    ScalarFieldBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

// L2 norm.
template <typename Derived>
auto L2Norm(const ScalarFieldBase<Derived>& u) {
  return std::sqrt(std::abs(L2InnerProduct(u, u)));
}

template <typename Derived>
auto L2Norm(ScalarFieldBase<Derived>&& u) {
  return L2Norm(u);
}

// Equality operator (based on L2 norm).
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(const ScalarFieldBase<Derived1>& u1,
                const ScalarFieldBase<Derived2>& u2) {
  assert(u1.FieldSize() == u2.FieldSize());
  return L2Norm(u1 - u2) <
         std::numeric_limits<typename Derived1::Real>::epsilon();
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(ScalarFieldBase<Derived1>&& u1,
                const ScalarFieldBase<Derived2>& u2) {
  return u1 == u2;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(const ScalarFieldBase<Derived1>& u1,
                ScalarFieldBase<Derived2>&& u2) {
  return u1 == u2;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(ScalarFieldBase<Derived1>&& u1,
                ScalarFieldBase<Derived2>&& u2) {
  return u1 == u2;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_EXPRESSIONS_GUARD_H