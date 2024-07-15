#ifndef GSH_TRANS_SCALAR_FIELD_OVERLOADS_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_OVERLOADS_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../CoefficientBase.h"
#include "../Concepts.h"
#include "../GridBase.h"
#include "ScalarFieldCoefficientBase.h"
#include "ScalarFieldCoefficientUnary.h"

namespace GSHTrans {

// Negation.
template <typename Derived>
auto operator-(const ScalarFieldCoefficientBase<Derived>& u) {
  return ScalarFieldCoefficientUnary(u, [](auto x) { return -x; });
}

template <typename Derived>
auto operator-(ScalarFieldCoefficientBase<Derived>&& u) {
  return -u;
}

/*

// Real part.
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(const ScalarFieldCoefficientBase<Derived>& u) {
  return ScalarFieldCoefficientUnary(u, [](auto x) { return std::real(x); });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(ScalarFieldCoefficientBase<Derived>&& u) {
  return real(u);
}

// Imaginary part.
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto imag(const ScalarFieldCoefficientBase<Derived>& u) {
  return ScalarFieldCoefficientUnary(u, [](auto x) { return std::imag(x); });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto imag(ScalarFieldCoefficientBase<Derived>&& u) {
  return imag(u);
}

// Real to complex
template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(const ScalarFieldCoefficientBase<Derived>& u) {
  return ScalarFieldCoefficientUnary(
      u, [](auto x) -> typename Derived::Complex { return x; });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(ScalarFieldCoefficientBase<Derived>&& u) {
  return complex(u);
}

// Complex conjugation.
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(const ScalarFieldCoefficientBase<Derived>& u) {
  return ScalarFieldCoefficientUnary(u, [](auto x) { return std::conj(x); });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(ScalarFieldCoefficientBase<Derived>&& u) {
  return conj(u);
}

*/

/*

// Scalar multiplication.
template <typename Derived>
auto operator*(const ScalarFieldCoefficientBase<Derived>& u,
               typename Derived::Scalar s) {
  return ScalarFieldCoefficientBinaryWithScalar(u, s, std::multiplies<>());
}

template <typename Derived>
auto operator*(typename Derived::Scalar s,
               const ScalarFieldCoefficientBase<Derived>& u) {
  return u * s;
}

template <typename Derived>
auto operator*(ScalarFieldCoefficientBase<Derived>&& u,
               typename Derived::Scalar s) {
  return u * s;
}

template <typename Derived>
auto operator*(typename Derived::Scalar s,
               ScalarFieldCoefficientBase<Derived>&& u) {
  return u * s;
}

// Scalar addition.
template <typename Derived>
auto operator+(const ScalarFieldCoefficientBase<Derived>& u,
               typename Derived::Scalar s) {
  return ScalarFieldCoefficientBinaryWithScalar(u, s, std::plus<>());
}

template <typename Derived>
auto operator+(typename Derived::Scalar s,
               const ScalarFieldCoefficientBase<Derived>& u) {
  return u + s;
}

template <typename Derived>
auto operator+(ScalarFieldCoefficientBase<Derived>&& u,
               typename Derived::Scalar s) {
  return u + s;
}

template <typename Derived>
auto operator+(typename Derived::Scalar s,
               ScalarFieldCoefficientBase<Derived>&& u) {
  return u + s;
}

// Scalar subtraction.
template <typename Derived>
auto operator-(const ScalarFieldCoefficientBase<Derived>& u,
               typename Derived::Scalar s) {
  return ScalarFieldCoefficientBinaryWithScalar(u, s, std::minus<>());
}

template <typename Derived>
auto operator-(ScalarFieldCoefficientBase<Derived>&& u,
               typename Derived::Scalar s) {
  return u - s;
}

// Scalar division.
template <typename Derived>
auto operator/(const ScalarFieldCoefficientBase<Derived>& u,
               typename Derived::Scalar s) {
  return ScalarFieldCoefficientBinaryWithScalar(u, s, std::divides<>());
}

template <typename Derived>
auto operator/(ScalarFieldCoefficientBase<Derived>&& u,
               typename Derived::Scalar s) {
  return u / s;
}

// Addition.
template <typename Derived1, typename Derived2>
auto operator+(const ScalarFieldCoefficientBase<Derived1>& u1,
               const ScalarFieldCoefficientBase<Derived2>& u2) {
  return ScalarFieldCoefficientBinary(u1, u2, std::plus<>());
}

template <typename Derived1, typename Derived2>
auto operator+(const ScalarFieldCoefficientBase<Derived1>& u1,
               ScalarFieldCoefficientBase<Derived2>&& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(ScalarFieldCoefficientBase<Derived1>&& u1,
               const ScalarFieldCoefficientBase<Derived2>& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(ScalarFieldCoefficientBase<Derived1>&& u1,
               ScalarFieldCoefficientBase<Derived2>&& u2) {
  return u1 + u2;
}

// Subtraction.
template <typename Derived1, typename Derived2>
auto operator-(const ScalarFieldCoefficientBase<Derived1>& u1,
               const ScalarFieldCoefficientBase<Derived2>& u2) {
  return ScalarFieldCoefficientBinary(u1, u2, std::minus<>());
}

template <typename Derived1, typename Derived2>
auto operator-(const ScalarFieldCoefficientBase<Derived1>& u1,
               ScalarFieldCoefficientBase<Derived2>&& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(ScalarFieldCoefficientBase<Derived1>&& u1,
               const ScalarFieldCoefficientBase<Derived2>& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(ScalarFieldCoefficientBase<Derived1>&& u1,
               ScalarFieldCoefficientBase<Derived2>&& u2) {
  return u1 - u2;
}


*/

/*
// Integation.
template <typename Derived>
auto Integrate(const ScalarFieldCoefficientBase<Derived>& u) {
  using Scalar = typename Derived::Scalar;
  auto w = u.Weights();
  auto i = 0;
  auto sum = Scalar{0};
  for (auto [iTheta, iPhi] : u.PointIndices()) {
    sum += u[iTheta, iPhi] * w[i++];
  }
  return sum;
}

template <typename Derived>
auto Integrate(ScalarFieldCoefficientBase<Derived>&& u) {
  return Integrate(u);
}

// L2 inner product.
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(const ScalarFieldCoefficientBase<Derived1>& u1,
                    const ScalarFieldCoefficientBase<Derived2>& u2) {
  if constexpr (std::same_as<typename Derived1::Value, RealValued>) {
    return Integrate(u1 * u2);
  } else {
    return Integrate(conj(u1) * u2);
  }
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(ScalarFieldCoefficientBase<Derived1>&& u1,
                    const ScalarFieldCoefficientBase<Derived2>& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(const ScalarFieldCoefficientBase<Derived1>& u1,
                    ScalarFieldCoefficientBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(ScalarFieldCoefficientBase<Derived1>&& u1,
                    ScalarFieldCoefficientBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

// L2 norm.
template <typename Derived>
auto L2Norm(const ScalarFieldCoefficientBase<Derived>& u) {
  return std::sqrt(std::abs(L2InnerProduct(u, u)));
}

template <typename Derived>
auto L2Norm(ScalarFieldCoefficientBase<Derived>&& u) {
  return L2Norm(u);
}

// Equality operator (based on L2 norm).
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(const ScalarFieldCoefficientBase<Derived1>& u1,
                const ScalarFieldCoefficientBase<Derived2>& u2) {
  assert(u1.FieldSize() == u2.FieldSize());
  auto scale = std::max(L2Norm(u1), L2Norm(u2));
  return L2Norm(u1 - u2) <
         std::numeric_limits<typename Derived1::Real>::epsilon() * scale;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(ScalarFieldCoefficientBase<Derived1>&& u1,
                const ScalarFieldCoefficientBase<Derived2>& u2) {
  return u1 == u2;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(const ScalarFieldCoefficientBase<Derived1>& u1,
                ScalarFieldCoefficientBase<Derived2>&& u2) {
  return u1 == u2;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(ScalarFieldCoefficientBase<Derived1>&& u1,
                ScalarFieldCoefficientBase<Derived2>&& u2) {
  return u1 == u2;
}

*/

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_EXPRESSIONS_GUARD_H