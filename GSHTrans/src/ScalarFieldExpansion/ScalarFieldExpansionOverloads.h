#ifndef GSH_TRANS_SCALAR_FIELD_EXPANSION_OVERLOADS_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_EXPANSION_OVERLOADS_GUARD_H

#include <FFTWpp/Core>
#include <algorithm>
#include <concepts>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "../Concepts.h"
#include "../ExpansionBase.h"
#include "../GridBase.h"
#include "ScalarFieldExpansionBase.h"
#include "ScalarFieldExpansionUnary.h"

namespace GSHTrans {

// Negation.
template <typename _Derived>
auto operator-(const ScalarFieldExpansionBase<_Derived>& u) {
  return ScalarFieldExpansionUnary(u, [](auto x) { return -x; });
}

template <typename _Derived>
auto operator-(ScalarFieldExpansionBase<_Derived>&& u) {
  return -u;
}

// Scalar multiplication.
template <typename _Derived>
auto operator*(const ScalarFieldExpansionBase<_Derived>& u,
               typename _Derived::Scalar s) {
  return ScalarFieldExpansionBinaryWithScalar(u, s, std::multiplies<>());
}

template <typename _Derived>
auto operator*(typename _Derived::Scalar s,
               const ScalarFieldExpansionBase<_Derived>& u) {
  return u * s;
}

template <typename _Derived>
auto operator*(ScalarFieldExpansionBase<_Derived>&& u,
               typename _Derived::Scalar s) {
  return u * s;
}

template <typename _Derived>
auto operator*(typename _Derived::Scalar s,
               ScalarFieldExpansionBase<_Derived>&& u) {
  return u * s;
}

// Scalar division.
template <typename _Derived>
auto operator/(const ScalarFieldExpansionBase<_Derived>& u,
               typename _Derived::Scalar s) {
  return ScalarFieldExpansionBinaryWithScalar(u, s, std::divides<>());
}

template <typename _Derived>
auto operator/(ScalarFieldExpansionBase<_Derived>&& u,
               typename _Derived::Scalar s) {
  return u / s;
}

// Addition.
template <typename _Derived1, typename _Derived2>
auto operator+(const ScalarFieldExpansionBase<_Derived1>& u1,
               const ScalarFieldExpansionBase<_Derived2>& u2) {
  return ScalarFieldExpansionBinary(u1, u2, std::plus<>());
}

template <typename _Derived1, typename _Derived2>
auto operator+(const ScalarFieldExpansionBase<_Derived1>& u1,
               ScalarFieldExpansionBase<_Derived2>&& u2) {
  return u1 + u2;
}

template <typename _Derived1, typename _Derived2>
auto operator+(ScalarFieldExpansionBase<_Derived1>&& u1,
               const ScalarFieldExpansionBase<_Derived2>& u2) {
  return u1 + u2;
}

template <typename _Derived1, typename _Derived2>
auto operator+(ScalarFieldExpansionBase<_Derived1>&& u1,
               ScalarFieldExpansionBase<_Derived2>&& u2) {
  return u1 + u2;
}

// Subtraction.
template <typename _Derived1, typename _Derived2>
auto operator-(const ScalarFieldExpansionBase<_Derived1>& u1,
               const ScalarFieldExpansionBase<_Derived2>& u2) {
  return ScalarFieldExpansionBinary(u1, u2, std::minus<>());
}

template <typename _Derived1, typename _Derived2>
auto operator-(const ScalarFieldExpansionBase<_Derived1>& u1,
               ScalarFieldExpansionBase<_Derived2>&& u2) {
  return u1 - u2;
}

template <typename _Derived1, typename _Derived2>
auto operator-(ScalarFieldExpansionBase<_Derived1>&& u1,
               const ScalarFieldExpansionBase<_Derived2>& u2) {
  return u1 - u2;
}

template <typename _Derived1, typename _Derived2>
auto operator-(ScalarFieldExpansionBase<_Derived1>&& u1,
               ScalarFieldExpansionBase<_Derived2>&& u2) {
  return u1 - u2;
}

// L2 inner product.
template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto L2InnerProduct(const ScalarFieldExpansionBase<_Derived1>& u1,
                    const ScalarFieldExpansionBase<_Derived2>& u2) {
  using Complex = typename _Derived1::Complex;
  auto sum = Complex(0);
  for (auto [l, m] : u1.Indices()) {
    if constexpr (std::same_as<typename _Derived1::Value, RealValued>) {
      auto scale = m == 0 ? static_cast<Complex>(1) : static_cast<Complex>(2);
      sum += scale * std::conj(u1[l, m]) * u2[l, m];
    } else {
      sum += std::conj(u1[l, m]) * u2[l, m];
    }
  }
  return sum;
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto L2InnerProduct(ScalarFieldExpansionBase<_Derived1>&& u1,
                    const ScalarFieldExpansionBase<_Derived2>& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto L2InnerProduct(const ScalarFieldExpansionBase<_Derived1>& u1,
                    ScalarFieldExpansionBase<_Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto L2InnerProduct(ScalarFieldExpansionBase<_Derived1>&& u1,
                    ScalarFieldExpansionBase<_Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

// L2 norm.
template <typename _Derived>
auto L2Norm(const ScalarFieldExpansionBase<_Derived>& u) {
  return std::sqrt(std::abs(L2InnerProduct(u, u)));
}

template <typename _Derived>
auto L2Norm(ScalarFieldExpansionBase<_Derived>&& u) {
  return L2Norm(u);
}

// Equality operator (based on L2 norm).
template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto operator==(const ScalarFieldExpansionBase<_Derived1>& u1,
                const ScalarFieldExpansionBase<_Derived2>& u2) {
  assert(u1.FieldSize() == u2.FieldSize());
  auto scale = std::max(L2Norm(u1), L2Norm(u2));
  return L2Norm(u1 - u2) <
         std::numeric_limits<typename _Derived1::Real>::epsilon() * scale;
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto operator==(ScalarFieldExpansionBase<_Derived1>&& u1,
                const ScalarFieldExpansionBase<_Derived2>& u2) {
  return u1 == u2;
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto operator==(const ScalarFieldExpansionBase<_Derived1>& u1,
                ScalarFieldExpansionBase<_Derived2>&& u2) {
  return u1 == u2;
}

template <typename _Derived1, typename _Derived2>
requires std::same_as<typename _Derived1::Value, typename _Derived2::Value>
auto operator==(ScalarFieldExpansionBase<_Derived1>&& u1,
                ScalarFieldExpansionBase<_Derived2>&& u2) {
  return u1 == u2;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_EXPRESSIONS_GUARD_H