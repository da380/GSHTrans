#ifndef GSH_TRANS_SCALAR_FIELD_EXPANSION_OVERLOADS_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_EXPANSION_OVERLOADS_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../ExpansionBase.h"
#include "../GridBase.h"
#include "ScalarFieldExpansionBase.h"
#include "ScalarFieldExpansionUnary.h"

namespace GSHTrans {

// Negation.
template <typename Derived>
auto operator-(const ScalarFieldExpansionBase<Derived>& u) {
  return ScalarFieldExpansionUnary(u, [](auto x) { return -x; });
}

template <typename Derived>
auto operator-(ScalarFieldExpansionBase<Derived>&& u) {
  return -u;
}

/*


// Real to complex
template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(const ScalarFieldExpansionBase<Derived>& u) {
  return ScalarFieldExpansionUnary(
      u, [](auto x) -> typename Derived::Complex { return x; });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(ScalarFieldExpansionBase<Derived>&& u) {
  return complex(u);
}

*/

// Scalar multiplication.
template <typename Derived>
auto operator*(const ScalarFieldExpansionBase<Derived>& u,
               typename Derived::Scalar s) {
  return ScalarFieldExpansionBinaryWithScalar(u, s, std::multiplies<>());
}

template <typename Derived>
auto operator*(typename Derived::Scalar s,
               const ScalarFieldExpansionBase<Derived>& u) {
  return u * s;
}

template <typename Derived>
auto operator*(ScalarFieldExpansionBase<Derived>&& u,
               typename Derived::Scalar s) {
  return u * s;
}

template <typename Derived>
auto operator*(typename Derived::Scalar s,
               ScalarFieldExpansionBase<Derived>&& u) {
  return u * s;
}

// Scalar division.
template <typename Derived>
auto operator/(const ScalarFieldExpansionBase<Derived>& u,
               typename Derived::Scalar s) {
  return ScalarFieldExpansionBinaryWithScalar(u, s, std::divides<>());
}

template <typename Derived>
auto operator/(ScalarFieldExpansionBase<Derived>&& u,
               typename Derived::Scalar s) {
  return u / s;
}

// Addition.
template <typename Derived1, typename Derived2>
auto operator+(const ScalarFieldExpansionBase<Derived1>& u1,
               const ScalarFieldExpansionBase<Derived2>& u2) {
  return ScalarFieldExpansionBinary(u1, u2, std::plus<>());
}

template <typename Derived1, typename Derived2>
auto operator+(const ScalarFieldExpansionBase<Derived1>& u1,
               ScalarFieldExpansionBase<Derived2>&& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(ScalarFieldExpansionBase<Derived1>&& u1,
               const ScalarFieldExpansionBase<Derived2>& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(ScalarFieldExpansionBase<Derived1>&& u1,
               ScalarFieldExpansionBase<Derived2>&& u2) {
  return u1 + u2;
}

// Subtraction.
template <typename Derived1, typename Derived2>
auto operator-(const ScalarFieldExpansionBase<Derived1>& u1,
               const ScalarFieldExpansionBase<Derived2>& u2) {
  return ScalarFieldExpansionBinary(u1, u2, std::minus<>());
}

template <typename Derived1, typename Derived2>
auto operator-(const ScalarFieldExpansionBase<Derived1>& u1,
               ScalarFieldExpansionBase<Derived2>&& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(ScalarFieldExpansionBase<Derived1>&& u1,
               const ScalarFieldExpansionBase<Derived2>& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(ScalarFieldExpansionBase<Derived1>&& u1,
               ScalarFieldExpansionBase<Derived2>&& u2) {
  return u1 - u2;
}

/*
// Integation.
template <typename Derived>
auto Integrate(const ScalarFieldExpansionBase<Derived>& u) {
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
auto Integrate(ScalarFieldExpansionBase<Derived>&& u) {
  return Integrate(u);
}

// L2 inner product.
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(const ScalarFieldExpansionBase<Derived1>& u1,
                    const ScalarFieldExpansionBase<Derived2>& u2) {
  if constexpr (std::same_as<typename Derived1::Value, RealValued>) {
    return Integrate(u1 * u2);
  } else {
    return Integrate(conj(u1) * u2);
  }
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(ScalarFieldExpansionBase<Derived1>&& u1,
                    const ScalarFieldExpansionBase<Derived2>& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(const ScalarFieldExpansionBase<Derived1>& u1,
                    ScalarFieldExpansionBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(ScalarFieldExpansionBase<Derived1>&& u1,
                    ScalarFieldExpansionBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

// L2 norm.
template <typename Derived>
auto L2Norm(const ScalarFieldExpansionBase<Derived>& u) {
  return std::sqrt(std::abs(L2InnerProduct(u, u)));
}

template <typename Derived>
auto L2Norm(ScalarFieldExpansionBase<Derived>&& u) {
  return L2Norm(u);
}

// Equality operator (based on L2 norm).
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(const ScalarFieldExpansionBase<Derived1>& u1,
                const ScalarFieldExpansionBase<Derived2>& u2) {
  assert(u1.FieldSize() == u2.FieldSize());
  auto scale = std::max(L2Norm(u1), L2Norm(u2));
  return L2Norm(u1 - u2) <
         std::numeric_limits<typename Derived1::Real>::epsilon() * scale;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(ScalarFieldExpansionBase<Derived1>&& u1,
                const ScalarFieldExpansionBase<Derived2>& u2) {
  return u1 == u2;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(const ScalarFieldExpansionBase<Derived1>& u1,
                ScalarFieldExpansionBase<Derived2>&& u2) {
  return u1 == u2;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(ScalarFieldExpansionBase<Derived1>&& u1,
                ScalarFieldExpansionBase<Derived2>&& u2) {
  return u1 == u2;
}

*/

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_EXPRESSIONS_GUARD_H