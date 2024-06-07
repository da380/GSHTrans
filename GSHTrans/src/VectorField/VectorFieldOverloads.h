#ifndef GSH_TRANS_VECTOR_FIELD_OVERLOADS_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_OVERLOADS_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "VectorFieldBase.h"
#include "VectorFieldBinary.h"
#include "VectorFieldBinaryWithScalar.h"
#include "VectorFieldComplexToImag.h"
#include "VectorFieldComplexToReal.h"
#include "VectorFieldRealToComplex.h"
#include "VectorFieldUnary.h"

namespace GSHTrans {

// Negation.
template <typename Derived>
auto operator-(const VectorFieldBase<Derived>& u) {
  return VectorFieldUnary(u, [](auto x) { return -x; });
}

template <typename Derived>
auto operator-(VectorFieldBase<Derived>&& u) {
  return -u;
}

// Map real to complex.
template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(const VectorFieldBase<Derived>& u) {
  return VectorFieldRealToComplex(u);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(VectorFieldBase<Derived>&& u) {
  return complex(u);
}

// Real part.
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(const VectorFieldBase<Derived>& u) {
  return VectorFieldComplexToReal(u);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(VectorFieldBase<Derived>&& u) {
  return real(u);
}

// Imaginary part.
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto imag(const VectorFieldBase<Derived>& u) {
  return VectorFieldComplexToImag(u);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto imag(VectorFieldBase<Derived>&& u) {
  return imag(u);
}

/*

// Complex conjugation.
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(const VectorFieldBase<Derived>& u) {
  return VectorFieldUnary(u, [](auto x) { return std::conj(x); });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(VectorFieldBase<Derived>&& u) {
  return conj(u);
}

*/

// Scalar multiplication.
template <typename Derived>
auto operator*(const VectorFieldBase<Derived>& u, typename Derived::Scalar s) {
  return VectorFieldBinaryWithScalar(u, std::multiplies<>(), s);
}

template <typename Derived>
auto operator*(typename Derived::Scalar s, const VectorFieldBase<Derived>& u) {
  return u * s;
}

template <typename Derived>
auto operator*(VectorFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return u * s;
}

template <typename Derived>
auto operator*(typename Derived::Scalar s, VectorFieldBase<Derived>&& u) {
  return u * s;
}

// Scalar division.
template <typename Derived>
auto operator/(const VectorFieldBase<Derived>& u, typename Derived::Scalar s) {
  return VectorFieldBinaryWithScalar(u, std::divides<>(), s);
}

template <typename Derived>
auto operator/(VectorFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return u / s;
}

// Scalar field multiplication.
template <typename Derived0, typename Derived1>
auto operator*(const VectorFieldBase<Derived0>& u,
               const ScalarFieldBase<Derived1>& f) {
  return VectorFieldBinaryWithScalarField(u, f, std::multiplies<>());
}

template <typename Derived0, typename Derived1>
auto operator*(VectorFieldBase<Derived0>&& u,
               const ScalarFieldBase<Derived1>& f) {
  return u * f;
}

template <typename Derived0, typename Derived1>
auto operator*(const VectorFieldBase<Derived0>& u,
               ScalarFieldBase<Derived1>&& f) {
  return u * f;
}

template <typename Derived0, typename Derived1>
auto operator*(VectorFieldBase<Derived0>&& u, ScalarFieldBase<Derived1>&& f) {
  return u * f;
}

template <typename Derived0, typename Derived1>
auto operator*(const ScalarFieldBase<Derived1>& f,
               const VectorFieldBase<Derived0>& u) {
  return u * f;
}

template <typename Derived0, typename Derived1>
auto operator*(ScalarFieldBase<Derived1>&& f,
               const VectorFieldBase<Derived0>& u) {
  return u * f;
}

template <typename Derived0, typename Derived1>
auto operator*(const ScalarFieldBase<Derived1>& f,
               VectorFieldBase<Derived0>&& u) {
  return u * f;
}

template <typename Derived0, typename Derived1>
auto operator*(ScalarFieldBase<Derived1>&& f, VectorFieldBase<Derived0>&& u) {
  return u * f;
}

// Scalar field division.
template <typename Derived0, typename Derived1>
auto operator/(const VectorFieldBase<Derived0>& u,
               const ScalarFieldBase<Derived1>& f) {
  return VectorFieldBinaryWithScalarField(u, f, std::divides<>());
}

template <typename Derived0, typename Derived1>
auto operator/(VectorFieldBase<Derived0>&& u,
               const ScalarFieldBase<Derived1>& f) {
  return u / f;
}

template <typename Derived0, typename Derived1>
auto operator/(const VectorFieldBase<Derived0>& u,
               ScalarFieldBase<Derived1>&& f) {
  return u / f;
}

template <typename Derived0, typename Derived1>
auto operator/(VectorFieldBase<Derived0>&& u, ScalarFieldBase<Derived1>&& f) {
  return u / f;
}

/*

// Addition.
template <typename Derived1, typename Derived2>
auto operator+(const VectorFieldBase<Derived1>& u1,
               const VectorFieldBase<Derived2>& u2) {
  return VectorFieldBinary(u1, u2, std::plus<>());
}

template <typename Derived1, typename Derived2>
auto operator+(const VectorFieldBase<Derived1>& u1,
               VectorFieldBase<Derived2>&& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(VectorFieldBase<Derived1>&& u1,
               const VectorFieldBase<Derived2>& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(VectorFieldBase<Derived1>&& u1, VectorFieldBase<Derived2>&& u2) {
  return u1 + u2;
}

// Subtraction.
template <typename Derived1, typename Derived2>
auto operator-(const VectorFieldBase<Derived1>& u1,
               const VectorFieldBase<Derived2>& u2) {
  return VectorFieldBinary(u1, u2, std::minus<>());
}

template <typename Derived1, typename Derived2>
auto operator-(const VectorFieldBase<Derived1>& u1,
               VectorFieldBase<Derived2>&& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(VectorFieldBase<Derived1>&& u1,
               const VectorFieldBase<Derived2>& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(VectorFieldBase<Derived1>&& u1, VectorFieldBase<Derived2>&& u2) {
  return u1 - u2;
}

// Multiplication.
template <typename Derived1, typename Derived2>
auto operator*(const VectorFieldBase<Derived1>& u1,
               const VectorFieldBase<Derived2>& u2) {
  return VectorFieldBinary(u1, u2, std::multiplies<>());
}

template <typename Derived1, typename Derived2>
auto operator*(const VectorFieldBase<Derived1>& u1,
               VectorFieldBase<Derived2>&& u2) {
  return u1 * u2;
}

template <typename Derived1, typename Derived2>
auto operator*(VectorFieldBase<Derived1>&& u1,
               const VectorFieldBase<Derived2>& u2) {
  return u1 * u2;
}

template <typename Derived1, typename Derived2>
auto operator*(VectorFieldBase<Derived1>&& u1, VectorFieldBase<Derived2>&& u2) {
  return u1 * u2;
}

// Division.
template <typename Derived1, typename Derived2>
auto operator/(const VectorFieldBase<Derived1>& u1,
               const VectorFieldBase<Derived2>& u2) {
  return VectorFieldBinary(u1, u2, std::divides<>());
}

template <typename Derived1, typename Derived2>
auto operator/(const VectorFieldBase<Derived1>& u1,
               VectorFieldBase<Derived2>&& u2) {
  return u1 / u2;
}

template <typename Derived1, typename Derived2>
auto operator/(VectorFieldBase<Derived1>&& u1,
               const VectorFieldBase<Derived2>& u2) {
  return u1 / u2;
}

template <typename Derived1, typename Derived2>
auto operator/(VectorFieldBase<Derived1>&& u1, VectorFieldBase<Derived2>&& u2) {
  return u1 / u2;
}

// Integation.
template <typename Derived>
auto Integrate(const VectorFieldBase<Derived>& u) {
  using Vector = typename Derived::Vector;
  auto w = u.Weights();
  auto i = 0;
  auto sum = Vector{0};
  for (auto [iTheta, iPhi] : u.PointIndices()) {
    sum += u(iTheta, iPhi) * w[i++];
  }
  return sum;
}

template <typename Derived>
auto Integrate(VectorFieldBase<Derived>&& u) {
  return Integrate(u);
}

// L2 inner product.
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(const VectorFieldBase<Derived1>& u1,
                    const VectorFieldBase<Derived2>& u2) {
  if constexpr (std::same_as<typename Derived1::Value, RealValued>) {
    return Integrate(u1 * u2);
  } else {
    return Integrate(conj(u1) * u2);
  }
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(VectorFieldBase<Derived1>&& u1,
                    const VectorFieldBase<Derived2>& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(const VectorFieldBase<Derived1>& u1,
                    VectorFieldBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(VectorFieldBase<Derived1>&& u1,
                    VectorFieldBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

// L2 norm.
template <typename Derived>
auto L2Norm(const VectorFieldBase<Derived>& u) {
  return std::sqrt(std::abs(L2InnerProduct(u, u)));
}

template <typename Derived>
auto L2Norm(VectorFieldBase<Derived>&& u) {
  return L2Norm(u);
}

// Equality operator (based on L2 norm).
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(const VectorFieldBase<Derived1>& u1,
                const VectorFieldBase<Derived2>& u2) {
  assert(u1.FieldSize() == u2.FieldSize());
  return L2Norm(u1 - u2) <
         std::numeric_limits<typename Derived1::Real>::epsilon();
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(VectorFieldBase<Derived1>&& u1,
                const VectorFieldBase<Derived2>& u2) {
  return u1 == u2;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(const VectorFieldBase<Derived1>& u1,
                VectorFieldBase<Derived2>&& u2) {
  return u1 == u2;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(VectorFieldBase<Derived1>&& u1,
                VectorFieldBase<Derived2>&& u2) {
  return u1 == u2;
}

*/

}  // namespace GSHTrans

#endif