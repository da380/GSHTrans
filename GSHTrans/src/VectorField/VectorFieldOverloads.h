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
#include "VectorFieldBinaryToScalarField.h"
#include "VectorFieldBinaryWithScalar.h"
#include "VectorFieldComplexToImag.h"
#include "VectorFieldComplexToReal.h"
#include "VectorFieldConjugate.h"
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

// Complex conjugation.
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(const VectorFieldBase<Derived>& u) {
  return VectorFieldConjugate(u);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(VectorFieldBase<Derived>&& u) {
  return conj(u);
}

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

// Pointwise inner product.
template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto InnerProduct(const VectorFieldBase<Derived0>& u0,
                  const VectorFieldBase<Derived1>& u1) {
  return VectorFieldBinaryToScalarField(
      u0, u1, [](auto m0, auto z0, auto p0, auto m1, auto z1, auto p1) {
        if constexpr (std::same_as<typename Derived0::Value, RealValued>) {
          return 2 * m0 * m1 + z0 * z1 + 2 * p0 * p1;
        } else {
          return std::conj(m0) * m1 + std::conj(z0) * z1 + std::conj(p0) * p1;
        }
      });
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto InnerProduct(VectorFieldBase<Derived0>&& u0,
                  const VectorFieldBase<Derived1>& u1) {
  return InnerProduct(u0, u1);
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto InnerProduct(const VectorFieldBase<Derived0>& u0,
                  VectorFieldBase<Derived1>&& u1) {
  return InnerProduct(u0, u1);
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto InnerProduct(VectorFieldBase<Derived0>&& u0,
                  VectorFieldBase<Derived1>&& u1) {
  return InnerProduct(u0, u1);
}

// Pointwise duality product.
template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto DualityProduct(const VectorFieldBase<Derived0>& u0,
                    const VectorFieldBase<Derived1>& u1) {
  return VectorFieldBinaryToScalarField(
      u0, u1, [](auto m0, auto z0, auto p0, auto m1, auto z1, auto p1) {
        if constexpr (std::same_as<typename Derived0::Value, RealValued>) {
          return -2 * m0 * m1 + z0 * z1 - 2 * p0 * p1;
        } else {
          return -m0 * p1 + z0 * z1 - p0 * m1;
        }
      });
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto DualityProduct(VectorFieldBase<Derived0>&& u0,
                    const VectorFieldBase<Derived1>& u1) {
  return DualityProduct(u0, u1);
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto DualityProduct(const VectorFieldBase<Derived0>& u0,
                    VectorFieldBase<Derived1>&& u1) {
  return DualityProduct(u0, u1);
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto DualityProduct(VectorFieldBase<Derived0>&& u0,
                    VectorFieldBase<Derived1>&& u1) {
  return DualityProduct(u0, u1);
}

// L2 inner product.
template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto L2InnerProduct(const VectorFieldBase<Derived0>& u0,
                    const VectorFieldBase<Derived1>& u1) {
  return Integrate(InnerProduct(u0, u1));
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto L2InnerProduct(VectorFieldBase<Derived0>&& u0,
                    const VectorFieldBase<Derived1>& u1) {
  return L2InnerProduct(u0, u1);
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto L2InnerProduct(const VectorFieldBase<Derived0>& u0,
                    VectorFieldBase<Derived1>&& u1) {
  return L2InnerProduct(u0, u1);
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto L2InnerProduct(VectorFieldBase<Derived0>&& u0,
                    VectorFieldBase<Derived1>&& u1) {
  return L2InnerProduct(u0, u1);
}

// L2 duality product.
template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto L2DualityProduct(const VectorFieldBase<Derived0>& u0,
                      const VectorFieldBase<Derived1>& u1) {
  return Integrate(DualityProduct(u0, u1));
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto L2DualityProduct(VectorFieldBase<Derived0>&& u0,
                      const VectorFieldBase<Derived1>& u1) {
  return L2DualityProduct(u0, u1);
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto L2DualityProduct(const VectorFieldBase<Derived0>& u0,
                      VectorFieldBase<Derived1>&& u1) {
  return L2DualityProduct(u0, u1);
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto L2DualityProduct(VectorFieldBase<Derived0>&& u0,
                      VectorFieldBase<Derived1>&& u1) {
  return L2DualityProduct(u0, u1);
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
template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto operator==(const VectorFieldBase<Derived0>& u0,
                const VectorFieldBase<Derived1>& u1) {
  assert(u0.FieldSize() == u1.FieldSize());
  auto scale = std::max(L2Norm(u0), L2Norm(u1));
  return L2Norm(u0 - u1) <
         std::numeric_limits<typename Derived1::Real>::epsilon() * scale;
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto operator==(VectorFieldBase<Derived0>&& u0,
                const VectorFieldBase<Derived1>& u1) {
  return u0 == u1;
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto operator==(const VectorFieldBase<Derived0>& u0,
                VectorFieldBase<Derived1>&& u1) {
  return u0 == u1;
}

template <typename Derived0, typename Derived1>
requires std::same_as<typename Derived0::Value, typename Derived1::Value>
auto operator==(VectorFieldBase<Derived0>&& u0,
                VectorFieldBase<Derived1>&& u1) {
  return u0 == u1;
}

}  // namespace GSHTrans

#endif