#ifndef GSH_TRANS_VECTOR_FIELD_OVERLOADS_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_OVERLOADS_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "../ScalarField/ScalarFieldBase.h"
#include "VectorFieldBase.h"

namespace GSHTrans {

//-----------------------------------------------------//
//              VectorField -> VectorField             //
//-----------------------------------------------------//
template <typename Derived>
auto operator-(const VectorFieldBase<Derived>& u) {
  return VectorFieldPointwiseUnary(u, [](auto x) { return -x; });
}

template <typename Derived>
auto operator-(VectorFieldBase<Derived>&& u) {
  return -u;
}

//-----------------------------------------------------//
//         ComplexVectorField -> RealVectorField       //
//-----------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(const VectorFieldBase<Derived>& u) {
  return RealifiedVectorField(u);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(VectorFieldBase<Derived>&& u) {
  return real(u);
}

//-----------------------------------------------------//
//         RealVectorField -> ComplexVectorField       //
//-----------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(const VectorFieldBase<Derived>& u) {
  return ComplexifiedVectorField(u);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(VectorFieldBase<Derived>&& u) {
  return complex(u);
}

//-----------------------------------------------------//
//       ComplexVectorField -> ComplexVectorField      //
//-----------------------------------------------------//
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

//-----------------------------------------------------//
//          VectorField x Scalar -> VectorField        //
//-----------------------------------------------------//
template <typename Derived>
auto operator*(const VectorFieldBase<Derived>& u, typename Derived::Scalar s) {
  return VectorFieldPointwiseUnaryWithScalar(u, std::multiplies<>(), s);
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

template <typename Derived>
auto operator/(const VectorFieldBase<Derived>& u, typename Derived::Scalar s) {
  return VectorFieldPointwiseUnaryWithScalar(u, std::divides<>(), s);
}

template <typename Derived>
auto operator/(VectorFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return u / s;
}

//-----------------------------------------------------//
//      VectorField x VectorField -> VectorField       //
//-----------------------------------------------------//

template <typename Derived1, typename Derived2>
auto operator+(const VectorFieldBase<Derived1>& u1,
               const VectorFieldBase<Derived2>& u2) {
  return VectorFieldPointwiseBinary(u1, u2, std::plus<>());
}

template <typename Derived1, typename Derived2>
auto operator+(VectorFieldBase<Derived1>&& u1,
               const VectorFieldBase<Derived2>& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(const VectorFieldBase<Derived1>& u1,
               VectorFieldBase<Derived2>&& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(VectorFieldBase<Derived1>&& u1, VectorFieldBase<Derived2>&& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator-(const VectorFieldBase<Derived1>& u1,
               const VectorFieldBase<Derived2>& u2) {
  return VectorFieldPointwiseBinary(u1, u2, std::minus<>());
}

template <typename Derived1, typename Derived2>
auto operator-(VectorFieldBase<Derived1>&& u1,
               const VectorFieldBase<Derived2>& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(const VectorFieldBase<Derived1>& u1,
               VectorFieldBase<Derived2>&& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(VectorFieldBase<Derived1>&& u1, VectorFieldBase<Derived2>&& u2) {
  return u1 - u2;
}

//-----------------------------------------------------//
//      VectorField x ScalarField -> VectorField       //
//-----------------------------------------------------//
template <typename Derived1, typename Derived2>
auto operator*(const VectorFieldBase<Derived1>& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return VectorFieldProductScalarField(u1, u2);
}

template <typename Derived1, typename Derived2>
auto operator*(VectorFieldBase<Derived1>&& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return u1 * u2;
}

template <typename Derived1, typename Derived2>
auto operator*(const VectorFieldBase<Derived1>& u1,
               ScalarFieldBase<Derived2>&& u2) {
  return u1 * u2;
}

template <typename Derived1, typename Derived2>
auto operator*(VectorFieldBase<Derived1>&& u1, ScalarFieldBase<Derived2>&& u2) {
  return u1 * u2;
}

template <typename Derived1, typename Derived2>
auto operator*(const ScalarFieldBase<Derived1>& u1,
               const VectorFieldBase<Derived2>& u2) {
  return u2 * u1;
}

template <typename Derived1, typename Derived2>
auto operator*(ScalarFieldBase<Derived1>&& u1,
               const VectorFieldBase<Derived2>& u2) {
  return u2 * u1;
}

template <typename Derived1, typename Derived2>
auto operator*(const ScalarFieldBase<Derived1>& u1,
               VectorFieldBase<Derived2>&& u2) {
  return u2 * u1;
}

template <typename Derived1, typename Derived2>
auto operator*(ScalarFieldBase<Derived1>&& u1, VectorFieldBase<Derived2>&& u2) {
  return u2 * u1;
}

//-----------------------------------------------------//
//      VectorField x VectorField -> ScalarField       //
//-----------------------------------------------------//
template <typename Derived1, typename Derived2>
auto InnerProduct(const VectorFieldBase<Derived1>& u1,
                  const VectorFieldBase<Derived2>& u2) {
  return VectorFieldInnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
auto InnerProduct(const VectorFieldBase<Derived1>& u1,
                  VectorFieldBase<Derived2>&& u2) {
  return InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
auto InnerProduct(VectorFieldBase<Derived1>&& u1,
                  const VectorFieldBase<Derived2>& u2) {
  return InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
auto InnerProduct(VectorFieldBase<Derived1>&& u1,
                  VectorFieldBase<Derived2>&& u2) {
  return InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
auto DualityProduct(const VectorFieldBase<Derived1>& u1,
                    const VectorFieldBase<Derived2>& u2) {
  return VectorFieldDualityProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
auto DualityProduct(const VectorFieldBase<Derived1>& u1,
                    VectorFieldBase<Derived2>&& u2) {
  return DualityProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
auto DualityProduct(VectorFieldBase<Derived1>&& u1,
                    const VectorFieldBase<Derived2>& u2) {
  return DualityProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
auto DualityProduct(VectorFieldBase<Derived1>&& u1,
                    VectorFieldBase<Derived2>&& u2) {
  return DualityProduct(u1, u2);
}

//-----------------------------------------------------//
//          VectorField x VectorField -> Scalar        //
//-----------------------------------------------------//
template <typename Derived1, typename Derived2>
auto L2InnerProduct(const VectorFieldBase<Derived1>& u1,
                    const VectorFieldBase<Derived2>& u2) {
  return Integrate(InnerProduct(u1, u2));
}

template <typename Derived1, typename Derived2>
auto L2InnerProduct(VectorFieldBase<Derived1>&& u1,
                    const VectorFieldBase<Derived2>& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
auto L2InnerProduct(const VectorFieldBase<Derived1>& u1,
                    VectorFieldBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
auto L2InnerProduct(VectorFieldBase<Derived1>&& u1,
                    VectorFieldBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

//-----------------------------------------------------//
//                VectorField  -> Real                 //
//-----------------------------------------------------//
template <typename Derived>
auto L2Norm(const VectorFieldBase<Derived>& u) {
  return std::sqrt(std::abs(L2InnerProduct(u, u)));
}

template <typename Derived>
auto L2Norm(VectorFieldBase<Derived>&& u) {
  return L2Norm(u);
}

//-------------------------------------------------------//
//           VectorField x VectorField -> bool           //
//-------------------------------------------------------//
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(const VectorFieldBase<Derived1>& u1,
                const VectorFieldBase<Derived2>& u2) {
  assert(u1.ComponentSize() == u2.ComponentSize());
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

}  // namespace GSHTrans

#endif  // GSH_TRANS_VECTOR_FIELD_OVERLOADS_GUARD_H