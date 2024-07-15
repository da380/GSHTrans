#ifndef GSH_TRANS_MATRIX_FIELD_OVERLOADS_GUARD_H
#define GSH_TRANS_MATRIX_FIELD_OVERLOADS_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "../ScalarField/ScalarFieldBase.h"
#include "../VectorField/VectorFieldBase.h"
#include "MatrixFieldBase.h"

namespace GSHTrans {

//  Negative operator.
template <typename Derived>
auto operator-(const MatrixFieldBase<Derived>& A) {
  return MatrixFieldNegative(A);
}

template <typename Derived>
auto operator-(MatrixFieldBase<Derived>&& A) {
  return -A;
}

// Adjoint operator.
template <typename Derived>
auto Adjoint(const MatrixFieldBase<Derived>& A) {
  return MatrixFieldAdjoint(A);
}

template <typename Derived>
auto Adjoint(MatrixFieldBase<Derived>&& A) {
  return Adjoint(A);
}

/*

//-----------------------------------------------------//
//              MatrixField -> MatrixField             //
//-----------------------------------------------------//
template <typename Derived>
auto operator-(const MatrixFieldBase<Derived>& A) {
  return MatrixFieldPointwiseUnary(A, [](auto x) { return -x; });
}

template <typename Derived>
auto operator-(MatrixFieldBase<Derived>&& A) {
  return -A;
}

template <typename Derived>
auto Transpose(const MatrixFieldBase<Derived>& A) {
  return MatrixFieldTranspose(A);
}

template <typename Derived>
auto Transpose(MatrixFieldBase<Derived>&& A) {
  return Transpose(A);
}

template <typename Derived>
auto Adjoint(const MatrixFieldBase<Derived>& A) {
  return MatrixFieldAdjoint(A);
}

template <typename Derived>
auto Adjoint(MatrixFieldBase<Derived>&& A) {
  return Adjoint(A);
}

//-----------------------------------------------------//
//         RealMatrixField -> ComplexMatrixField       //
//-----------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(const MatrixFieldBase<Derived>& A) {
  return ComplexifiedMatrixField(A);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(MatrixFieldBase<Derived>&& A) {
  return complex(A);
}

//-----------------------------------------------------//
//         ComplexMatrixField -> RealMatrixField       //
//-----------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(const MatrixFieldBase<Derived>& A) {
  return RealifiedMatrixField(A);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(MatrixFieldBase<Derived>&& A) {
  return real(A);
}

//-----------------------------------------------------//
//       ComplexMatrixField -> ComplexMatrixField      //
//-----------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(const MatrixFieldBase<Derived>& A) {
  return MatrixFieldConjugate(A);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(MatrixFieldBase<Derived>&& A) {
  return conj(A);
}

//-----------------------------------------------------//
//               MatrixField -> ScalarField            //
//-----------------------------------------------------//
template <typename Derived>
auto Tr(const MatrixFieldBase<Derived>& A) {
  return MatrixFieldTrace(A);
}

template <typename Derived>
auto Tr(MatrixFieldBase<Derived>&& A) {
  return Tr(A);
}

template <typename Derived>
auto Det(const MatrixFieldBase<Derived>& A) {
  return MatrixFieldDeterminant(A);
}

template <typename Derived>
auto Det(MatrixFieldBase<Derived>&& A) {
  return Det(A);
}

//-----------------------------------------------------//
//          MatrixField x Scalar -> MatrixField        //
//-----------------------------------------------------//
template <typename Derived>
auto operator*(const MatrixFieldBase<Derived>& A, typename Derived::Scalar s) {
  return MatrixFieldPointwiseUnaryWithScalar(A, std::multiplies<>(), s);
}

template <typename Derived>
auto operator*(typename Derived::Scalar s, const MatrixFieldBase<Derived>& A) {
  return A * s;
}

template <typename Derived>
auto operator*(MatrixFieldBase<Derived>&& A, typename Derived::Scalar s) {
  return A * s;
}

template <typename Derived>
auto operator*(typename Derived::Scalar s, MatrixFieldBase<Derived>&& A) {
  return A * s;
}

template <typename Derived>
auto operator/(const MatrixFieldBase<Derived>& A, typename Derived::Scalar s) {
  return MatrixFieldPointwiseUnaryWithScalar(A, std::divides<>(), s);
}

template <typename Derived>
auto operator/(MatrixFieldBase<Derived>&& A, typename Derived::Scalar s) {
  return A / s;
}

//-----------------------------------------------------//
//      MatrixField x MatrixField -> MatrixField       //
//-----------------------------------------------------//

template <typename Derived1, typename Derived2>
auto operator+(const MatrixFieldBase<Derived1>& A1,
               const MatrixFieldBase<Derived2>& A2) {
  return MatrixFieldBinary(A1, A2, std::plus<>());
}

template <typename Derived1, typename Derived2>
auto operator+(MatrixFieldBase<Derived1>&& A1,
               const MatrixFieldBase<Derived2>& A2) {
  return A1 + A2;
}

template <typename Derived1, typename Derived2>
auto operator+(const MatrixFieldBase<Derived1>& A1,
               MatrixFieldBase<Derived2>&& A2) {
  return A1 + A2;
}

template <typename Derived1, typename Derived2>
auto operator+(MatrixFieldBase<Derived1>&& A1, MatrixFieldBase<Derived2>&& A2) {
  return A1 + A2;
}

template <typename Derived1, typename Derived2>
auto operator-(const MatrixFieldBase<Derived1>& A1,
               const MatrixFieldBase<Derived2>& A2) {
  return MatrixFieldBinary(A1, A2, std::minus<>());
}

template <typename Derived1, typename Derived2>
auto operator-(MatrixFieldBase<Derived1>&& A1,
               const MatrixFieldBase<Derived2>& A2) {
  return A1 - A2;
}

template <typename Derived1, typename Derived2>
auto operator-(const MatrixFieldBase<Derived1>& A1,
               MatrixFieldBase<Derived2>&& A2) {
  return A1 - A2;
}

template <typename Derived1, typename Derived2>
auto operator-(MatrixFieldBase<Derived1>&& A1, MatrixFieldBase<Derived2>&& A2) {
  return A1 - A2;
}

*/

}  // namespace GSHTrans

#endif  // GSH_TRANS_MATRIX_FIELD_OVERLOADS_GUARD_H