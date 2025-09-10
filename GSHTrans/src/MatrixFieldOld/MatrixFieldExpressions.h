#ifndef GSH_TRANS_MATRIX_FIELD_EXPRESSIONS_GUARD_H
#define GSH_TRANS_MATRIX_FIELD_EXPRESSIONS_GUARD_H

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

//-----------------------------------------------------//
//                 Negative expression                 //
//-----------------------------------------------------//
template <typename Derived>
class MatrixFieldNegative
    : public MatrixFieldBase<MatrixFieldNegative<Derived>> {
 public:
  using Int = typename MatrixFieldBase<MatrixFieldNegative<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using SelfAdjoint = typename Derived::SelfAdjoint;
  using Diagonal = typename Derived::Diagonal;
  using Isotropic = typename Derived::Isotropic;

  // Methods needed to inherit from MatrixFieldBase.
  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const -> Scalar {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha, beta);

    if constexpr (Diagonal::value) {
      return alpha + beta == 0 ? -_A(alpha, beta, iTheta, iPhi) : 0;
    } else {
      return -_A(-alpha, -beta, iTheta, iPhi);
    }
  }
  auto operator()(Int alpha, Int beta) const {
    this->CheckCanonicalIndices(alpha, beta);
    return MatrixFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  MatrixFieldNegative() = delete;
  MatrixFieldNegative(const MatrixFieldBase<Derived>& A) : _A{A} {}

  MatrixFieldNegative(const MatrixFieldNegative&) = default;
  MatrixFieldNegative(MatrixFieldNegative&&) = default;

  // Assignment.
  MatrixFieldNegative& operator=(MatrixFieldNegative&) = default;
  MatrixFieldNegative& operator=(MatrixFieldNegative&&) = default;

 private:
  const MatrixFieldBase<Derived>& _A;
};

//-----------------------------------------------------//
//                  Adjoint expression                 //
//-----------------------------------------------------//
template <typename Derived>
class MatrixFieldAdjoint : public MatrixFieldBase<MatrixFieldAdjoint<Derived>> {
 public:
  using Int = typename MatrixFieldBase<MatrixFieldAdjoint<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using SelfAdjoint = typename Derived::SelfAdjoint;
  using Diagonal = typename Derived::Diagonal;
  using Isotropic = typename Derived::Isotropic;

  // Methods needed to inherit from MatrixFieldBase.
  auto GetGrid() const { return _A.GetGrid(); }

  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const -> Scalar
  requires std::same_as<Value, RealValued>
  {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha, beta);
    if constexpr (Diagonal::value) {
      return alpha + beta == 0 ? _A(alpha, beta, iTheta, iPhi) : 0;
    } else {
      return _A(beta, alpha, iTheta, iPhi);
    }
  }

  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const -> Scalar
  requires std::same_as<Value, ComplexValued>
  {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha, beta);
    if constexpr (Diagonal::value) {
      if constexpr (SelfAdjoint::value) {
        return alpha + beta == 0 ? _A(alpha, beta, iTheta, iPhi) : 0;
      } else {
        return alpha + beta == 0 ? std::conj(_A(alpha, beta, iTheta, iPhi)) : 0;
      }
    }
    if constexpr (SelfAdjoint::value) {
      return _A(alpha, beta, iTheta, iPhi);
    } else {
      return std::conj(_A(beta, alpha, iTheta, iPhi));
    }
  }

  auto operator()(Int alpha, Int beta) const {
    this->CheckCanonicalIndices(alpha, beta);
    return MatrixFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  MatrixFieldAdjoint() = delete;
  MatrixFieldAdjoint(const MatrixFieldBase<Derived>& A) : _A{A} {}

  MatrixFieldAdjoint(const MatrixFieldAdjoint&) = default;
  MatrixFieldAdjoint(MatrixFieldAdjoint&&) = default;

  // Assignment.
  MatrixFieldAdjoint& operator=(MatrixFieldAdjoint&) = default;
  MatrixFieldAdjoint& operator=(MatrixFieldAdjoint&&) = default;

 private:
  const MatrixFieldBase<Derived>& _A;
};

//-------------------------------------------------//
//     Complexification of real matrix field       //
//-------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
class ComplexifiedMatrixField
    : public MatrixFieldBase<ComplexifiedMatrixField<Derived>> {
 public:
  using Int = typename MatrixFieldBase<ComplexifiedMatrixField<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Value = ComplexValued;
  using Scalar = Complex;
  using SelfAdjoint = typename Derived::SelfAdjoint;
  using Diagonal = typename Derived::Diagonal;
  using Isotropic = typename Derived::Isotropic;

  // Methods needed to inherit from MatrixField Base.
  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const -> Complex {
    this->CheckCanonicalIndices(alpha, beta);
    this->CheckPointIndices(iTheta, iPhi);
    constexpr auto ii = Complex(0, 1);
    auto index = 3 * alpha + beta;
    if (index < 0) {
      return static_cast<Scalar>(MinusOneToPower(alpha + beta)) *
             (_A(-alpha, -beta, iTheta, iPhi) -
              ii * _A(alpha, beta, iTheta, iPhi));
    } else if (index == 0) {
      return _A(0, 0, iTheta, iPhi);
    } else {
      return _A(alpha, beta, iTheta, iPhi) +
             ii * _A(-alpha, -beta, iTheta, iPhi);
    }
  }

  auto operator()(Int alpha, Int beta) const {
    this->CheckCanonicalIndices(alpha, beta);
    return MatrixFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  ComplexifiedMatrixField(const MatrixFieldBase<Derived>& A) : _A{A} {}

  ComplexifiedMatrixField(const ComplexifiedMatrixField&) = default;
  ComplexifiedMatrixField(ComplexifiedMatrixField&&) = default;

  // Assignment.
  ComplexifiedMatrixField& operator=(const ComplexifiedMatrixField&) = default;
  ComplexifiedMatrixField& operator=(ComplexifiedMatrixField&&) = default;

 private:
  const MatrixFieldBase<Derived>& _A;
};

//-------------------------------------------------//
//     Realification of complex matrix field       //
//-------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
class RealifiedMatrixField
    : public MatrixFieldBase<RealifiedMatrixField<Derived>> {
 public:
  using Int = typename MatrixFieldBase<RealifiedMatrixField<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Value = RealValued;
  using Scalar = Real;
  using SelfAdjoint = typename Derived::SelfAdjoint;
  using Diagonal = typename Derived::Diagonal;
  using Isotropic = typename Derived::Isotropic;

  // Methods needed to inherit from MatrixField Base.
  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    this->CheckCanonicalIndices(alpha, beta);
    this->CheckPointIndices(iTheta, iPhi);
    constexpr auto half = static_cast<Real>(1) / static_cast<Real>(2);
    auto index = 3 * alpha + beta;
    auto fac = static_cast<Scalar>(MinusOneToPower(alpha + beta));
    if (index < 0) {
      return half * std::imag(_A(-alpha, -beta, iTheta, iPhi) -
                              fac * _A(alpha, beta, iTheta, iPhi));
    } else if (index == 0) {
      return std::real(_A(0, 0, iTheta, iPhi));
    } else {
      return half * std::real(_A(alpha, beta, iTheta, iPhi) +
                              fac * _A(-alpha, -beta, iTheta, iPhi));
    }
  }

  auto operator()(Int alpha, Int beta) const {
    this->CheckCanonicalIndices(alpha, beta);
    return MatrixFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  RealifiedMatrixField(const MatrixFieldBase<Derived>& A) : _A{A} {}

  RealifiedMatrixField(const RealifiedMatrixField&) = default;
  RealifiedMatrixField(RealifiedMatrixField&&) = default;

  // Assignment.
  RealifiedMatrixField& operator=(const RealifiedMatrixField&) = default;
  RealifiedMatrixField& operator=(RealifiedMatrixField&&) = default;

 private:
  const MatrixFieldBase<Derived>& _A;
};

//-------------------------------------------------//
//            Pointwise Binary expression          //
//-------------------------------------------------//
template <typename Derived1, typename Derived2, typename Function>
requires requires() {
  requires std::same_as<typename Derived1::Real, typename Derived2::Real>;
  requires std::same_as<typename Derived1::Value, typename Derived2::Value>;
  requires std::invocable<Function, typename Derived1::Scalar,
                          typename Derived2::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived1::Scalar,
                           typename Derived1::Scalar>,
      typename Derived1::Scalar>;
}
class MatrixFieldBinary
    : public MatrixFieldBase<MatrixFieldBinary<Derived1, Derived2, Function>> {
 public:
  using Int = typename MatrixFieldBase<
      MatrixFieldBinary<Derived1, Derived2, Function>>::Int;
  using Grid = typename Derived1::Grid;
  using Scalar = typename Derived1::Scalar;
  using Value = typename Derived1::Value;
  using Real = typename Derived1::Real;
  using Complex = typename Derived1::Complex;

  // Methods needed to inherit from MatrixFieldBase.
  auto GetGrid() const { return _A1.GetGrid(); }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha, beta);
    return _f(_A1(alpha, beta, iTheta, iPhi), _A2(alpha, beta, iTheta, iPhi));
  }
  auto operator()(Int alpha, Int beta) const {
    this->CheckCanonicalIndices(alpha, beta);
    return MatrixFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  MatrixFieldBinary() = delete;
  MatrixFieldBinary(const MatrixFieldBase<Derived1>& A1,
                    const MatrixFieldBase<Derived2>& A2, Function f)
      : _A1{A1}, _A2{A2}, _f{f} {
    assert(_A1.FieldSize() == _A2.FieldSize());
  }

  MatrixFieldBinary(const MatrixFieldBinary&) = default;
  MatrixFieldBinary(MatrixFieldBinary&&) = default;

  // Assignment.
  MatrixFieldBinary& operator=(MatrixFieldBinary&) = default;
  MatrixFieldBinary& operator=(MatrixFieldBinary&&) = default;

 private:
  const MatrixFieldBase<Derived1>& _A1;
  const MatrixFieldBase<Derived2>& _A2;
  Function _f;
};

//-------------------------------------------------//
//             Matrix trace expression             //
//-------------------------------------------------//
template <typename Derived>
class MatrixFieldTrace : public ScalarFieldBase<MatrixFieldTrace<Derived>> {
 public:
  using Int = typename ScalarFieldBase<MatrixFieldTrace<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);

    if constexpr (std::same_as<Value, ComplexValued>) {
      return -_A(-1, 1, iTheta, iPhi) + _A(0, 0, iTheta, iPhi) -
             _A(+1, 1, iTheta, iPhi);
    } else {
      return -2 * _A(1, -1, iTheta, iPhi) + _A(0, 0, iTheta, iPhi);
    }
  }

  // Constructors.
  MatrixFieldTrace() = delete;
  MatrixFieldTrace(const MatrixFieldBase<Derived>& A) : _A{A} {}

  MatrixFieldTrace(const MatrixFieldTrace&) = default;
  MatrixFieldTrace(MatrixFieldTrace&&) = default;

  // Assignment.
  MatrixFieldTrace& operator=(MatrixFieldTrace&) = default;
  MatrixFieldTrace& operator=(MatrixFieldTrace&&) = default;

 private:
  const MatrixFieldBase<Derived>& _A;
};

//-------------------------------------------------//
//          Matrix determinant expression          //
//-------------------------------------------------//
template <typename Derived>
class MatrixFieldDeterminant
    : public ScalarFieldBase<MatrixFieldDeterminant<Derived>> {
 public:
  using Int = typename ScalarFieldBase<MatrixFieldDeterminant<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);

    if constexpr (std::same_as<Value, ComplexValued>) {
      return -_A(-1, -1, iTheta, iPhi) *
                 (_A(0, 0, iTheta, iPhi) * _A(+1, +1, iTheta, iPhi) -
                  _A(0, +1, iTheta, iPhi) * _A(+1, 0, iTheta, iPhi)) +
             _A(-1, 0, iTheta, iPhi) *
                 (_A(0, -1, iTheta, iPhi) * _A(+1, +1, iTheta, iPhi) -
                  _A(0, +1, iTheta, iPhi) * _A(+1, -1, iTheta, iPhi)) -
             _A(-1, +1, iTheta, iPhi) *
                 (_A(0, -1, iTheta, iPhi) * _A(+1, 0, iTheta, iPhi) -
                  _A(0, 0, iTheta, iPhi) * _A(+1, -1, iTheta, iPhi));
    } else {
      return _A(0, 0, iTheta, iPhi) *
                 (_A(-1, +1, iTheta, iPhi) * _A(-1, +1, iTheta, iPhi) +
                  _A(+1, -1, iTheta, iPhi) * _A(+1, -1, iTheta, iPhi) -
                  _A(-1, -1, iTheta, iPhi) * _A(-1, -1, iTheta, iPhi) -
                  _A(+1, +1, iTheta, iPhi) * _A(+1, +1, iTheta, iPhi)) +
             2 * _A(+1, -1, iTheta, iPhi) *
                 (_A(-1, 0, iTheta, iPhi) * _A(0, -1, iTheta, iPhi) +
                  _A(+1, 0, iTheta, iPhi) * _A(0, +1, iTheta, iPhi)) +
             2 * _A(-1, +1, iTheta, iPhi) *
                 (_A(-1, 0, iTheta, iPhi) * _A(0, +1, iTheta, iPhi) -
                  _A(+1, 0, iTheta, iPhi) * _A(0, -1, iTheta, iPhi)) +
             2 * _A(-1, -1, iTheta, iPhi) *
                 (_A(-1, 0, iTheta, iPhi) * _A(0, +1, iTheta, iPhi) +
                  _A(+1, 0, iTheta, iPhi) * _A(0, -1, iTheta, iPhi)) +
             2 * _A(+1, +1, iTheta, iPhi) *
                 (_A(+1, 0, iTheta, iPhi) * _A(0, +1, iTheta, iPhi) -
                  _A(-1, 0, iTheta, iPhi) * _A(0, -1, iTheta, iPhi));
    }
  }

  // Constructors.
  MatrixFieldDeterminant() = delete;
  MatrixFieldDeterminant(const MatrixFieldBase<Derived>& A) : _A{A} {}

  MatrixFieldDeterminant(const MatrixFieldDeterminant&) = default;
  MatrixFieldDeterminant(MatrixFieldDeterminant&&) = default;

  // Assignment.
  MatrixFieldDeterminant& operator=(MatrixFieldDeterminant&) = default;
  MatrixFieldDeterminant& operator=(MatrixFieldDeterminant&&) = default;

 private:
  const MatrixFieldBase<Derived>& _A;
};

//-------------------------------------------------//
//             Matrix action expression            //
//-------------------------------------------------//
template <typename Derived1, typename Derived2>
requires requires() {
  requires std::same_as<typename Derived1::Real, typename Derived2::Real>;
  requires std::same_as<typename Derived1::Value, typename Derived2::Value>;
}
class MatrixFieldAction
    : public VectorFieldBase<MatrixFieldAction<Derived1, Derived2>> {
 public:
  using Int =
      typename VectorFieldBase<MatrixFieldAction<Derived1, Derived2>>::Int;
  using Grid = typename Derived1::Grid;
  using Scalar = typename Derived1::Scalar;
  using Value = typename Derived1::Value;
  using Real = typename Derived1::Real;
  using Complex = typename Derived1::Complex;

  // Methods needed to inherit from VectorFieldBase.
  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    if constexpr (std::same_as<Value, Complex>) {
      return -_A(alpha, -1, iTheta, iPhi) * _u(+1, iTheta, iPhi) +
             _A(alpha, 0, iTheta, iPhi) * _u(0, iTheta, iPhi) -
             _A(alpha, 1, iTheta, iPhi) * _u(-1, iTheta, iPhi);
    } else {
      switch (alpha) {
        case -1: {
          return -(_A(+1, -1, iTheta, iPhi) + _A(+1, +1, iTheta, iPhi)) *
                     _u(-1, iTheta, iPhi) +
                 _A(-1, 0, iTheta, iPhi) * _u(0, iTheta, iPhi) +
                 (_A(-1, -1, iTheta, iPhi) - _A(-1, +1, iTheta, iPhi)) *
                     _u(+1, iTheta, iPhi);
        }
        case 0: {
          return 2 * _A(0, -1, iTheta, iPhi) * _u(-1, iTheta, iPhi) +
                 _A(0, 0, iTheta, iPhi) * _u(0, iTheta, iPhi) +
                 2 * _A(0, +1, iTheta, iPhi) * _u(+1, iTheta, iPhi);
        }
        case +1: {
          return (_A(-1, -1, iTheta, iPhi) + _A(-1, +1, iTheta, iPhi)) *
                     _u(-1, iTheta, iPhi) +
                 _A(+1, 0, iTheta, iPhi) * _u(0, iTheta, iPhi) +
                 (_A(+1, +1, iTheta, iPhi) - _A(+1, -1, iTheta, iPhi)) *
                     _u(+1, iTheta, iPhi);
        }
        default:
          return 0;
      }
    }
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  MatrixFieldAction() = delete;
  MatrixFieldAction(const MatrixFieldBase<Derived1>& A,
                    const VectorFieldBase<Derived2>& u)
      : _A{A}, _u{u} {
    assert(_A.FieldSize() == _u.FieldSize());
  }

  MatrixFieldAction(const MatrixFieldAction&) = default;
  MatrixFieldAction(MatrixFieldAction&&) = default;

  // Assignment.
  MatrixFieldAction& operator=(MatrixFieldAction&) = default;
  MatrixFieldAction& operator=(MatrixFieldAction&&) = default;

 private:
  const MatrixFieldBase<Derived1>& _A;
  const VectorFieldBase<Derived2>& _u;
};

//-------------------------------------------------//
//         Matrix multiplication expression        //
//-------------------------------------------------//
template <typename Derived1, typename Derived2>
requires requires() {
  requires std::same_as<typename Derived1::Real, typename Derived2::Real>;
  requires std::same_as<typename Derived1::Value, typename Derived2::Value>;
}
class MatrixFieldMultiplication
    : public MatrixFieldBase<MatrixFieldMultiplication<Derived1, Derived2>> {
 public:
  using Int = typename VectorFieldBase<
      MatrixFieldMultiplication<Derived1, Derived2>>::Int;
  using Grid = typename Derived1::Grid;
  using Scalar = typename Derived1::Scalar;
  using Value = typename Derived1::Value;
  using Real = typename Derived1::Real;
  using Complex = typename Derived1::Complex;

  // Methods needed to inherit from VectorFieldBase.
  auto GetGrid() const { return _A1.GetGrid(); }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha, beta);
    if constexpr (std::same_as<Value, Complex>) {
      return -_A1(alpha, +1, iTheta, iPhi) * _A2(-1, beta, iTheta, iPhi) +
             _A1(alpha, 0, iTheta, iPhi) * _A2(0, beta, iTheta, iPhi) -
             _A1(alpha, -1, iTheta, iPhi) * _A2(+1, beta, iTheta, iPhi);
    } else {
      return 0;
    }
  }
  auto operator()(Int alpha, Int beta) const {
    this->CheckCanonicalIndices(alpha, beta);
    return VectorFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  MatrixFieldMultiplication() = delete;
  MatrixFieldMultiplication(const MatrixFieldBase<Derived1>& A1,
                            const MatrixFieldBase<Derived2>& A2)
      : _A1{A1}, _A2{A2} {
    assert(_A1.FieldSize() == _A2.FieldSize());
  }

  MatrixFieldMultiplication(const MatrixFieldMultiplication&) = default;
  MatrixFieldMultiplication(MatrixFieldMultiplication&&) = default;

  // Assignment.
  MatrixFieldMultiplication& operator=(MatrixFieldMultiplication&) = default;
  MatrixFieldMultiplication& operator=(MatrixFieldMultiplication&&) = default;

 private:
  const MatrixFieldBase<Derived1>& _A1;
  const MatrixFieldBase<Derived2>& _A2;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_MATRIX_FIELD_EXPRESSIONS_GUARD_H
