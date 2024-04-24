#ifndef GSH_TRANS_VECTOR_FIELD_EXPRESSIONS_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_EXPRESSIONS_GUARD_H

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
//            Complex conjugate expression             //
//-----------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
class VectorFieldConjugate
    : public VectorFieldBase<VectorFieldConjugate<Derived>> {
 public:
  using Int = typename VectorFieldBase<VectorFieldConjugate<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from VectorFieldBase.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return alpha == 0 ? std::conj(_u(0, iTheta, iPhi))
                      : -std::conj(_u(-alpha, iTheta, iPhi));
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  VectorFieldConjugate() = delete;
  VectorFieldConjugate(const VectorFieldBase<Derived>& u) : _u{u} {}

  VectorFieldConjugate(const VectorFieldConjugate&) = default;
  VectorFieldConjugate(VectorFieldConjugate&&) = default;

  // Assignment.
  VectorFieldConjugate& operator=(VectorFieldConjugate&) = default;
  VectorFieldConjugate& operator=(VectorFieldConjugate&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
};

//-------------------------------------------------//
//     Complexification of real vector field       //
//-------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
class ComplexifiedVectorField
    : public VectorFieldBase<ComplexifiedVectorField<Derived>> {
 public:
  using Int = typename VectorFieldBase<ComplexifiedVectorField<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Value = ComplexValued;
  using Scalar = Complex;

  // Methods needed to inherit from VectorField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const -> Complex {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    constexpr auto ii = Complex(0, 1);
    switch (alpha) {
      case -1:
        return -_u(1, iTheta, iPhi) + ii * _u(-1, iTheta, iPhi);
      case 0:
        return _u(0, iTheta, iPhi);
      case 1:
        return _u(1, iTheta, iPhi) + ii * _u(-1, iTheta, iPhi);
      default:
        return 0;
    }
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  ComplexifiedVectorField(const VectorFieldBase<Derived>& u) : _u{u} {}

  ComplexifiedVectorField(const ComplexifiedVectorField&) = default;
  ComplexifiedVectorField(ComplexifiedVectorField&&) = default;

  // Assignment.
  ComplexifiedVectorField& operator=(const ComplexifiedVectorField&) = default;
  ComplexifiedVectorField& operator=(ComplexifiedVectorField&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
};

//-------------------------------------------------//
//      Realification of complex vector field      //
//-------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
class RealifiedVectorField
    : public VectorFieldBase<RealifiedVectorField<Derived>> {
 public:
  using Int = typename VectorFieldBase<RealifiedVectorField<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Value = RealValued;
  using Scalar = Real;

  // Methods needed to inherit from VectorField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const -> Real {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    constexpr auto half = static_cast<Real>(1) / static_cast<Real>(2);
    switch (alpha) {
      case -1:
        return half * std::imag(_u(1, iTheta, iPhi) + _u(-1, iTheta, iPhi));
      case 0:
        return std::real(_u(0, iTheta, iPhi));
      case 1:
        return half * std::real(_u(1, iTheta, iPhi) - _u(-1, iTheta, iPhi));
      default:
        return 0;
    }
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  RealifiedVectorField(const VectorFieldBase<Derived>& u) : _u{u} {}

  RealifiedVectorField(const RealifiedVectorField&) = default;
  RealifiedVectorField(RealifiedVectorField&&) = default;

  // Assignment.
  RealifiedVectorField& operator=(const RealifiedVectorField&) = default;
  RealifiedVectorField& operator=(RealifiedVectorField&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
};

//-------------------------------------------------//
//              PointwiseUnary expression          //
//-------------------------------------------------//
template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class VectorFieldPointwiseUnary
    : public VectorFieldBase<VectorFieldPointwiseUnary<Derived, Function>> {
 public:
  using Int = typename VectorFieldBase<
      VectorFieldPointwiseUnary<Derived, Function>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from VectorField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return _f(_u(alpha, iTheta, iPhi));
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  VectorFieldPointwiseUnary() = delete;
  VectorFieldPointwiseUnary(const VectorFieldBase<Derived>& u, Function f)
      : _u{u}, _f{f} {}

  VectorFieldPointwiseUnary(const VectorFieldPointwiseUnary&) = default;
  VectorFieldPointwiseUnary(VectorFieldPointwiseUnary&&) = default;

  // Assignment.
  VectorFieldPointwiseUnary& operator=(VectorFieldPointwiseUnary&) = default;
  VectorFieldPointwiseUnary& operator=(VectorFieldPointwiseUnary&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
  Function _f;
};

//-------------------------------------------------//
//       PointwiseUnary expression with Scalar     //
//-------------------------------------------------//
template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar,
                          typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar,
                           typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class VectorFieldPointwiseUnaryWithScalar
    : public VectorFieldBase<
          VectorFieldPointwiseUnaryWithScalar<Derived, Function>> {
 public:
  using Int = typename VectorFieldBase<
      VectorFieldPointwiseUnaryWithScalar<Derived, Function>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from VectorFieldBase.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return _f(_u(alpha, iTheta, iPhi), _s);
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  VectorFieldPointwiseUnaryWithScalar() = delete;
  VectorFieldPointwiseUnaryWithScalar(const VectorFieldBase<Derived>& u,
                                      Function f, Scalar s)
      : _u{u}, _f{f}, _s{s} {}

  VectorFieldPointwiseUnaryWithScalar(
      const VectorFieldPointwiseUnaryWithScalar&) = default;
  VectorFieldPointwiseUnaryWithScalar(VectorFieldPointwiseUnaryWithScalar&&) =
      default;

  // Assignment.
  VectorFieldPointwiseUnaryWithScalar& operator=(
      VectorFieldPointwiseUnaryWithScalar&) = default;
  VectorFieldPointwiseUnaryWithScalar& operator=(
      VectorFieldPointwiseUnaryWithScalar&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
  Function _f;
  Scalar _s;
};

//-------------------------------------------------//
//            PointwiseBinary expression           //
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
class VectorFieldPointwiseBinary
    : public VectorFieldBase<
          VectorFieldPointwiseBinary<Derived1, Derived2, Function>> {
 public:
  using Int = typename VectorFieldBase<
      VectorFieldPointwiseBinary<Derived1, Derived2, Function>>::Int;
  using Grid = typename Derived1::Grid;
  using Scalar = typename Derived1::Scalar;
  using Value = typename Derived1::Value;
  using Real = typename Derived1::Real;
  using Complex = typename Derived1::Complex;

  // Methods needed to inherit from VectorFieldBase.
  auto GetGrid() const { return _u1.GetGrid(); }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return _f(_u1(alpha, iTheta, iPhi), _u2(alpha, iTheta, iPhi));
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  VectorFieldPointwiseBinary() = delete;
  VectorFieldPointwiseBinary(const VectorFieldBase<Derived1>& u1,
                             const VectorFieldBase<Derived2>& u2, Function f)
      : _u1{u1}, _u2{u2}, _f{f} {
    assert(_u1.FieldSize() == _u2.FieldSize());
  }

  VectorFieldPointwiseBinary(const VectorFieldPointwiseBinary&) = default;
  VectorFieldPointwiseBinary(VectorFieldPointwiseBinary&&) = default;

  // Assignment.
  VectorFieldPointwiseBinary& operator=(VectorFieldPointwiseBinary&) = default;
  VectorFieldPointwiseBinary& operator=(VectorFieldPointwiseBinary&&) = default;

 private:
  const VectorFieldBase<Derived1>& _u1;
  const VectorFieldBase<Derived2>& _u2;
  Function _f;
};

//-----------------------------------------------------//
//      Pointwise inner/duality product expressions    //
//-----------------------------------------------------//
template <typename Derived1, typename Derived2>
requires requires() {
  requires std::same_as<typename Derived1::Real, typename Derived2::Real>;
  requires std::same_as<typename Derived1::Value, typename Derived2::Value>;
}
class VectorFieldInnerProduct
    : public ScalarFieldBase<VectorFieldInnerProduct<Derived1, Derived2>> {
 public:
  using Int = typename ScalarFieldBase<
      VectorFieldInnerProduct<Derived1, Derived2>>::Int;
  using Grid = typename Derived1::Grid;
  using Scalar = typename Derived1::Scalar;
  using Value = typename Derived1::Value;
  using Real = typename Derived1::Real;
  using Complex = typename Derived1::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u1.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    if constexpr (std::same_as<Value, RealValued>) {
      return 2 * _u1(-1, iTheta, iPhi) * _u2(-1, iTheta, iPhi) +
             _u1(0, iTheta, iPhi) * _u2(0, iTheta, iPhi) +
             2 * _u1(1, iTheta, iPhi) * _u2(1, iTheta, iPhi);
    } else {
      return std::conj(_u1(-1, iTheta, iPhi)) * _u2(-1, iTheta, iPhi) +
             std::conj(_u1(0, iTheta, iPhi)) * _u2(0, iTheta, iPhi) +
             std::conj(_u1(1, iTheta, iPhi)) * _u2(1, iTheta, iPhi);
    }
  }

  // Constructors.
  VectorFieldInnerProduct() = delete;
  VectorFieldInnerProduct(const VectorFieldBase<Derived1>& u1,
                          const VectorFieldBase<Derived2>& u2)
      : _u1{u1}, _u2{u2} {
    assert(_u1.FieldSize() == _u2.FieldSize());
  }

  VectorFieldInnerProduct(const VectorFieldInnerProduct&) = default;
  VectorFieldInnerProduct(VectorFieldInnerProduct&&) = default;

  // Assignment.
  VectorFieldInnerProduct& operator=(VectorFieldInnerProduct&) = default;
  VectorFieldInnerProduct& operator=(VectorFieldInnerProduct&&) = default;

 private:
  const VectorFieldBase<Derived1>& _u1;
  const VectorFieldBase<Derived2>& _u2;
};

template <typename Derived1, typename Derived2>
requires requires() {
  requires std::same_as<typename Derived1::Real, typename Derived2::Real>;
  requires std::same_as<typename Derived1::Value, typename Derived2::Value>;
}
class VectorFieldDualityProduct
    : public ScalarFieldBase<VectorFieldInnerProduct<Derived1, Derived2>> {
 public:
  using Int = typename ScalarFieldBase<
      VectorFieldInnerProduct<Derived1, Derived2>>::Int;
  using Grid = typename Derived1::Grid;
  using Scalar = typename Derived1::Scalar;
  using Value = typename Derived1::Value;
  using Real = typename Derived1::Real;
  using Complex = typename Derived1::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u1.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    if constexpr (std::same_as<Value, RealValued>) {
      return -2 * _u1(-1, iTheta, iPhi) * _u2(1, iTheta, iPhi) +
             _u1(0, iTheta, iPhi) * _u2(0, iTheta, iPhi) -
             2 * _u1(1, iTheta, iPhi) * _u2(-1, iTheta, iPhi);
    } else {
      return -_u1(-1, iTheta, iPhi) * _u2(1, iTheta, iPhi) +
             _u1(0, iTheta, iPhi) * _u2(0, iTheta, iPhi) -
             _u1(1, iTheta, iPhi) * _u2(-1, iTheta, iPhi);
    }
  }

  // Constructors.
  VectorFieldDualityProduct() = delete;
  VectorFieldDualityProduct(const VectorFieldBase<Derived1>& u1,
                            const VectorFieldBase<Derived2>& u2)
      : _u1{u1}, _u2{u2} {
    assert(_u1.FieldSize() == _u2.FieldSize());
  }

  VectorFieldDualityProduct(const VectorFieldDualityProduct&) = default;
  VectorFieldDualityProduct(VectorFieldDualityProduct&&) = default;

  // Assignment.
  VectorFieldDualityProduct& operator=(VectorFieldDualityProduct&) = default;
  VectorFieldDualityProduct& operator=(VectorFieldDualityProduct&&) = default;

 private:
  const VectorFieldBase<Derived1>& _u1;
  const VectorFieldBase<Derived2>& _u2;
};

//-----------------------------------------------------//
//   VectorField product with ScalarField  expression  //
//-----------------------------------------------------//
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
class VectorFieldProductScalarField
    : public VectorFieldBase<
          VectorFieldProductScalarField<Derived1, Derived2>> {
 public:
  using Int = typename VectorFieldBase<
      VectorFieldProductScalarField<Derived1, Derived2>>::Int;
  using Grid = typename Derived1::Grid;
  using Scalar = typename Derived1::Scalar;
  using Value = typename Derived1::Value;
  using Real = typename Derived1::Real;
  using Complex = typename Derived1::Complex;

  // Methods needed to inherit from VectorFieldBase.
  auto GetGrid() const { return _u1.GetGrid(); }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    this->CheckCanonicalIndices(alpha);
    this->CheckPointIndice(iTheta, iPhi);
    return _u1(alpha, iTheta, iPhi) * _u2(iTheta, iPhi);
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  VectorFieldProductScalarField() = delete;
  VectorFieldProductScalarField(const VectorFieldBase<Derived1>& u1,
                                const ScalarFieldBase<Derived2>& u2)
      : _u1{u1}, _u2{u2} {
    assert(_u1.FieldSize() == _u2.FieldSize());
  }

  VectorFieldProductScalarField(const VectorFieldProductScalarField&) = default;
  VectorFieldProductScalarField(VectorFieldProductScalarField&&) = default;

  // Assignment.
  VectorFieldProductScalarField& operator=(
      const VectorFieldProductScalarField&) = default;

  VectorFieldProductScalarField& operator=(VectorFieldProductScalarField&&) =
      default;

 private:
  const VectorFieldBase<Derived1>& _u1;
  const ScalarFieldBase<Derived2>& _u2;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_VECTOR_FIELD_EXPRESSIONS_GUARD_H