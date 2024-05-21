#ifndef GSH_TRANS_SCALAR_FIELD_EXPRESSIONS_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_EXPRESSIONS_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "ScalarFieldBase.h"

namespace GSHTrans {

// Complexification of a real field.
template <typename Derived>
requires requires() {
  requires std::same_as<typename Derived::Value, RealValued>;
  requires std::same_as<typename Derived::Grid::MRange, All> &&
               std::same_as<typename Derived::Grid::NRange, All>;
}
class ComplexifiedScalarField
    : public ScalarFieldBase<ComplexifiedScalarField<Derived>> {
 public:
  using Int = typename ScalarFieldBase<ComplexifiedScalarField<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar = Complex;
  using Value = ComplexValued;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) -> Complex {
    this->CheckPointIndices(iTheta, iPhi);
    return _u(iTheta, iPhi);
  }

  // Constructors.
  ComplexifiedScalarField(const ScalarFieldBase<Derived>& u) : _u{u} {}

  ComplexifiedScalarField(const ComplexifiedScalarField&) = default;
  ComplexifiedScalarField(ComplexifiedScalarField&&) = default;

  // Assignment.
  ComplexifiedScalarField& operator=(const ComplexifiedScalarField&) = default;
  ComplexifiedScalarField& operator=(ComplexifiedScalarField&&) = default;

 private:
  const ScalarFieldBase<Derived>& _u;
};

// Pointwise unary expressions.
template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
               std::invoke_result_t<Function, typename Derived::Scalar>,
               typename Derived::Scalar> ||
               (std::convertible_to<
                    std::invoke_result_t<Function, typename Derived::Scalar>,
                    typename Derived::Complex> &&
                std::same_as<typename Derived::Grid::MRange, All> &&
                std::same_as<typename Derived::Grid::MRange, All>);
}
class ScalarFieldPointwiseUnary
    : public ScalarFieldBase<ScalarFieldPointwiseUnary<Derived, Function>> {
 public:
  using Int = typename ScalarFieldBase<
      ScalarFieldPointwiseUnary<Derived, Function>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = std::invoke_result_t<Function, typename Derived::Scalar>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u(iTheta, iPhi));
  }

  // Constructors.
  ScalarFieldPointwiseUnary() = delete;
  ScalarFieldPointwiseUnary(const ScalarFieldBase<Derived>& u, Function f)
      : _u{u}, _f{f} {}

  ScalarFieldPointwiseUnary(const ScalarFieldPointwiseUnary&) = default;
  ScalarFieldPointwiseUnary(ScalarFieldPointwiseUnary&&) = default;

  // Assignment.
  ScalarFieldPointwiseUnary& operator=(ScalarFieldPointwiseUnary&) = default;
  ScalarFieldPointwiseUnary& operator=(ScalarFieldPointwiseUnary&&) = default;

 private:
  const ScalarFieldBase<Derived>& _u;
  Function _f;
};

// Pointwise unary expression parameterised by a scalar.
template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar,
                          typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar,
                           typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class ScalarFieldPointwiseUnaryWithScalar
    : public ScalarFieldBase<
          ScalarFieldPointwiseUnaryWithScalar<Derived, Function>> {
 public:
  using Int = typename ScalarFieldBase<
      ScalarFieldPointwiseUnaryWithScalar<Derived, Function>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = std::invoke_result_t<Function, typename Derived::Scalar,
                                      typename Derived::Scalar>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u(iTheta, iPhi), _s);
  }

  // Constructors.
  ScalarFieldPointwiseUnaryWithScalar() = delete;
  ScalarFieldPointwiseUnaryWithScalar(const ScalarFieldBase<Derived>& u,
                                      Function f, Scalar s)
      : _u{u}, _f{f}, _s{s} {}

  ScalarFieldPointwiseUnaryWithScalar(
      const ScalarFieldPointwiseUnaryWithScalar&) = default;
  ScalarFieldPointwiseUnaryWithScalar(ScalarFieldPointwiseUnaryWithScalar&&) =
      default;

  // Assignment.
  ScalarFieldPointwiseUnaryWithScalar& operator=(
      ScalarFieldPointwiseUnaryWithScalar&) = default;
  ScalarFieldPointwiseUnaryWithScalar& operator=(
      ScalarFieldPointwiseUnaryWithScalar&&) = default;

 private:
  const ScalarFieldBase<Derived>& _u;
  Function _f;
  Scalar _s;
};

// Pointwise binary expressions.
template <typename Derived1, typename Derived2, typename Function>
requires requires() {
  requires std::same_as<typename Derived1::Real, typename Derived2::Real>;
  requires std::same_as<typename Derived1::Value, typename Derived2::Value>;
  requires std::invocable<Function, typename Derived1::Scalar,
                          typename Derived2::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived1::Scalar,
                           typename Derived2::Scalar>,
      typename Derived1::Scalar>;
}
class ScalarFieldPointwiseBinary
    : public ScalarFieldBase<
          ScalarFieldPointwiseBinary<Derived1, Derived2, Function>> {
 public:
  using Int = typename ScalarFieldBase<
      ScalarFieldPointwiseBinary<Derived1, Derived2, Function>>::Int;
  using Grid = typename Derived1::Grid;
  using Scalar = typename Derived1::Scalar;
  using Value = typename Derived1::Value;
  using Real = typename Derived1::Real;
  using Complex = typename Derived1::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u1.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u1(iTheta, iPhi), _u2(iTheta, iPhi));
  }

  // Constructors.
  ScalarFieldPointwiseBinary() = delete;
  ScalarFieldPointwiseBinary(const ScalarFieldBase<Derived1>& u1,
                             const ScalarFieldBase<Derived2>& u2, Function f)
      : _u1{u1}, _u2{u2}, _f{f} {
    assert(_u1.FieldSize() == _u2.FieldSize());
  }

  ScalarFieldPointwiseBinary(const ScalarFieldPointwiseBinary&) = default;
  ScalarFieldPointwiseBinary(ScalarFieldPointwiseBinary&&) = default;

  // Assignment.
  ScalarFieldPointwiseBinary& operator=(ScalarFieldPointwiseBinary&) = default;
  ScalarFieldPointwiseBinary& operator=(ScalarFieldPointwiseBinary&&) = default;

 private:
  const ScalarFieldBase<Derived1>& _u1;
  const ScalarFieldBase<Derived2>& _u2;
  Function _f;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_EXPRESSIONS_GUARD_H