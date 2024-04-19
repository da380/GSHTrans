#ifndef GSH_TRANS_MATRIX_FIELD_GUARD_H
#define GSH_TRANS_MATRIX_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "Concepts.h"
#include "GridBase.h"
#include "Indexing.h"
#include "ScalarField.h"
#include "VectorField.h"

namespace GSHTrans {

//-------------------------------------------------//
//               Define the base class             //
//-------------------------------------------------//
template <typename _Derived>
class MatrixFieldBase : public FieldBase<MatrixFieldBase<_Derived>> {
 public:
  using Int = typename FieldBase<MatrixFieldBase<_Derived>>::Int;

  // Methods related to the grid.
  auto GetGrid() const { return Derived().GetGrid(); }

  // Methods related to the data.
  auto ComponentSize() const { return GetGrid().FieldSize(); }

  auto CanonicalIndices() const {
    return std::ranges::views::cartesian_product(
        std::ranges::views::iota(-1, 2), std::ranges::views::iota(-1, 2));
  }

  void CheckCanonicalIndices(Int alpha, Int beta) const {
    assert(std::abs(alpha) <= 1);
    assert(std::abs(beta) <= 1);
  }

  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    return Derived().operator()(alpha, beta, iTheta, iPhi);
  }
  auto operator()(Int alpha, Int beta) const {
    return Derived().operator()(alpha, beta);
  }

  void Print() const {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      for (auto alpha = -1; alpha <= 1; alpha++) {
        for (auto beta = -1; beta <= 1; beta++) {
          std::cout << operator()(alpha, beta, iTheta, iPhi) << " ";
        }
        std::cout << std::endl;
      }
      std::cout << "-------------" << std::endl;
    }
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

//-------------------------------------------------//
//              Matrix component view              //
//-------------------------------------------------//
template <typename Derived>
class MatrixFieldComponentView
    : public ScalarFieldBase<MatrixFieldComponentView<Derived>> {
 public:
  using Int = typename ScalarFieldBase<MatrixFieldComponentView<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    return _A(_alpha, _beta, iTheta, iPhi);
  }

  // Constructors.
  MatrixFieldComponentView(const MatrixFieldBase<Derived>& A, Int alpha,
                           Int beta)
      : _A{A}, _alpha{alpha}, _beta{beta} {}

  MatrixFieldComponentView(const MatrixFieldComponentView&) = default;
  MatrixFieldComponentView(MatrixFieldComponentView&&) = default;

  // Assignment.
  MatrixFieldComponentView& operator=(const MatrixFieldComponentView&) =
      default;
  MatrixFieldComponentView& operator=(MatrixFieldComponentView&&) = default;

 private:
  const VectorFieldBase<Derived>& _A;
  Int _alpha;
  Int _beta;
};

//-------------------------------------------------//
//       Matrix field that stores its data         //
//-------------------------------------------------//
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class MatrixField : public MatrixFieldBase<MatrixField<_Grid, _Value>> {
 public:
  using Int = typename MatrixFieldBase<MatrixField<_Grid, _Value>>::Int;
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;

  // Methods needed to inherit from MatrixField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    this->CheckCanonicalIndices(alpha, beta);
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(alpha, beta, iTheta, iPhi)];
  }
  auto operator()(Int alpha, Int beta) const {
    this->CheckCanonicalIndices(alpha, beta);
    return MatrixFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  MatrixField() = default;

  MatrixField(_Grid grid)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(size())} {}

  MatrixField(_Grid grid, std::array<Scalar, 9>&& A) : MatrixField(grid) {
    for (auto [i, index] :
         std::ranges::views::enumerate(this->CanonicalIndices())) {
      auto [alpha, beta] = index;
      this->operator()(alpha, beta) = A[i];
    }
  }

  // Assignment.
  MatrixField& operator=(const MatrixField&) = default;
  MatrixField& operator=(MatrixField&&) = default;

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(const MatrixFieldBase<Derived>& other) {
    assert(this->ComponentSize() == other.ComponentSize());
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(MatrixFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Methods to make it a range.
  auto size() const { return 9 * this->ComponentSize(); }
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  // Value assignment.
  auto& operator()(Int alpha, Int beta, Int iTheta, Int iPhi) {
    this->CheckCanonicalIndices(alpha, beta);
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(alpha, beta, iTheta, iPhi)];
  }
  auto operator()(Int alpha, Int beta) {
    this->CheckCanonicalIndices(alpha, beta);
    auto start = std::next(begin(), Offset(alpha, beta));
    auto finish = std::next(start, this->ComponentSize());
    auto data = std::ranges::subrange(start, finish);
    return ScalarFieldView(_grid, data);
  }

 private:
  _Grid _grid;
  FFTWpp::vector<Scalar> _data;

  auto Offset(Int alpha, Int beta) const {
    return (3 * (alpha + 1) + (beta + 1)) * (this->ComponentSize());
  }

  template <typename Derived>
  void CopyValues(const MatrixFieldBase<Derived>& other) {
    for (auto [alpha, beta] : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator()(alpha, beta, iTheta, iPhi) =
            other(alpha, beta, iTheta, iPhi);
      }
    }
  }

  auto Index(Int alpha, Int beta, Int iTheta, int iPhi) const {
    return Offset(alpha, beta) + iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

//-----------------------------------------------------//
//            Complex conjugate expression             //
//-----------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
class MatrixFieldConjugate
    : public MatrixFieldBase<MatrixFieldConjugate<Derived>> {
 public:
  using Int = typename MatrixFieldBase<MatrixFieldConjugate<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from MatrixFieldBase.
  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha, beta);
    return static_cast<Real>(MinusOneToPower(alpha + beta)) *
           std::conj(_A(-alpha, -beta, iTheta, iPhi));
  }
  auto operator()(Int alpha, Int beta) const {
    this->CheckCanonicalIndices(alpha, beta);
    return MatrixFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  MatrixFieldConjugate() = delete;
  MatrixFieldConjugate(const MatrixFieldBase<Derived>& A) : _A{A} {}

  MatrixFieldConjugate(const MatrixFieldConjugate&) = default;
  MatrixFieldConjugate(MatrixFieldConjugate&&) = default;

  // Assignment.
  MatrixFieldConjugate& operator=(MatrixFieldConjugate&) = default;
  MatrixFieldConjugate& operator=(MatrixFieldConjugate&&) = default;

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

  // Methods needed to inherit from MatrixFieldBase.
  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha, beta);
    if (std::same_as<Value, RealValued>) {
      return _A(beta, alpha, iTheta, iPhi);
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
//      Pointswise PointwiseUnary expression       //
//-------------------------------------------------//
template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class MatrixFieldPointwiseUnary
    : public MatrixFieldBase<MatrixFieldPointwiseUnary<Derived, Function>> {
 public:
  using Int = typename MatrixFieldBase<
      MatrixFieldPointwiseUnary<Derived, Function>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = std::invoke_result_t<Function, typename Derived::Scalar>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from MatrixField Base.
  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    return _f(_A(alpha, beta, iTheta, iPhi));
  }
  auto operator()(Int alpha, Int beta) const {
    return MatrixFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  MatrixFieldPointwiseUnary() = delete;
  MatrixFieldPointwiseUnary(const MatrixFieldBase<Derived>& A, Function f)
      : _A{A}, _f{f} {}

  MatrixFieldPointwiseUnary(const MatrixFieldPointwiseUnary&) = default;
  MatrixFieldPointwiseUnary(MatrixFieldPointwiseUnary&&) = default;

  // Assignment.
  MatrixFieldPointwiseUnary& operator=(MatrixFieldPointwiseUnary&) = default;
  MatrixFieldPointwiseUnary& operator=(MatrixFieldPointwiseUnary&&) = default;

 private:
  const MatrixFieldBase<Derived>& _A;
  Function _f;
};

//-------------------------------------------------//
//      PointwiseUnary expression with Scalar      //
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
class MatrixFieldPointwiseUnaryWithScalar
    : public MatrixFieldBase<
          MatrixFieldPointwiseUnaryWithScalar<Derived, Function>> {
 public:
  using Int = typename MatrixFieldBase<
      MatrixFieldPointwiseUnaryWithScalar<Derived, Function>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from MatrixFieldBase.
  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha, beta);
    return _f(_A(alpha, beta, iTheta, iPhi), _s);
  }
  auto operator()(Int alpha, Int beta) const {
    this->CheckCanonicalIndices(alpha);
    return MatrixFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  MatrixFieldPointwiseUnaryWithScalar() = delete;
  MatrixFieldPointwiseUnaryWithScalar(const MatrixFieldBase<Derived>& A,
                                      Function f, Scalar s)
      : _A{A}, _f{f}, _s{s} {}

  MatrixFieldPointwiseUnaryWithScalar(
      const MatrixFieldPointwiseUnaryWithScalar&) = default;
  MatrixFieldPointwiseUnaryWithScalar(MatrixFieldPointwiseUnaryWithScalar&&) =
      default;

  // Assignment.
  MatrixFieldPointwiseUnaryWithScalar& operator=(
      MatrixFieldPointwiseUnaryWithScalar&) = default;
  MatrixFieldPointwiseUnaryWithScalar& operator=(
      MatrixFieldPointwiseUnaryWithScalar&&) = default;

 private:
  const MatrixFieldBase<Derived>& _A;
  Function _f;
  Scalar _s;
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
    assert(_A1.ComponentSize() == _A2.ComponentSize());
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
      auto u = -_A(-1, 1) + _A(0, 0) - _A(+1, 1);
      return u(iTheta, iPhi);
    } else {
      auto u = -2 * _A(1, -1) + _A(0, 0);
      return u(iTheta, iPhi);
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
      return 0;
    } else {
      auto u =
          _A(0, 0) * (_A(-1, +1) * _A(-1, +1) + _A(+1, -1) * _A(+1, -1) -
                      _A(-1, -1) * _A(-1, -1) - _A(+1, +1) * _A(+1, +1)) +
          2 * _A(+1, -1) * (_A(-1, 0) * _A(0, -1) + _A(+1, 0) * _A(0, +1)) +
          2 * _A(-1, +1) * (_A(-1, 0) * _A(0, +1) - _A(+1, 0) * _A(0, -1)) +
          2 * _A(-1, -1) * (_A(-1, 0) * _A(0, +1) + _A(+1, 0) * _A(0, -1)) +
          2 * _A(+1, +1) * (_A(+1, 0) * _A(0, +1) - _A(-1, 0) * _A(0, -1));
      return u(iTheta, iPhi);
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
      auto v = -_A(alpha, -1) * _u(+1) + _A(alpha, 0) * _u(0) -
               _A(alpha, 1) * _u(-1);
      return v(iTheta, iPhi);
    } else {
      switch (alpha) {
        case -1: {
          auto v = -(_A(+1, -1) + _A(+1, +1)) * _u(-1) + _A(-1, 0) * _u(0) +
                   (_A(-1, -1) - _A(-1, +1)) * _u(+1);
          return v(iTheta, iPhi);
        }
        case 0: {
          auto v = 2 * _A(0, -1) * _u(-1) + _A(0, 0) * _u(0) +
                   2 * _A(0, +1) * _u(+1);
          return v(iTheta, iPhi);
        }
        case +1: {
          auto v = (_A(-1, -1) + _A(-1, +1)) * _u(-1) + _A(+1, 0) * _u(0) +
                   (_A(+1, +1) - _A(+1, -1)) * _u(+1);
          return v(iTheta, iPhi);
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
    assert(_A.ComponentSize() == _u.ComponentSize());
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

}  // namespace GSHTrans

#endif  // GSH_TRANS_MATRIX_FIELD_GUARD_H