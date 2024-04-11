#ifndef GSH_TRANS_VECTOR_FIELD_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "Concepts.h"
#include "GridBase.h"
#include "Indexing.h"
#include "ScalarField.h"
#include "Utility.h"

namespace GSHTrans {

//-------------------------------------------------//
//               Define the base class             //
//-------------------------------------------------//
template <typename _Derived>
class VectorFieldBase : public FieldBase<VectorFieldBase<_Derived>> {
  using Int = std::ptrdiff_t;

 public:
  // Methods related to the grid.
  auto GetGrid() const { return Derived().GetGrid(); }

  // Methods related to the data.
  auto size() const { return 3 * GetGrid().FieldSize(); }
  auto ComponentSize() const { return GetGrid().FieldSize(); }

  auto CanonicalIndices() const { return std::ranges::views::iota(-1, 2); }

  void CheckCanonicalIndices(Int alpha) const { assert(std::abs(alpha) <= 1); }

  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    return Derived().operator()(alpha, iTheta, iPhi);
  }
  auto operator()(Int alpha) const { return Derived().operator()(alpha); }

  void Print() const {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      std::cout << iTheta << " " << iPhi << " ";
      for (auto alpha : CanonicalIndices()) {
        std::cout << operator()(alpha, iTheta, iPhi) << " ";
      }
      std::cout << std::endl;
    }
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

//-------------------------------------------------//
//              Vector component view              //
//-------------------------------------------------//
template <typename Derived>
class VectorFieldComponentView
    : public ScalarFieldBase<VectorFieldComponentView<Derived>> {
  using Int = std::ptrdiff_t;

 public:
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _u(_alpha, iTheta, iPhi);
  }

  // Constructors.
  VectorFieldComponentView(const VectorFieldBase<Derived>& u, Int alpha)
      : _u{u}, _alpha{alpha} {
    this->CheckCanonicalIndices(alpha);
  }

  VectorFieldComponentView(const VectorFieldComponentView&) = default;
  VectorFieldComponentView(VectorFieldComponentView&&) = default;

  // Assignment.
  VectorFieldComponentView& operator=(const VectorFieldComponentView&) =
      default;
  VectorFieldComponentView& operator=(VectorFieldComponentView&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
  Int _alpha;
};

//-------------------------------------------------//
//       Vector field that stores its data         //
//-------------------------------------------------//
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class VectorField : public VectorFieldBase<VectorField<_Grid, _Value>> {
  using Int = std::ptrdiff_t;

 public:
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;

  // Methods needed to inherit from VectorField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return _data[Index(alpha, iTheta, iPhi)];
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  VectorField() = default;

  VectorField(_Grid grid)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->size())} {}

  VectorField(_Grid grid, std::array<Scalar, 3>&& u) : VectorField(grid) {
    for (auto [i, alpha] :
         std::ranges::views::enumerate(this->CanonicalIndices())) {
      this->operator()(alpha) = u[i];
    }
  }

  // Assignment.
  VectorField& operator=(const VectorField&) = default;
  VectorField& operator=(VectorField&&) = default;

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(const VectorFieldBase<Derived>& other) {
    assert(this->size() == other.size());
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(VectorFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Iterators.
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  // Value assignment.
  auto& operator()(Int alpha, Int iTheta, Int iPhi) {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return _data[Index(alpha, iTheta, iPhi)];
  }
  auto operator()(Int alpha) {
    this->CheckCanonicalIndices(alpha);
    auto start = std::next(begin(), Offset(alpha));
    auto finish = std::next(start, this->ComponentSize());
    auto data = std::ranges::subrange(start, finish);
    return ScalarFieldView(_grid, data);
  }

 private:
  _Grid _grid;
  FFTWpp::vector<Scalar> _data;

  template <typename Derived>
  void CopyValues(const VectorFieldBase<Derived>& other) {
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator()(alpha, iTheta, iPhi) = other(alpha, iTheta, iPhi);
      }
    }
  }

  auto Offset(Int alpha) const { return (alpha + 1) * (this->ComponentSize()); }

  auto Index(Int alpha, Int iTheta, int iPhi) const {
    return Offset(alpha) + iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

//-------------------------------------------------//
//      Vector field with a view to its data       //
//-------------------------------------------------//
template <typename _Grid, std::ranges::view _View>
requires requires() {
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>;
  requires std::ranges::random_access_range<_View>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
}
class VectorFieldView : public VectorFieldBase<VectorFieldView<_Grid, _View>> {
  using Int = std::ptrdiff_t;

 public:
  using Grid = _Grid;
  using Scalar = std::ranges::range_value_t<_View>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return _data[Index(alpha, iTheta, iPhi)];
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  VectorFieldView() = default;

  VectorFieldView(_Grid grid, _View data) : _grid{grid}, _data{data} {
    assert(this->size() == _data.size());
  }

  VectorFieldView(const VectorFieldView&) = default;
  VectorFieldView(VectorFieldView&&) = default;

  // Assignment.
  VectorFieldView& operator=(const VectorFieldView&) = default;
  VectorFieldView& operator=(VectorFieldView&&) = default;

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(const VectorFieldBase<Derived>& other) {
    assert(this->size() == other.size());
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(VectorFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Iterators.
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  // Value assignement.
  auto& operator()(Int alpha, Int iTheta, Int iPhi) {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return _data[Index(alpha, iTheta, iPhi)];
  }
  auto operator()(Int alpha) {
    this->CheckCanonicalIndices(alpha);
    auto start = std::next(begin(), Offset(alpha));
    auto finish = std::next(start, this->ComponentSize());
    auto data = std::ranges::subrange(start, finish);
    return ScalarFieldView(_grid, data);
  }

 private:
  _Grid _grid;
  _View _data;

  template <typename Derived>
  void CopyValues(const VectorFieldBase<Derived>& other) {
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator()(alpha, iTheta, iPhi) = other(alpha, iTheta, iPhi);
      }
    }
  }

  auto Offset(Int alpha) const { return (alpha + 1) * (this->ComponentSize()); }

  auto Index(Int alpha, Int iTheta, int iPhi) const {
    return Offset(alpha) + iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

//-----------------------------------------------------//
//            Complex conjugate expression             //
//-----------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
class VectorFieldConjugate
    : public VectorFieldBase<VectorFieldConjugate<Derived>> {
  using Int = std::ptrdiff_t;

 public:
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
  using Int = std::ptrdiff_t;

 public:
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
  using Int = std::ptrdiff_t;

 public:
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
//                 Unary expression                //
//-------------------------------------------------//
template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class VectorFieldUnary
    : public VectorFieldBase<VectorFieldUnary<Derived, Function>> {
  using Int = std::ptrdiff_t;

 public:
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
  VectorFieldUnary() = delete;
  VectorFieldUnary(const VectorFieldBase<Derived>& u, Function f)
      : _u{u}, _f{f} {}

  VectorFieldUnary(const VectorFieldUnary&) = default;
  VectorFieldUnary(VectorFieldUnary&&) = default;

  // Assignment.
  VectorFieldUnary& operator=(VectorFieldUnary&) = default;
  VectorFieldUnary& operator=(VectorFieldUnary&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
  Function _f;
};

//-------------------------------------------------//
//            Unary expression with Scalar         //
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
class VectorFieldUnaryWithScalar
    : public VectorFieldBase<VectorFieldUnaryWithScalar<Derived, Function>> {
  using Int = std::ptrdiff_t;

 public:
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
  VectorFieldUnaryWithScalar() = delete;
  VectorFieldUnaryWithScalar(const VectorFieldBase<Derived>& u, Function f,
                             Scalar s)
      : _u{u}, _f{f}, _s{s} {}

  VectorFieldUnaryWithScalar(const VectorFieldUnaryWithScalar&) = default;
  VectorFieldUnaryWithScalar(VectorFieldUnaryWithScalar&&) = default;

  // Assignment.
  VectorFieldUnaryWithScalar& operator=(VectorFieldUnaryWithScalar&) = default;
  VectorFieldUnaryWithScalar& operator=(VectorFieldUnaryWithScalar&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
  Function _f;
  Scalar _s;
};

//-------------------------------------------------//
//                 Binary expression               //
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
class VectorFieldBinary
    : public VectorFieldBase<VectorFieldBinary<Derived1, Derived2, Function>> {
  using Int = std::ptrdiff_t;

 public:
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
  VectorFieldBinary() = delete;
  VectorFieldBinary(const VectorFieldBase<Derived1>& u1,
                    const VectorFieldBase<Derived2>& u2, Function f)
      : _u1{u1}, _u2{u2}, _f{f} {
    assert(_u1.size() == _u2.size());
  }

  VectorFieldBinary(const VectorFieldBinary&) = default;
  VectorFieldBinary(VectorFieldBinary&&) = default;

  // Assignment.
  VectorFieldBinary& operator=(VectorFieldBinary&) = default;
  VectorFieldBinary& operator=(VectorFieldBinary&&) = default;

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
  using Int = std::ptrdiff_t;

 public:
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
    assert(_u1.size() == _u2.size());
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
  using Int = std::ptrdiff_t;

 public:
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
    assert(_u1.size() == _u2.size());
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
  using Int = std::ptrdiff_t;

 public:
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
    assert(_u1.ComponentSize() == _u2.size());
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

//-----------------------------------------------------//
//              VectorField -> VectorField             //
//-----------------------------------------------------//
template <typename Derived>
auto operator-(const VectorFieldBase<Derived>& u) {
  return VectorFieldUnary(u, [](auto x) { return -x; });
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
  return VectorFieldUnaryWithScalar(u, std::multiplies<>(), s);
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
  return VectorFieldUnaryWithScalar(u, std::divides<>(), s);
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
  return VectorFieldBinary(u1, u2, std::plus<>());
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
  return VectorFieldBinary(u1, u2, std::minus<>());
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
  assert(u1.size() == u2.size());
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

#endif  // GSH_TRANS_VECTOR_FIELD_GUARD_H