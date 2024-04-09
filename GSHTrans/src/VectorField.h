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
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    return Derived().operator()(alpha, iTheta, iPhi);
  }
  auto operator()(Int alpha) const { return Derived().operator()(alpha); }

  void Print() const {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      std::cout << iTheta << " " << iPhi << " " << operator()(-1, iTheta, iPhi)
                << " " << operator()(0, iTheta, iPhi) << " "
                << operator()(1, iTheta, iPhi) << std::endl;
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
    return _u(_alpha, iTheta, iPhi);
  }

  // Constructors.
  VectorFieldComponentView(const VectorFieldBase<Derived>& u, Int alpha)
      : _u{u}, _alpha{alpha} {}

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
    return _data[Index(alpha, iTheta, iPhi)];
  }
  auto operator()(Int alpha) const {
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  VectorField() = default;

  VectorField(_Grid grid)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->size())} {}

  VectorField(_Grid grid, Scalar minus, Scalar zero, Scalar plus)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->size())} {
    this->operator()(-1) = minus;
    this->operator()(0) = zero;
    this->operator()(1) = plus;
  }

  // Assignment.
  VectorField& operator=(const VectorField&) = default;
  VectorField& operator=(VectorField&&) = default;

  // Iterators.
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  // Value assignment.
  auto& operator()(Int alpha, Int iTheta, Int iPhi) {
    return _data[Index(alpha, iTheta, iPhi)];
  }
  auto operator()(Int alpha) {
    const auto size = this->ComponentSize();
    auto start = std::next(begin(), size * (alpha + 1));
    auto finish = std::next(start, size);
    auto data = std::ranges::subrange(start, finish);
    return ScalarFieldView(_grid, data);
  }

 private:
  _Grid _grid;
  FFTWpp::vector<Scalar> _data;

  auto Index(Int alpha, Int iTheta, int iPhi) const {
    return (this->ComponentSize()) * (alpha + 1) +
           iTheta * this->NumberOfLongitudes() + iPhi;
  }
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
    constexpr auto ii = Complex(0, 1);
    switch (alpha) {
      case -1:
        return _u(1, iTheta, iPhi) - ii * _u(-1, iTheta, iPhi);
      case 0:
        return _u(0, iTheta, iPhi);
      case 1:
        return _u(1, iTheta, iPhi) + ii * _u(-1, iTheta, iPhi);
      default:
        return 0;
    }
  }
  auto operator()(Int alpha) const {
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

template <typename Derived>
auto Complexify(const VectorFieldBase<Derived>& u) {
  return ComplexifiedVectorField(u);
}

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
    constexpr auto half = static_cast<Real>(1) / static_cast<Real>(2);
    switch (alpha) {
      case -1:
        return half * std::imag(_u(1, iTheta, iPhi) - _u(-1, iTheta, iPhi));
      case 0:
        return std::real(_u(0, iTheta, iPhi));
      case 1:
        return half * std::real(_u(1, iTheta, iPhi) + _u(-1, iTheta, iPhi));
      default:
        return 0;
    }
  }
  auto operator()(Int alpha) const {
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
               typename Derived::Real> ||
               std::convertible_to<
                   std::invoke_result_t<Function, typename Derived::Scalar>,
                   typename Derived::Complex>;
}
class VectorFieldUnary
    : public VectorFieldBase<VectorFieldUnary<Derived, Function>> {
  using Int = std::ptrdiff_t;

 public:
  using Grid = typename Derived::Grid;
  using Scalar = std::invoke_result_t<Function, typename Derived::Scalar>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from VectorField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    return _f(_u(alpha, iTheta, iPhi));
  }
  auto operator()(Int alpha) const {
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
//          VectorField -> RealVectorField             //
//-----------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto abs(const VectorFieldBase<Derived>& u) {
  return real(u);
}

//-----------------------------------------------------//
//         ComplexVectorField -> RealVectorField       //
//-----------------------------------------------------//

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(const VectorFieldBase<Derived>& u) {
  return RealifiedVectorField(u);
}

//-----------------------------------------------------//
//         RealVectorField -> ComplexVectorField       //
//-----------------------------------------------------//

template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(const VectorFieldBase<Derived>& u) {
  return ComplexifiedVectorField(u);
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_VECTOR_FIELD_GUARD_H