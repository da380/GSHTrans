#ifndef GSH_TRANS_SCALAR_FIELD_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "Concepts.h"
#include "GridBase.h"
#include "Indexing.h"

namespace GSHTrans {

//-------------------------------------------------//
//               Define the base class             //
//-------------------------------------------------//
template <typename _Derived>
class ScalarFieldBase : public FieldBase<ScalarFieldBase<_Derived>> {
  using Int = std::ptrdiff_t;

 public:
  // Methods related to the grid.
  auto GetGrid() const { return Derived().GetGrid(); }

  // Methods related to the data.
  auto size() const { return GetGrid().FieldSize(); }
  auto operator()(Int iTheta, Int iPhi) const {
    return Derived().operator()(iTheta, iPhi);
  }

  void Print() const {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      std::cout << iTheta << " " << iPhi << " " << operator()(iTheta, iPhi)
                << std::endl;
    }
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

//-------------------------------------------------//
//              Constant scalar field              //
//-------------------------------------------------//
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ConstantScalarField
    : public ScalarFieldBase<ConstantScalarField<_Grid, _Value>> {
  using Int = std::ptrdiff_t;

 public:
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _s;
  }

  // Constructors.
  ConstantScalarField() = default;

  ConstantScalarField(_Grid grid, Scalar s) : _grid{grid}, _s{s} {}

  // Assignment.
  ConstantScalarField& operator=(const ConstantScalarField&) = default;
  ConstantScalarField& operator=(ConstantScalarField&&) = default;

  auto& operator=(Scalar s) {
    _s = s;
    return *this;
  }

 private:
  _Grid _grid;
  Scalar _s;
};

//-------------------------------------------------//
//       Scalar field that stores its data         //
//-------------------------------------------------//
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ScalarField : public ScalarFieldBase<ScalarField<_Grid, _Value>> {
  using Int = std::ptrdiff_t;

 public:
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

  // Constructors.
  ScalarField() = default;

  ScalarField(_Grid grid)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->size())} {}

  ScalarField(_Grid grid, Scalar s)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->size(), s)} {}

  template <typename Function>
  requires requires() {
    requires std::invocable<Function, Real, Real>;
    requires std::convertible_to<std::invoke_result_t<Function, Real, Real>,
                                 Scalar>;
  }
  ScalarField(_Grid grid, Function f) : ScalarField(grid) {
    std::ranges::copy(_grid.InterpolateFunction(f), _data.begin());
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  ScalarField(const ScalarFieldBase<Derived>& other)
      : ScalarField(other.GetGrid()) {
    assert(this->size() == other.size());
    CopyValues(other);
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  ScalarField(ScalarFieldBase<Derived>&& other) : ScalarField(other) {}

  ScalarField(const ScalarField&) = default;
  ScalarField(ScalarField&&) = default;

  // Assignment.
  ScalarField& operator=(const ScalarField&) = default;
  ScalarField& operator=(ScalarField&&) = default;

  auto& operator=(Scalar s) {
    std::ranges::copy(std::ranges::views::repeat(s, this->size()),
                      _data.begin());
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator=(const ScalarFieldBase<Derived>& other) {
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator=(ScalarFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Iterators.
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  // Value assignement.
  auto& operator()(Int iTheta, Int iPhi) {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

 private:
  _Grid _grid;
  FFTWpp::vector<Scalar> _data;

  template <typename Derived>
  void CopyValues(const ScalarFieldBase<Derived>& other) {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator()(iTheta, iPhi) = other(iTheta, iPhi);
    }
  }

  auto Index(Int iTheta, int iPhi) const {
    return iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

//-------------------------------------------------//
//      Scalar field with a view to its data       //
//-------------------------------------------------//
template <typename _Grid, std::ranges::view _View>
requires requires() {
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>;
  requires std::ranges::random_access_range<_View>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
}
class ScalarFieldView : public ScalarFieldBase<ScalarFieldView<_Grid, _View>> {
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
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

  // Constructors.
  ScalarFieldView() = default;

  ScalarFieldView(_Grid grid, _View data) : _grid{grid}, _data{data} {
    assert(this->size() == _data.size());
  }

  ScalarFieldView(const ScalarFieldView&) = default;
  ScalarFieldView(ScalarFieldView&&) = default;

  // Assignment.
  ScalarFieldView& operator=(const ScalarFieldView&) = default;
  ScalarFieldView& operator=(ScalarFieldView&&) = default;

  auto& operator=(Scalar s) {
    std::ranges::copy(std::ranges::views::repeat(s, this->size()),
                      _data.begin());
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(const ScalarFieldBase<Derived>& other) {
    assert(this->size() == other.size());
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(ScalarFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Iterators.
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  // Value assignement.
  auto& operator()(Int iTheta, Int iPhi) {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

 private:
  _Grid _grid;
  _View _data;

  template <typename Derived>
  void CopyValues(const ScalarFieldBase<Derived>& other) {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator()(iTheta, iPhi) = other(iTheta, iPhi);
    }
  }

  auto Index(Int iTheta, int iPhi) const {
    return iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

//-------------------------------------------------//
//        Complexification of a real field         //
//-------------------------------------------------//

template <typename Derived>
requires requires() {
  requires std::same_as<typename Derived::Value, RealValued>;
  requires std::same_as<typename Derived::Grid::MRange, All> &&
               std::same_as<typename Derived::Grid::NRange, All>;
}
class ComplexifiedScalarField
    : public ScalarFieldBase<ComplexifiedScalarField<Derived>> {
  using Int = std::ptrdiff_t;

 public:
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar = Complex;
  using Value = ComplexValued;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const -> Complex {
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

//-------------------------------------------------//
//                 Unary expression                //
//-------------------------------------------------//
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
class ScalarFieldUnary
    : public ScalarFieldBase<ScalarFieldUnary<Derived, Function>> {
  using Int = std::ptrdiff_t;

 public:
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
  ScalarFieldUnary() = delete;
  ScalarFieldUnary(const ScalarFieldBase<Derived>& u, Function f)
      : _u{u}, _f{f} {}

  ScalarFieldUnary(const ScalarFieldUnary&) = default;
  ScalarFieldUnary(ScalarFieldUnary&&) = default;

  // Assignment.
  ScalarFieldUnary& operator=(ScalarFieldUnary&) = default;
  ScalarFieldUnary& operator=(ScalarFieldUnary&&) = default;

 private:
  const ScalarFieldBase<Derived>& _u;
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
class ScalarFieldUnaryWithScalar
    : public ScalarFieldBase<ScalarFieldUnaryWithScalar<Derived, Function>> {
  using Int = std::ptrdiff_t;

 public:
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
  ScalarFieldUnaryWithScalar() = delete;
  ScalarFieldUnaryWithScalar(const ScalarFieldBase<Derived>& u, Function f,
                             Scalar s)
      : _u{u}, _f{f}, _s{s} {}

  ScalarFieldUnaryWithScalar(const ScalarFieldUnaryWithScalar&) = default;
  ScalarFieldUnaryWithScalar(ScalarFieldUnaryWithScalar&&) = default;

  // Assignment.
  ScalarFieldUnaryWithScalar& operator=(ScalarFieldUnaryWithScalar&) = default;
  ScalarFieldUnaryWithScalar& operator=(ScalarFieldUnaryWithScalar&&) = default;

 private:
  const ScalarFieldBase<Derived>& _u;
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
                           typename Derived2::Scalar>,
      typename Derived1::Scalar>;
}

class ScalarFieldBinary
    : public ScalarFieldBase<ScalarFieldBinary<Derived1, Derived2, Function>> {
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
    return _f(_u1(iTheta, iPhi), _u2(iTheta, iPhi));
  }

  // Constructors.
  ScalarFieldBinary() = delete;
  ScalarFieldBinary(const ScalarFieldBase<Derived1>& u1,
                    const ScalarFieldBase<Derived2>& u2, Function f)
      : _u1{u1}, _u2{u2}, _f{f} {
    assert(_u1.size() == _u2.size());
  }

  ScalarFieldBinary(const ScalarFieldBinary&) = default;
  ScalarFieldBinary(ScalarFieldBinary&&) = default;

  // Assignment.
  ScalarFieldBinary& operator=(ScalarFieldBinary&) = default;
  ScalarFieldBinary& operator=(ScalarFieldBinary&&) = default;

 private:
  const ScalarFieldBase<Derived1>& _u1;
  const ScalarFieldBase<Derived2>& _u2;
  Function _f;
};

//-----------------------------------------------------//
//               ScalarField -> ScalarField            //
//-----------------------------------------------------//

template <typename Derived>
auto operator-(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return -x; });
}

template <typename Derived>
auto operator-(ScalarFieldBase<Derived>&& u) {
  return -u;
}

template <typename Derived>
auto sqrt(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return std::sqrt(x); });
}

template <typename Derived>
auto sqrt(ScalarFieldBase<Derived>&& u) {
  return sqrt(u);
}

//-----------------------------------------------------//
//            ScalarField -> RealScalarField           //
//-----------------------------------------------------//

template <typename Derived>
auto abs(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return std::abs(x); });
}

template <typename Derived>
auto abs(ScalarFieldBase<Derived>&& u) {
  return abs(u);
}

//-----------------------------------------------------//
//        ComplexScalarField -> RealScalarField        //
//-----------------------------------------------------//

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return std::real(x); });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto real(ScalarFieldBase<Derived>&& u) {
  return real(u);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto imag(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return std::imag(x); });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto imag(ScalarFieldBase<Derived>&& u) {
  return imag(u);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return std::conj(x); });
}

template <typename Derived>
requires std::same_as<typename Derived::Value, ComplexValued>
auto conj(ScalarFieldBase<Derived>&& u) {
  return conj(u);
}

//-----------------------------------------------------//
//        RealScalarField -> ComplexScalarField        //
//-----------------------------------------------------//
template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(const ScalarFieldBase<Derived>& u) {
  return ComplexifiedScalarField(u);
}

template <typename Derived>
requires std::same_as<typename Derived::Value, RealValued>
auto complex(ScalarFieldBase<Derived>&& u) {
  return complex(u);
}

//-----------------------------------------------------//
//         ScalarField x Scalar -> ScalarField         //
//-----------------------------------------------------//
template <typename Derived>
auto operator*(const ScalarFieldBase<Derived>& u, typename Derived::Scalar s) {
  return ScalarFieldUnaryWithScalar(u, std::multiplies<>(), s);
}

template <typename Derived>
auto operator*(typename Derived::Scalar s, const ScalarFieldBase<Derived>& u) {
  return u * s;
}

template <typename Derived>
auto operator*(ScalarFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return u * s;
}

template <typename Derived>
auto operator*(typename Derived::Scalar s, ScalarFieldBase<Derived>&& u) {
  return u * s;
}

template <typename Derived>
auto operator+(const ScalarFieldBase<Derived>& u, typename Derived::Scalar s) {
  return ScalarFieldUnaryWithScalar(u, std::plus<>(), s);
}

template <typename Derived>
auto operator+(typename Derived::Scalar s, const ScalarFieldBase<Derived>& u) {
  return u + s;
}

template <typename Derived>
auto operator+(ScalarFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return u + s;
}

template <typename Derived>
auto operator+(typename Derived::Scalar s, ScalarFieldBase<Derived>&& u) {
  return u + s;
}

template <typename Derived>
auto operator-(const ScalarFieldBase<Derived>& u, typename Derived::Scalar s) {
  return ScalarFieldUnaryWithScalar(u, std::minus<>(), s);
}

template <typename Derived>
auto operator-(ScalarFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return u - s;
}

template <typename Derived>
auto operator/(const ScalarFieldBase<Derived>& u, typename Derived::Scalar s) {
  return ScalarFieldUnaryWithScalar(u, std::divides<>(), s);
}

template <typename Derived>
auto operator/(ScalarFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return u / s;
}

template <typename Derived>
auto pow(const ScalarFieldBase<Derived>& u, typename Derived::Scalar s) {
  return ScalarFieldUnaryWithScalar(
      u, [](auto x, auto y) { return std::pow(x, y); }, s);
}

template <typename Derived>
auto pow(ScalarFieldBase<Derived>&& u, typename Derived::Scalar s) {
  return pow(u, s);
}

//-----------------------------------------------------//
//      ScalarField x ScalarField -> ScalarField       //
//-----------------------------------------------------//

template <typename Derived1, typename Derived2>
auto operator+(const ScalarFieldBase<Derived1>& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return ScalarFieldBinary(u1, u2, std::plus<>());
}

template <typename Derived1, typename Derived2>
auto operator+(const ScalarFieldBase<Derived1>& u1,
               ScalarFieldBase<Derived2>&& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(ScalarFieldBase<Derived1>&& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator+(ScalarFieldBase<Derived1>&& u1, ScalarFieldBase<Derived2>&& u2) {
  return u1 + u2;
}

template <typename Derived1, typename Derived2>
auto operator-(const ScalarFieldBase<Derived1>& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return ScalarFieldBinary(u1, u2, std::minus<>());
}

template <typename Derived1, typename Derived2>
auto operator-(const ScalarFieldBase<Derived1>& u1,
               ScalarFieldBase<Derived2>&& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(ScalarFieldBase<Derived1>&& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator-(ScalarFieldBase<Derived1>&& u1, ScalarFieldBase<Derived2>&& u2) {
  return u1 - u2;
}

template <typename Derived1, typename Derived2>
auto operator*(const ScalarFieldBase<Derived1>& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return ScalarFieldBinary(u1, u2, std::multiplies<>());
}

template <typename Derived1, typename Derived2>
auto operator*(const ScalarFieldBase<Derived1>& u1,
               ScalarFieldBase<Derived2>&& u2) {
  return u1 * u2;
}

template <typename Derived1, typename Derived2>
auto operator*(ScalarFieldBase<Derived1>&& u1,
               const ScalarFieldBase<Derived2>& u2) {
  return u1 * u2;
}

template <typename Derived1, typename Derived2>
auto operator*(ScalarFieldBase<Derived1>&& u1, ScalarFieldBase<Derived2>&& u2) {
  return u1 * u2;
}

//------------------------------------------------------//
//                 ScalarField -> Scalar                //
//------------------------------------------------------//
template <typename Derived>
auto Integrate(const ScalarFieldBase<Derived>& u) {
  using Scalar = typename Derived::Scalar;
  auto w = u.Weights();
  auto i = 0;
  auto sum = Scalar{0};
  for (auto [iTheta, iPhi] : u.PointIndices()) {
    sum += u(iTheta, iPhi) * w[i++];
  }
  return sum;
}

template <typename Derived>
auto Integrate(ScalarFieldBase<Derived>&& u) {
  return Integrate(u);
}

//------------------------------------------------------//
//           ScalarField x ScalarField -> Scalar        //
//------------------------------------------------------//
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(const ScalarFieldBase<Derived1>& u1,
                    const ScalarFieldBase<Derived2>& u2) {
  if constexpr (std::same_as<typename Derived1::Value, RealValued>) {
    return Integral(u1 * u2);
  } else {
    return Integral(conj(u1) * u2);
  }
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(ScalarFieldBase<Derived1>&& u1,
                    const ScalarFieldBase<Derived2>& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(const ScalarFieldBase<Derived1>& u1,
                    ScalarFieldBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto L2InnerProduct(ScalarFieldBase<Derived1>&& u1,
                    ScalarFieldBase<Derived2>&& u2) {
  return L2InnerProduct(u1, u2);
}

//------------------------------------------------------//
//                  ScalarField -> Real                 //
//------------------------------------------------------//
template <typename Derived>
auto L2Norm(const ScalarFieldBase<Derived>& u) {
  return std::sqrt(std::abs(L2InnerProduct(u, u)));
}

template <typename Derived>
auto L2Norm(ScalarFieldBase<Derived>&& u) {
  return L2Norm(u);
}

//-------------------------------------------------------//
//           ScalarField x ScalarField -> bool           //
//-------------------------------------------------------//
template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(const ScalarFieldBase<Derived1>& u1,
                const ScalarFieldBase<Derived2>& u2) {
  assert(u1.size() == u2.size());
  return L2Norm(u1 - u2) <
         std::numeric_limits<typename Derived1::Real>::epsilon();
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(ScalarFieldBase<Derived1>&& u1,
                const ScalarFieldBase<Derived2>& u2) {
  return u1 == u2;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(const ScalarFieldBase<Derived1>& u1,
                ScalarFieldBase<Derived2>&& u2) {
  return u1 == u2;
}

template <typename Derived1, typename Derived2>
requires std::same_as<typename Derived1::Value, typename Derived2::Value>
auto operator==(ScalarFieldBase<Derived1>&& u1,
                ScalarFieldBase<Derived2>&& u2) {
  return u1 == u2;
}

}  // namespace GSHTrans

#endif  // namespace GSHTrans
