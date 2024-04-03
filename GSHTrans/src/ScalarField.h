#ifndef GSH_TRANS_SCALAR_FIELD_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "Concepts.h"
#include "GridBase.h"

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
  auto size() const { return this->GetGrid().FieldSize(); }
  auto operator[](Int i) const { return Derived()[i]; }
  auto operator()(Int iTheta, Int iPhi) const {
    return Derived().operator[](this->Index(iTheta, iPhi));
  }

  void Print() const {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      std::cout << operator()(iTheta, iPhi) << std::endl;
    }
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
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
  auto operator[](Int i) const { return _data[i]; }

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

  // Additional methods.
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }
  auto& operator()(Int iTheta, Int iPhi) {
    return _data[this->Index(iTheta, iPhi)];
  }

 private:
  _Grid _grid;
  FFTWpp::vector<Scalar> _data;

  template <typename Derived>
  void CopyValues(const ScalarFieldBase<Derived>& other) {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      this->operator()(iTheta, iPhi) = other(iTheta, iPhi);
    }
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
  auto operator[](Int i) const { return _data[i]; }

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
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator=(const ScalarFieldBase<Derived>& other) {
    assert(this->size() == other.size());
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator=(ScalarFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Additional methods.
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }
  auto& operator()(Int iTheta, Int iPhi) {
    return _data[this->Index(iTheta, iPhi)];
  }

 private:
  _Grid _grid;
  _View _data;

  template <typename Derived>
  void CopyValues(const ScalarFieldBase<Derived>& other) {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      this->operator()(iTheta, iPhi) = other(iTheta, iPhi);
    }
  }
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
  auto operator[](Int i) const { return _f(_u[i]); }

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
               typename Derived::Real> ||
               std::convertible_to<
                   std::invoke_result_t<Function, typename Derived::Scalar,
                                        typename Derived::Scalar>,
                   typename Derived::Complex>;
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
  auto operator[](Int i) const { return _f(_u[i], _s); }

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
//                 Binary expression                //
//-------------------------------------------------//
template <typename Derived1, typename Derived2, typename Function>
requires requires() {
  requires std::same_as<typename Derived1::Real, typename Derived2::Real>;
  requires std::invocable<Function, typename Derived1::Scalar,
                          typename Derived2::Scalar>;
  requires std::convertible_to<
               std::invoke_result_t<Function, typename Derived1::Scalar,
                                    typename Derived2::Scalar>,
               typename Derived1::Real> ||
               std::convertible_to<
                   std::invoke_result_t<Function, typename Derived1::Scalar,
                                        typename Derived2::Scalar>,
                   typename Derived1::Complex>;
}

class ScalarFieldBinary
    : public ScalarFieldBase<ScalarFieldBinary<Derived1, Derived2, Function>> {
  using Int = std::ptrdiff_t;
  using Grid1 = typename Derived1::Grid;
  using Grid2 = typename Derived2::Grid;
  using Scalar1 = typename Derived1::Scalar;
  using Scalar2 = typename Derived2::Scalar;

 public:
  using Grid = std::conditional_t<
      ComplexFloatingPoint<Scalar1>, Grid1,
      std::conditional_t<ComplexFloatingPoint<Scalar2>, Grid2, Grid1>>;
  using Scalar = std::invoke_result_t<Function, typename Derived1::Scalar,
                                      typename Derived2::Scalar>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Real = typename Derived1::Real;
  using Complex = typename Derived1::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u1.GetGrid(); }
  auto operator[](Int i) const { return _f(_u1[i], _u2[i]); }

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
//                   Field -> Field                    //
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
auto abs(const ScalarFieldBase<Derived>& u) {
  return ScalarFieldUnary(u, [](auto x) { return std::abs(x); });
}

template <typename Derived>
auto abs(ScalarFieldBase<Derived>&& u) {
  return abs(u);
}

//-----------------------------------------------------//
//               Field x Scalar -> Field               //
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

//-----------------------------------------------------//
//                Field x Field -> Field               //
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

}  // namespace GSHTrans

#endif  // namespace GSHTrans
