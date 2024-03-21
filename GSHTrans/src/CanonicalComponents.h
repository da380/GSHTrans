#ifndef GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H

#include <FFTWpp/Core>
#include <algorithm>
#include <complex>
#include <iostream>
#include <ranges>

#include "Concepts.h"
#include "GridBase.h"

namespace GSHTrans {

//----------------------------------------------------------------//
//                      Define the base class                     //
//----------------------------------------------------------------//

template <typename Derived>
class CanonicalComponentBase {
  using Int = std::ptrdiff_t;

 public:
  // Data access methods.
  auto View() const { return _Derived().View(); }
  auto View() { return _Derived().View(); }

  auto begin() { return _Derived().View().begin(); }
  auto end() { return _Derived().View().end(); }

  auto operator()(Int iTheta, Int iPhi) const {
    auto i = iTheta * this->NumberOfLongitudes() + iPhi;
    return _Derived()[i];
  }

  auto& operator()(Int iTheta, Int iPhi)
  requires std::ranges::output_range<
      typename Derived::view_type,
      std::ranges::range_value_t<typename Derived::view_type>>
  {
    auto i = iTheta * this->NumberOfLongitudes() + iPhi;
    return _Derived()[i];
  }

  // Grid access methods.
  auto Grid() const { return _Derived().Grid(); }

  auto UpperIndex() const { return _Derived().UpperIndex(); }
  auto MaxDegree() const { return Grid().MaxDegree(); }

  auto NumberOfCoLatitudes() const {
    return _Derived().Grid().NumberOfCoLatitudes();
  }
  auto NumberOfLongitudes() const {
    return _Derived().Grid().NumberOfLongitudes();
  }
  auto CoLatitudes() const { return Grid().CoLatitudes(); }
  auto Longitudes() const { return Grid().Longitudes(); }
  auto Points() const { return Grid().Points(); }

  auto CoLatitudeIndices() const { return Grid().CoLatitudeIndices(); }
  auto LongitudeIndices() const { return Grid().LongitudeIndices(); }
  auto PointIndices() const { return Grid().PointIndices(); }

 private:
  auto& _Derived() const { return static_cast<const Derived&>(*this); }
  auto& _Derived() { return static_cast<Derived&>(*this); }
};

//----------------------------------------------------------------//
//             Canonical component storing its data               //
//----------------------------------------------------------------//

template <typename _Grid, RealOrComplexValued Type>
requires std::derived_from<_Grid, GridBase<_Grid>>
class CanonicalComponent
    : public CanonicalComponentBase<CanonicalComponent<_Grid, Type>> {
  using Int = std::ptrdiff_t;
  using Scalar = std::conditional_t<std::same_as<Type, RealValued>,
                                    typename _Grid::real_type,
                                    typename _Grid::complex_type>;
  using Vector = FFTWpp::vector<Scalar>;

 public:
  using grid_type = _Grid;
  using real_or_complex_type = Type;
  using view_type = std::ranges::views::all_t<Vector>;

  // Methods required for CanonicalComponentBase.
  auto UpperIndex() const { return _n; }
  auto Grid() const { return _grid; }
  auto View() { return std::ranges::views::all(_data); }
  auto View() const { return std::ranges::views::all(_data); }
  auto operator[](Int i) const { return _data[i]; }
  auto& operator[](Int i) { return _data[i]; }

  // Constructors and assignement operators.
  CanonicalComponent() = default;

  CanonicalComponent(_Grid grid, Int n)
      : _grid{grid}, _n{n}, _data{Vector(_grid.ComponentSize())} {
    assert(std::ranges::contains(this->Grid().UpperIndices(), _n));
  }

  CanonicalComponent(_Grid grid, Int n, Scalar s)
      : _grid{grid}, _n{n}, _data{Vector(_grid.ComponentSize(), s)} {
    assert(std::ranges::contains(this->Grid().UpperIndices(), _n));
  }

  template <typename Function>
  requires ScalarFunction2D<Function, typename _Grid::real_type, Scalar>
  CanonicalComponent(_Grid grid, Int n, Function f)
      : CanonicalComponent(grid, n) {
    std::ranges::copy(_grid.InterpolateFunction(f), _data.begin());
  }

  CanonicalComponent(const CanonicalComponent&) = default;
  CanonicalComponent(CanonicalComponent&&) = default;

  template <typename Derived>
  CanonicalComponent(const CanonicalComponentBase<Derived>& other)
      : _grid{other.Grid()},
        _n{other.UpperIndex()},
        _data{Vector(other.View().cbegin(), other.View().cend())} {}

  template <typename Derived>
  CanonicalComponent(CanonicalComponentBase<Derived>&& other)
      : CanonicalComponent(other) {}

  CanonicalComponent& operator=(const CanonicalComponent&) = default;
  CanonicalComponent& operator=(CanonicalComponent&&) = default;

  template <typename Derived>
  requires std::convertible_to<std::ranges::range_value_t<Derived>, Scalar>
  CanonicalComponent& operator=(const CanonicalComponentBase<Derived>& other) {
    assert(other.View().size() == this->size());
    std::ranges::copy(other.View() | std::ranges::views::transform(
                                         [](auto x) -> Scalar { return x; }),
                      _data.begin());
    return *this;
  }

  template <typename Derived>
  CanonicalComponent& operator=(CanonicalComponentBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  CanonicalComponent& operator=(Scalar s) {
    std::ranges::copy(std::ranges::views::repeat(s, this->size()),
                      _data.begin());
    return *this;
  }

 private:
  _Grid _grid;
  Int _n;
  Vector _data;
};

// Deduction guide for construction from another base type.
template <typename Derived>
CanonicalComponent(CanonicalComponentBase<Derived>& other)
    -> CanonicalComponent<typename Derived::grid_type,
                          typename Derived::real_or_complex_type>;

template <typename Derived>
CanonicalComponent(CanonicalComponentBase<Derived>&& other)
    -> CanonicalComponent<typename Derived::grid_type,
                          typename Derived::real_or_complex_type>;

// Type aliases for real and complex valued components.
template <typename _Grid>
using RealCanonicalComponent = CanonicalComponent<_Grid, RealValued>;

template <typename _Grid>
using ComplexCanonicalComponent = CanonicalComponent<_Grid, ComplexValued>;

//----------------------------------------------------------------//
//            Canonical component with view to its data           //
//----------------------------------------------------------------//

template <typename _Grid, std::ranges::view _View>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::ranges::input_range<_View>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::real_type>;
}
class CanonicalComponentView
    : public CanonicalComponentBase<CanonicalComponentView<_Grid, _View>>,
      public std::ranges::view_interface<CanonicalComponentView<_Grid, _View>> {
  using Int = std::ptrdiff_t;
  using Real = typename _Grid::real_type;
  using Scalar = std::ranges::range_value_t<_View>;

 public:
  using view_type = _View;
  using real_or_complex_type =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using grid_type = _Grid;

  // Methods required for CanonicalComponentBase.
  auto UpperIndex() const { return _n; }
  auto Grid() const { return _grid; }
  auto View() { return _view; }
  auto View() const { return _view; }
  // auto operator[](Int i) const { return _view[i]; }
  // auto& operator[](Int i) { return _view[i]; }

  // Constructors and assigment operators.
  CanonicalComponentView() = default;

  CanonicalComponentView(_Grid grid, Int n, _View view)
      : _grid{std::move(grid)}, _n{n}, _view{std::move(view)} {
    assert(std::ranges::contains(this->Grid().UpperIndices(), _n));
  }

  CanonicalComponentView(const CanonicalComponentView&) = default;
  CanonicalComponentView(CanonicalComponentView&&) = default;

  CanonicalComponentView& operator=(const CanonicalComponentView&) = default;
  CanonicalComponentView& operator=(CanonicalComponentView&&) = default;

  template <typename Derived>
  requires requires() {
    requires std::ranges::output_range<_View, Scalar>;
    requires std::convertible_to<std::ranges::range_value_t<Derived>, Scalar>;
  }
  CanonicalComponentView& operator=(
      const CanonicalComponentBase<Derived>& other) {
    assert(other.size() == this->size());
    assert(
        std::ranges::contains(this->Grid().UpperIndices(), other.UpperIndex()));
    std::ranges::copy(other.View() | std::ranges::views::transform(
                                         [](auto x) -> Scalar { return x; }),
                      _view.begin());

    return *this;
  }

  template <typename Derived>
  CanonicalComponentView& operator=(CanonicalComponentBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  CanonicalComponentView& operator=(Scalar s)
  requires std::ranges::output_range<_View, Scalar>
  {
    std::ranges::copy(std::ranges::views::repeat(s, this->size()),
                      _view.begin());
    return *this;
  }

 private:
  _Grid _grid;
  Int _n;
  _View _view;
};

// Deduction guide to allow construction from ranges.
template <typename _Grid, std::ranges::viewable_range R>
CanonicalComponentView(_Grid, std::ptrdiff_t, R&&)
    -> CanonicalComponentView<_Grid, std::ranges::views::all_t<R>>;

// Range adaptors to form CanonicalComponentView from a view.
template <typename _Grid>
requires std::derived_from<_Grid, GridBase<_Grid>>
class FormCanonicalComponentView : public std::ranges::range_adaptor_closure<
                                       FormCanonicalComponentView<_Grid>> {
  using Int = std::ptrdiff_t;

 public:
  FormCanonicalComponentView(_Grid grid, Int n)
      : _grid{std::move(grid)}, _n{n} {
    assert(std::ranges::contains(_grid.UpperIndices(), _n));
  }

  template <std::ranges::view View>
  auto operator()(View v) {
    return CanonicalComponentView(_grid, _n, v);
  }

 private:
  _Grid _grid;
  Int _n;
};

// Range adaptors to form a constant CanonicalComponentView from a view.
template <typename _Grid>
requires std::derived_from<_Grid, GridBase<_Grid>>
class FormConstantCanonicalComponentView
    : public std::ranges::range_adaptor_closure<
          FormConstantCanonicalComponentView<_Grid>> {
  using Int = std::ptrdiff_t;

 public:
  FormConstantCanonicalComponentView(_Grid grid, Int n)
      : _grid{std::move(grid)}, _n{n} {
    assert(std::ranges::contains(_grid.UpperIndices(), _n));
  }

  template <std::ranges::view View>
  auto operator()(View v) {
    return CanonicalComponentView(_grid, _n, v | std::ranges::views::as_const);
  }

 private:
  _Grid _grid;
  Int _n;
};

//--------------------------------------------------------//
//           Functions defined on the base class          //
//--------------------------------------------------------//

namespace CanonicalComponentDetails {
template <Field S, Field T>
struct ReturnValueHelper {
  using type =
      std::conditional_t<ComplexFloatingPoint<T>, T,
                         std::conditional_t<ComplexFloatingPoint<S>, S, T>>;
};

template <Field S, Field T>
using ReturnValue = typename ReturnValueHelper<S, T>::type;

}  // namespace CanonicalComponentDetails

// Unary operations.

template <typename Derived, typename Function>
auto CanonicalComponentViewUnary(const CanonicalComponentBase<Derived>& v,
                                 Function f) {
  return v.View() | std::ranges::views::transform(f) |
         FormCanonicalComponentView(v.Grid(), v.UpperIndex());
}

template <typename Derived>
auto operator-(CanonicalComponentBase<Derived>& v) {
  return CanonicalComponentViewUnary(v, [](auto x) { return -x; });
}

template <typename Derived>
auto operator-(CanonicalComponentBase<Derived>&& v) {
  return -v;
}

template <typename Derived>
auto real(CanonicalComponentBase<Derived>& v) {
  using T = std::ranges::range_value_t<CanonicalComponentBase<Derived>>;
  if constexpr (RealFloatingPoint<T>) {
    return CanonicalComponentView(v.Grid(), v.UpperIndex(), v.View());
  } else {
    return CanonicalComponentViewUnary(v, [](auto x) { return std::real(x); });
  }
}

template <typename Derived>
auto real(CanonicalComponentBase<Derived>&& v) {
  return real(v);
}

template <typename Derived>
auto imag(CanonicalComponentBase<Derived>& v) {
  using T = std::ranges::range_value_t<CanonicalComponentBase<Derived>>;
  if constexpr (RealFloatingPoint<T>) {
    return std::ranges::views::repeat(T{0}, v.size()) |
           FormCanonicalComponentView(v.Grid(), v.UpperIndex());
  } else {
    return CanonicalComponentViewUnary(v, [](auto x) { return std::imag(x); });
  }
}

template <typename Derived>
auto imag(CanonicalComponentBase<Derived>&& v) {
  return imag(v);
}

template <typename Derived>
auto conj(const CanonicalComponentBase<Derived>& v) {
  using T = std::ranges::range_value_t<CanonicalComponentBase<Derived>>;
  if constexpr (RealFloatingPoint<T>) {
    return CanonicalComponentView(v.Grid(), v.UpperIndex(), v.View());
  } else {
    return CanonicalComponentViewUnary(v, [](auto x) { return std::conj(x); });
  }
}

template <typename Derived>
auto conj(CanonicalComponentBase<Derived>&& v) {
  return conj(v);
}

template <typename Derived>
auto abs(CanonicalComponentBase<Derived>& v) {
  return CanonicalComponentViewUnary(v, [](auto x) { return std::abs(x); });
}

template <typename Derived>
auto abs(CanonicalComponentBase<Derived>&& v) {
  return abs(v);
}

// Unary operations involving a scalar.

template <typename Derived, Field S>
auto operator*(CanonicalComponentBase<Derived>& v, S s) {
  using View = CanonicalComponentBase<Derived>;
  using T =
      CanonicalComponentDetails::ReturnValue<S,
                                             std::ranges::range_value_t<View>>;
  auto t = T(s);
  std::cout << "----------------\n";
  std::cout << t << std::endl;
  return CanonicalComponentViewUnary(v, [&t](auto x) { return t * x; });
}

template <typename Derived, Field S>
auto operator*(CanonicalComponentBase<Derived>&& v, S s) {
  return v * s;
}

template <typename Derived, Field S>
auto operator*(S s, CanonicalComponentBase<Derived>& v) {
  return v * s;
}

template <typename Derived, Field S>
auto operator*(S s, CanonicalComponentBase<Derived>&& v) {
  return v * s;
}

template <typename Derived, Field S>
auto operator+(CanonicalComponentBase<Derived>& v, S s) {
  using View = CanonicalComponentBase<Derived>;
  using T =
      CanonicalComponentDetails::ReturnValue<S,
                                             std::ranges::range_value_t<View>>;
  auto t = T(s);
  return CanonicalComponentViewUnary(v, [t](auto x) { return t + x; });
}

template <typename Derived, Field S>
auto operator+(CanonicalComponentBase<Derived>&& v, S s) {
  return v + s;
}

template <typename Derived, Field S>
auto operator+(S s, CanonicalComponentBase<Derived>& v) {
  return v + s;
}

template <typename Derived, Field S>
auto operator+(S s, CanonicalComponentBase<Derived>&& v) {
  return v + s;
}

template <typename Derived, Field S>
auto operator-(CanonicalComponentBase<Derived>& v, S s) {
  using View = CanonicalComponentBase<Derived>;
  using T =
      CanonicalComponentDetails::ReturnValue<S,
                                             std::ranges::range_value_t<View>>;
  auto t = T(s);
  return CanonicalComponentViewUnary(v, [t](auto x) { return x - t; });
}

template <typename Derived, Field S>
auto operator-(CanonicalComponentBase<Derived>&& v, S s) {
  return v - s;
}

template <typename Derived, Field S>
auto operator/(CanonicalComponentBase<Derived>& v, S s) {
  using View = CanonicalComponentBase<Derived>;
  using T =
      CanonicalComponentDetails::ReturnValue<S,
                                             std::ranges::range_value_t<View>>;
  auto t = T(s);
  return CanonicalComponentViewUnary(v, [t](auto x) { return x / t; });
}

template <typename Derived, Field S>
auto operator/(CanonicalComponentBase<Derived>&& v, S s) {
  return v / s;
}

template <typename Derived, Field S>
auto pow(CanonicalComponentBase<Derived>& v, S s) {
  using View = CanonicalComponentBase<Derived>;
  using T =
      CanonicalComponentDetails::ReturnValue<S,
                                             std::ranges::range_value_t<View>>;
  auto t = T(s);
  return CanonicalComponentViewUnary(v, [t](auto x) { return std::pow(x, t); });
}

template <typename Derived, Field S>
auto pow(CanonicalComponentBase<Derived>&& v, S s) {
  return pow(v, s);
}

template <typename Derived1, typename Derived2, typename Function>
auto CanonicalComponentViewBinary(Function f,
                                  CanonicalComponentBase<Derived1>& v1,
                                  CanonicalComponentBase<Derived2>& v2) {
  return std::ranges::zip_transform_view(f, v1.View(), v2.View()) |
         FormCanonicalComponentView(v1.Grid(), v1.UpperIndex());
}

template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>& v1,
               CanonicalComponentBase<Derived2>& v2) {
  return CanonicalComponentViewBinary(std::plus<>(), v1, v2);
}

template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>&& v1,
               CanonicalComponentBase<Derived2>& v2) {
  return v1 + v2;
}

template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>& v1,
               CanonicalComponentBase<Derived2>&& v2) {
  return v1 + v2;
}

template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>&& v1,
               CanonicalComponentBase<Derived2>&& v2) {
  return v1 + v2;
}

template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>& v1,
               CanonicalComponentBase<Derived2>& v2) {
  return CanonicalComponentViewBinary(std::minus<>(), v1, v2);
}

template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>&& v1,
               CanonicalComponentBase<Derived2>& v2) {
  return v1 - v2;
}

template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>& v1,
               CanonicalComponentBase<Derived2>&& v2) {
  return v1 - v2;
}

template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>&& v1,
               CanonicalComponentBase<Derived2>&& v2) {
  return v1 - v2;
}

template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>& v1,
               CanonicalComponentBase<Derived2>& v2) {
  return CanonicalComponentViewBinary(std::multiplies<>(), v1, v2);
}

template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>&& v1,
               CanonicalComponentBase<Derived2>& v2) {
  return v1 * v2;
}

template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>& v1,
               CanonicalComponentBase<Derived2>&& v2) {
  return v1 * v2;
}

template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>&& v1,
               CanonicalComponentBase<Derived2>&& v2) {
  return v1 * v2;
}

// Integrate component over S2.
template <typename Derived>
auto Integrate(const CanonicalComponentBase<Derived>& v) {
  using T = std::ranges::range_value_t<CanonicalComponentBase<Derived>>;
  return std::ranges::fold_left(
      std::ranges::views::zip(v.Grid().Weights(), v.View()), T{0},
      [](auto acc, auto term) {
        auto [w, val] = term;
        return acc + val * w;
      });
}

template <typename Derived>
auto Integrate(CanonicalComponentBase<Derived>&& v) {
  return Integrate(v);
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H
