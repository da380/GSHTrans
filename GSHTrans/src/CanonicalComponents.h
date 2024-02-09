#ifndef GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H

#include <FFTWpp/All>
#include <algorithm>
#include <complex>
#include <functional>
#include <iostream>
#include <memory>
#include <ranges>
#include <type_traits>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

//----------------------------------------------------------------//
//                      Define the base class                     //
//----------------------------------------------------------------//

template <typename Derived>
class CanonicalComponentBase
    : public std::ranges::view_interface<CanonicalComponentBase<Derived>> {
  using Int = std::ptrdiff_t;

 public:
  // Return upper index.
  auto UpperIndex() const { return _Derived()._UpperIndex(); }

  // Return grid information,
  auto GridPointer() const { return _Derived()._grid; }

  auto View() const { return _Derived()._View(); }
  auto View() { return _Derived()._View(); }

  auto NumberOfCoLatitudes() const {
    return GridPointer()->NumberOfCoLatitudes();
  }
  auto NumberOfLongitudes() const {
    return GridPointer()->NumberOfLongitudes();
  }
  auto CoLatitudes() const { return GridPointer()->CoLatitudes(); }
  auto Longitudes() const { return GridPointer()->Longitudes(); }
  auto Points() const { return GridPointer()->Points(); }

  // Data access functions.
  auto begin() { return _Derived()._begin(); }
  auto end() { return _Derived()._end(); }

  auto operator()(Int iTheta, Int iPhi) const {
    auto i = iTheta * NumberOfLongitudes() + iPhi;
    return this->operator[](i);
  }

  auto& operator()(Int iTheta, Int iPhi)
  requires std::ranges::output_range<
      typename Derived::view_type,
      std::ranges::range_value_t<typename Derived::view_type>>
  {
    auto i = iTheta * NumberOfLongitudes() + iPhi;
    return this->operator[](i);
  }

 private:
  auto& _Derived() const { return static_cast<const Derived&>(*this); }
  auto& _Derived() { return static_cast<Derived&>(*this); }
};

//----------------------------------------------------------------//
//             Canonical component storing its data               //
//----------------------------------------------------------------//

template <typename Grid, RealOrComplexValued Type>
class CanonicalComponent
    : public CanonicalComponentBase<CanonicalComponent<Grid, Type>> {
  using Int = std::ptrdiff_t;
  using Scalar =
      std::conditional_t<std::same_as<Type, RealValued>,
                         typename Grid::real_type, typename Grid::complex_type>;
  using Vector = FFTWpp::vector<Scalar>;

 public:
  using grid_type = Grid;
  using view_type = std::ranges::views::all_t<Vector>;

  CanonicalComponent() = default;

  CanonicalComponent(std::shared_ptr<Grid> grid, Int n)
      : _grid{std::move(grid)}, _n{n}, _data{Vector(_grid->ComponentSize())} {
    assert(std::ranges::contains(this->GridPointer()->UpperIndices(), _n));
  }

  CanonicalComponent(std::shared_ptr<Grid> grid, Int n, Scalar s)
      : _grid{std::move(grid)},
        _n{n},
        _data{Vector(_grid->ComponentSize(), s)} {
    assert(std::ranges::contains(this->GridPointer()->UpperIndices(), _n));
  }

  template <typename Function>
  requires ScalarFunction2D<Function, typename Grid::real_type, Scalar>
  CanonicalComponent(std::shared_ptr<Grid> grid, Int n, Function f)
      : CanonicalComponent(grid, n) {
    auto fView = this->Points() |
                 std::ranges::views::transform([f](auto point) -> Scalar {
                   auto [theta, phi] = point;
                   return f(theta, phi);
                 });
    std::ranges::copy(fView, _data.begin());
  }

  CanonicalComponent(const CanonicalComponent&) = default;
  CanonicalComponent(CanonicalComponent&&) = default;

  template <typename Derived>
  CanonicalComponent(const CanonicalComponentBase<Derived>& other)
      : _grid{other.GridPointer()},
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
  std::shared_ptr<Grid> _grid;
  Int _n;
  Vector _data;

  auto _UpperIndex() const { return _n; }

  auto _begin() { return _data.begin(); }
  auto _end() { return _data.end(); }

  auto _View() { return std::ranges::views::all(_data); }
  auto _View() const {
    return std::ranges::views::all(_data) | std::ranges::views::as_const;
  }

  friend class CanonicalComponentBase<CanonicalComponent<Grid, Type>>;
};

// Deduction guide for construction from another base type.
template <typename Derived>
CanonicalComponent(CanonicalComponentBase<Derived>& other)
    -> CanonicalComponent<
        typename Derived::grid_type,
        std::conditional_t<RealFloatingPoint<std::ranges::range_value_t<
                               typename Derived::view_type>>,
                           RealValued, ComplexValued>>;

// Type aliases for real and complex valued components.
template <typename Grid>
using RealCanonicalComponent = CanonicalComponent<Grid, RealValued>;

template <typename Grid>
using ComplexCanonicalComponent = CanonicalComponent<Grid, ComplexValued>;

//----------------------------------------------------------------//
//            Canonical component with view to its data           //
//----------------------------------------------------------------//

template <typename Grid, std::ranges::view View>
requires requires() {
  requires std::ranges::input_range<View>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<View>>,
                        typename Grid::real_type>;
}
class CanonicalComponentView
    : public CanonicalComponentBase<CanonicalComponentView<Grid, View>> {
  using Int = std::ptrdiff_t;
  using Real = typename Grid::real_type;
  using Scalar = std::ranges::range_value_t<View>;

 public:
  using view_type = View;
  using grid_type = Grid;

  CanonicalComponentView() = default;

  CanonicalComponentView(std::shared_ptr<Grid> grid, Int n, View view)
      : _grid{std::move(grid)}, _n{n}, _view{std::move(view)} {
    assert(std::ranges::contains(this->GridPointer()->UpperIndices(), _n));
  }

  CanonicalComponentView(const CanonicalComponentView&) = default;
  CanonicalComponentView(CanonicalComponentView&&) = default;

  CanonicalComponentView& operator=(const CanonicalComponentView&) = default;
  CanonicalComponentView& operator=(CanonicalComponentView&&) = default;

  template <typename Derived>
  requires requires() {
    requires std::ranges::output_range<View, Scalar>;
    requires std::convertible_to<std::ranges::range_value_t<Derived>, Scalar>;
  }
  CanonicalComponentView& operator=(
      const CanonicalComponentBase<Derived>& other) {
    assert(other.View().size() == this->size());
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
  requires std::ranges::output_range<View, Scalar>
  {
    std::ranges::copy(std::ranges::views::repeat(s, this->size()),
                      _view.begin());

    return *this;
  }

 private:
  std::shared_ptr<Grid> _grid;
  Int _n;
  View _view;

  auto _UpperIndex() const { return _n; }

  auto _begin() { return _view.begin(); }
  auto _end() { return _view.end(); }

  auto _View() { return *this; }
  auto _View() const { return *this; }

  friend class CanonicalComponentBase<CanonicalComponentView<Grid, View>>;
};

// Deduction guide to allow construction from ranges.
template <typename Grid, std::ranges::viewable_range R>
CanonicalComponentView(std::shared_ptr<Grid>, R&&)
    -> CanonicalComponentView<Grid, std::ranges::views::all_t<R>>;

// Range adaptors to form CanonicalComponentView from a view.
template <typename Grid>
class FormCanonicalComponentView : public std::ranges::range_adaptor_closure<
                                       FormCanonicalComponentView<Grid>> {
  using Int = std::ptrdiff_t;

 public:
  FormCanonicalComponentView(std::shared_ptr<Grid> grid, Int n)
      : _grid{std::move(grid)}, _n{n} {
    assert(std::ranges::contains(_grid->UpperIndices(), _n));
  }

  template <std::ranges::view View>
  auto operator()(View v) {
    return CanonicalComponentView(_grid, _n, v);
  }

 private:
  std::shared_ptr<Grid> _grid;
  Int _n;
};

template <typename Grid>
class FormConstantCanonicalComponentView
    : public std::ranges::range_adaptor_closure<
          FormConstantCanonicalComponentView<Grid>> {
  using Int = std::ptrdiff_t;

 public:
  FormConstantCanonicalComponentView(std::shared_ptr<Grid> grid, Int n)
      : _grid{std::move(grid)}, _n{n} {
    assert(std::ranges::contains(_grid->UpperIndices(), _n));
  }

  template <std::ranges::view View>
  auto operator()(View v) {
    return CanonicalComponentView(_grid, _n, v | std::ranges::views::as_const);
  }

 private:
  std::shared_ptr<Grid> _grid;
  Int _n;
};

template <typename Grid, typename Function>
auto InterpolateCanonicalComponentView(std::shared_ptr<Grid> grid,
                                       std::ptrdiff_t n, Function f) {
  return grid->Points() | std::ranges::views::transform([f](auto point) {
           auto [theta, phi] = point;
           return f(theta, phi);
         }) |
         FormCanonicalComponentView(grid, n);
}

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
         FormCanonicalComponentView(v.GridPointer(), v.UpperIndex());
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
    return CanonicalComponentView(v.GridPointer(), v.UpperIndex(), v.View());
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
           FormCanonicalComponentView(v.GridPointer(), v.UpperIndex());
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
    return CanonicalComponentView(v.GridPointer(), v.UpperIndex(), v.View());
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
  return CanonicalComponentViewUnary(v, [t](auto x) { return t * x; });
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
         FormCanonicalComponentView(v1.GridPointer(), v1.UpperIndex());
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
  using View = CanonicalComponentBase<Derived>;
  using T = std::ranges::range_value_t<View>;
  auto dPhi = std::ranges::views::repeat(v.GridPointer()->LongitudeSpacing(),
                                         v.GridPointer()->NumberOfLongitudes());
  auto dTheta = v.GridPointer()->CoLatitudeWeights();
  auto dArea = std::ranges::views::cartesian_product(dTheta, dPhi) |
               std::ranges::views::transform([](auto pair) {
                 return std::get<0>(pair) * std::get<1>(pair);
               });
  return std::ranges::fold_left(std::ranges::views::zip(dArea, v.View()), T{0},
                                [](auto acc, auto term) {
                                  auto [dArea, val] = term;
                                  return acc + val * dArea;
                                });
}

template <typename Derived>
auto Integrate(CanonicalComponentBase<Derived>&& v) {
  return Integrate(v);
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H
