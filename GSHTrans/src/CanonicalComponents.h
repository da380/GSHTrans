#ifndef GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H

#include <FFTWpp/All>
#include <algorithm>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/tuple/tuple.hpp>
#include <complex>
#include <functional>
#include <iostream>
#include <memory>
#include <ranges>
#include <type_traits>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

template <typename Derived>
class CanonicalComponentBase
    : public std::ranges::view_interface<CanonicalComponentBase<Derived>> {
  using Int = std::ptrdiff_t;

 public:
  auto DataView() { return _Derived()._DataView(); }
  auto DataView() const { return _Derived()._DataView(); }

  auto begin() { return DataView().begin(); }
  auto end() { return DataView().end(); }

  auto& Grid() const { return _Derived()._Grid(); }

  auto& operator()(Int iTheta, Int iPhi) const {
    assert(iTheta >= 0 && iTheta < Grid().NumberOfCoLatitudes());
    assert(iPhi >= 0 && iPhi < Grid().NumberOfLongitudes());
    auto i = Grid().NumberOfLongitudes() * iTheta + iPhi;
    return this->operator[](i);
  }

  auto& operator()(Int iTheta, Int iPhi) {
    assert(iTheta >= 0 && iTheta < Grid().NumberOfCoLatitudes());
    assert(iPhi >= 0 && iPhi < Grid().NumberOfLongitudes());
    auto i = Grid().NumberOfLongitudes() * iTheta + iPhi;
    return this->operator[](i);
  }

  template <typename Function>
  void Interpolate(Function&& f) {
    auto i = Int{0};
    for (auto theta : Grid().CoLatitudes()) {
      for (auto phi : Grid().Longitudes()) {
        this->operator[](i++) = f(theta, phi);
      }
    }
  }

 private:
  auto& _Derived() { return static_cast<Derived&>(*this); }
  const auto& _Derived() const { return static_cast<const Derived&>(*this); }
};

template <typename GSHGrid, RealOrComplexValued CoefficientType>
class CanonicalComponent : public CanonicalComponentBase<
                               CanonicalComponent<GSHGrid, CoefficientType>> {
  using Int = std::ptrdiff_t;
  using Real = typename GSHGrid::real_type;
  using Complex = typename GSHGrid::complex_type;
  using Scalar = std::conditional_t<std::same_as<CoefficientType, RealValued>,
                                    Real, Complex>;
  using Vector = FFTWpp::vector<Scalar>;

 public:
  using value_type = Scalar;
  using coefficient_type = CoefficientType;

  CanonicalComponent(GSHGrid& grid)
      : _grid{grid}, _data{Vector(_grid.ComponentSize())} {}

  CanonicalComponent(GSHGrid& grid, Scalar value)
      : _grid{grid}, _data{Vector(_grid.ComponentSize(), value)} {}

  template <typename Function>
  requires ScalarFunction2D<Function, Real, Scalar>
  CanonicalComponent(GSHGrid& grid, Function&& f) : CanonicalComponent(grid) {
    this->Interpolate(f);
  }

  CanonicalComponent(CanonicalComponent&) = default;
  CanonicalComponent(CanonicalComponent&&) = default;

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, value_type>
  CanonicalComponent(CanonicalComponentBase<Derived>& other)
      : _grid{other.Grid()}, _data{Vector(other.cbegin(), other.cend())} {}

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, value_type>
  CanonicalComponent(CanonicalComponentBase<Derived>&& other)
      : CanonicalComponent(other) {}

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, value_type>
  auto& operator=(CanonicalComponentBase<Derived>& other) {
    assert(other.size() == this->size());
    std::ranges::copy(other, this->begin());
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, value_type>
  auto& operator=(CanonicalComponentBase<Derived>&& other) {
    *this = other;
    return *this;
  }

 private:
  GSHGrid& _grid;
  Vector _data;

  auto _DataView() { return std::ranges::views::all(_data); }

  auto _DataView() const {
    return std::ranges::views::all(_data) | std::ranges::views::as_const;
  }

  auto& _Grid() const { return _grid; }

  friend class CanonicalComponentBase<
      CanonicalComponent<GSHGrid, CoefficientType>>;
};

template <typename GSHGrid, std::ranges::random_access_range View>
requires std::same_as<RemoveComplex<std::ranges::range_value_t<View>>,
                      typename GSHGrid::real_type>
class CanonicalComponentView
    : public CanonicalComponentBase<CanonicalComponentView<GSHGrid, View>> {
  using Int = std::ptrdiff_t;

 public:
  using value_type = std::ranges::range_value_t<View>;
  using coefficient_type = std::conditional_t<RealFloatingPoint<value_type>,
                                              RealValued, ComplexValued>;

  CanonicalComponentView(GSHGrid& grid, View view) : _grid{grid}, _view{view} {}

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, value_type>
  auto& operator=(CanonicalComponentBase<Derived>& other) {
    assert(other.size() == this->size());
    std::ranges::copy(other, this->begin());
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, value_type>
  auto& operator=(CanonicalComponentBase<Derived>&& other) {
    *this = other;
    return *this;
  }

 private:
  GSHGrid& _grid;
  View _view;

  auto _DataView() { return _view; }

  auto _DataView() const { return _view | std::ranges::views::as_const; }

  auto& _Grid() const { return _grid; }

  friend class CanonicalComponentBase<CanonicalComponentView<GSHGrid, View>>;
};

//-------------------------------------------------------------//
//       Functions to transform canonical component views      //
//-------------------------------------------------------------//

template <typename T>
class AffineTransformation {
 public:
  AffineTransformation() : _scale{1}, _shift{0} {}
  AffineTransformation(T scale, T shift) : _scale{scale}, _shift{shift} {}
  T operator()(T x) const { return _scale * x + _shift; }

 private:
  T _scale;
  T _shift;
};

template <typename Derived, typename Function>
auto CanonicalComponentUnaryOperation(CanonicalComponentBase<Derived>& view,
                                      Function f) {
  return CanonicalComponentView(
      view.Grid(), view.DataView() | std::ranges::views::transform(f));
}

template <typename Derived, typename Function>
auto CanonicalComponentUnaryOperation(CanonicalComponentBase<Derived>&& view,
                                      Function f) {
  return CanonicalComponentUnaryOperation(view, f);
}

// View of the Complex conjugation.
template <typename Derived>
auto conj(CanonicalComponentBase<Derived>& view) {
  using T = typename Derived::value_type;
  if constexpr (ComplexFloatingPoint<T>) {
    return CanonicalComponentUnaryOperation(
        view, [](auto x) { return std::conj(x); });
  } else {
    return CanonicalComponentView(view.Grid(), view.DataView());
  }
}

template <typename Derived>
auto conj(CanonicalComponentBase<Derived>&& view) {
  return conj(view);
}

// View of the real part.
template <typename Derived>
auto real(CanonicalComponentBase<Derived>& view) {
  using T = typename Derived::value_type;
  if constexpr (ComplexFloatingPoint<T>) {
    return CanonicalComponentUnaryOperation(
        view, [](auto x) { return std::real(x); });
  } else {
    return CanonicalComponentView(view.Grid(), view.DataView());
  }
}

template <typename Derived>
auto real(CanonicalComponentBase<Derived>&& view) {
  return real(view);
}

// View of the imaginary part.
template <typename Derived>
auto imag(CanonicalComponentBase<Derived>& view) {
  using T = typename Derived::value_type;
  if constexpr (ComplexFloatingPoint<T>) {
    return CanonicalComponentUnaryOperation(
        view, [](auto x) { return std::imag(x); });
  } else {
    return CanonicalComponentUnaryOperation(view,
                                            [](auto x) -> T { return 0; });
  }
}

template <typename Derived>
auto imag(CanonicalComponentBase<Derived>&& view) {
  return imag(view);
}

// View of the absolute value.
template <typename Derived>
auto abs(CanonicalComponentBase<Derived>& view) {
  return CanonicalComponentUnaryOperation(view,
                                          [](auto x) { return std::abs(x); });
}

template <typename Derived>
auto abs(CanonicalComponentBase<Derived>&& view) {
  return abs(view);
}

// Overloads for unary minus.
template <typename Derived>
auto operator-(CanonicalComponentBase<Derived>& view) {
  return CanonicalComponentUnaryOperation(view, [](auto x) { return -x; });
}

template <typename Derived>
auto operator-(CanonicalComponentBase<Derived>&& view) {
  return -view;
}

// Transform a pair of views using a binary operation.
template <typename Derived1, typename Derived2, typename Function>
requires std::convertible_to<typename Derived1::value_type,
                             typename Derived2::value_type>
auto CanonicalComponentBinaryOperation(CanonicalComponentBase<Derived1>& view1,
                                       CanonicalComponentBase<Derived2>& view2,
                                       Function f) {
  assert(view1.size() == view2.size());
  auto view =
      std::ranges::views::zip_transform(f, view1.DataView(), view2.DataView()) |
      std::ranges::views::as_const;
  return CanonicalComponentView(view1.Grid(), view);
}

template <typename Derived1, typename Derived2, typename Function>
requires std::convertible_to<typename Derived1::value_type,
                             typename Derived2::value_type>
auto CanonicalComponentBinaryOperation(CanonicalComponentBase<Derived1>&& view1,
                                       CanonicalComponentBase<Derived2>& view2,
                                       Function f) {
  return CanonicalComponentBinaryOperation(view1, view2, f);
}

template <typename Derived1, typename Derived2, typename Function>
requires std::convertible_to<typename Derived1::value_type,
                             typename Derived2::value_type>
auto CanonicalComponentBinaryOperation(CanonicalComponentBase<Derived1>& view1,
                                       CanonicalComponentBase<Derived2>&& view2,
                                       Function f) {
  return CanonicalComponentBinaryOperation(view1, view2, f);
}

template <typename Derived1, typename Derived2, typename Function>
requires std::convertible_to<typename Derived1::value_type,
                             typename Derived2::value_type>
auto CanonicalComponentBinaryOperation(CanonicalComponentBase<Derived1>&& view1,
                                       CanonicalComponentBase<Derived2>&& view2,
                                       Function f) {
  return CanonicalComponentBinaryOperation(view1, view2, f);
}

// Overloads for addition.
template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>& view2) {
  return CanonicalComponentBinaryOperation(view1, view2, std::plus<>());
}

template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>& view2) {
  return view1 + view2;
}

template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>&& view2) {
  return view1 + view2;
}

template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>&& view2) {
  return view1 + view2;
}

// Overloads for subtraction.
template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>& view2) {
  return CanonicalComponentBinaryOperation(view1, view2, std::minus<>());
}

template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>& view2) {
  return view1 - view2;
}

template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>&& view2) {
  return view1 - view2;
}

template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>&& view2) {
  return view1 - view2;
}

// Overloads for mulitplication.
template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>& view2) {
  return CanonicalComponentBinaryOperation(view1, view2, std::multiplies<>());
}

template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>& view2) {
  return view1 * view2;
}

template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>&& view2) {
  return view1 * view2;
}

template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>&& view2) {
  return view1 * view2;
}

template <typename GSHGrid, std::ranges::random_access_range View,
          typename Scalar>
requires std::same_as<RemoveComplex<std::ranges::range_value_t<View>>,
                      typename GSHGrid::real_type>
class CanonicalComponentAffineView
    : public CanonicalComponentBase<
          CanonicalComponentAffineView<GSHGrid, View, Scalar>> {
  using Int = std::ptrdiff_t;

 public:
  using value_type = std::ranges::range_value_t<View>;
  using coefficient_type = std::conditional_t<RealFloatingPoint<value_type>,
                                              RealValued, ComplexValued>;

  CanonicalComponentAffineView(GSHGrid& grid, View view, Scalar scale,
                               Scalar shift)
      : _grid{grid}, _view{view}, _scale{scale}, _shift{shift} {}

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, value_type>
  auto& operator=(CanonicalComponentBase<Derived>& other) {
    assert(other.size() == this->size());
    std::ranges::copy(other, this->begin());
    return *this;
  }

 private:
  GSHGrid& _grid;
  View _view;
  Scalar _scale;
  Scalar _shift;

  value_type _AffineTransform(value_type x) { return _scale * x + _shift; }

  auto _DataView() {
    using T = value_type;
    std::function<T(T)> f = _AffineTransform;

    // return _view | std::ranges::views::transform(f);

    return _view;
  }

  auto _DataView() const {
    return _view | std::ranges::views::transform(
                       [this](auto x) -> value_type { return x; });
  }

  auto& _Grid() const { return _grid; }

  friend class CanonicalComponentBase<
      CanonicalComponentAffineView<GSHGrid, View, Scalar>>;
};

// Overloads for scalar multiplication.
template <typename Derived, typename Scalar>
auto operator*(CanonicalComponentBase<Derived>& view, Scalar a) {
  return CanonicalComponentAffineView(view.Grid(), view.DataView(), a,
                                      Scalar{0});
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H
