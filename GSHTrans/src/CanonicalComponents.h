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
#include <numeric>
#include <ranges>
#include <type_traits>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

//----------------------------------------------------------//
//             Forward declarations of the classes          //
//----------------------------------------------------------//
template <typename Derived>
class CanonicalComponentBase;

template <typename GSHGrid, RealOrComplexFloatingPoint Scalar>
requires std::same_as < RemoveComplex<Scalar>,
typename GSHGrid::real_type > class CanonicalComponent;

template <typename GSHGrid, std::ranges::common_range View>
requires std::same_as < RemoveComplex<std::ranges::range_value_t<View>>,
typename GSHGrid::real_type > class CanonicalComponentView;

//----------------------------------------------------------//
//              Forward declare some functions              //
//----------------------------------------------------------//
template <typename Derived, typename Function>
auto CanonicalComponentUnaryOperation(CanonicalComponentBase<Derived>& view,
                                      Function f);

template <typename Derived, typename Function>
auto CanonicalComponentUnaryOperation(CanonicalComponentBase<Derived>&& view,
                                      Function f);

template <typename Derived>
auto conj(CanonicalComponentBase<Derived>& view);

template <typename Derived>
auto conj(CanonicalComponentBase<Derived>&& view);

template <typename Derived>
auto real(CanonicalComponentBase<Derived>& view);

template <typename Derived>
auto real(CanonicalComponentBase<Derived>&& view);

template <typename Derived>
auto imag(CanonicalComponentBase<Derived>& view);

template <typename Derived>
auto imag(CanonicalComponentBase<Derived>&& view);

template <typename Derived>
auto abs(CanonicalComponentBase<Derived>& view);

template <typename Derived>
auto abs(CanonicalComponentBase<Derived>&& view);

template <typename Derived>
auto operator-(CanonicalComponentBase<Derived>& view);

template <typename Derived>
auto operator-(CanonicalComponentBase<Derived>&& view);

template <typename Derived, typename S, typename Function>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto CanonicalComponentScalarOperation(CanonicalComponentBase<Derived>& view,
                                       S a, Function f);

template <typename Derived, typename S, typename Function>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto CanonicalComponentScalarOperation(CanonicalComponentBase<Derived>&& view,
                                       S a, Function f);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator*(CanonicalComponentBase<Derived>& view, S a);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator*(CanonicalComponentBase<Derived>&& view, S a);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator*(S a, CanonicalComponentBase<Derived>& view);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator*(S a, CanonicalComponentBase<Derived>&& view);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator/(CanonicalComponentBase<Derived>& view, S a);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator/(CanonicalComponentBase<Derived>&& view, S a);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator+(CanonicalComponentBase<Derived>& view, S a);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator+(CanonicalComponentBase<Derived>&& view, S a);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator+(S a, CanonicalComponentBase<Derived>& view);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator+(S a, CanonicalComponentBase<Derived>&& view);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator-(CanonicalComponentBase<Derived>& view, S a);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator-(CanonicalComponentBase<Derived>&& view, S a);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator-(S a, CanonicalComponentBase<Derived>& view);

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator-(S a, CanonicalComponentBase<Derived>&& view);

template <typename Derived1, typename Derived2, typename Function>
auto CanonicalComponentBinaryOperation(CanonicalComponentBase<Derived1>& view1,
                                       CanonicalComponentBase<Derived2>& view2,
                                       Function f);

template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>& view2);

template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>& view2);

template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>&& view2);

template <typename Derived1, typename Derived2>
auto operator+(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>&& view2);

template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>& view2);

template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>& view2);

template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>&& view2);

template <typename Derived1, typename Derived2>
auto operator-(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>&& view2);

template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>& view2);

template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>& view2);

template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>& view1,
               CanonicalComponentBase<Derived2>&& view2);

template <typename Derived1, typename Derived2>
auto operator*(CanonicalComponentBase<Derived1>&& view1,
               CanonicalComponentBase<Derived2>&& view2);

// Convert range into view within the Boost Framwork.
template <std::ranges::common_range Range>
auto MakeView(Range&& range) {
  return boost::sub_range<Range>(range);
}

//-------------------------------------------------------//
//          Base class for canonical components          //
//-------------------------------------------------------//
template <typename Derived>
class CanonicalComponentBase {
  using Int = std::ptrdiff_t;

 public:
  auto DataView() { return GetDerived()._DataView(); }
  auto& Grid() const { return GetDerived()._Grid(); }

  auto begin() { return DataView().begin(); }
  auto end() { return DataView().end(); }
  auto size() { return DataView().size(); }

  auto& operator[](Int i) const {
    assert(i >= 0 && i < size());
    return DataView()[i];
  }

  auto& operator[](Int i) {
    assert(i >= 0 && i < size());
    return DataView()[i];
  }

  auto& operator()(Int iTheta, Int iPhi) const {
    assert(iTheta >= 0 && iTheta < Grid().NumberOfCoLatitudes());
    assert(iPhi >= 0 && iPhi < Grid().NumberOfLongitudes());
    auto i = Grid().NumberOfLongitudes() * iTheta + iPhi;
    return operator[](i);
  }

  auto& operator()(Int iTheta, Int iPhi) {
    assert(iTheta >= 0 && iTheta < Grid().NumberOfCoLatitudes());
    assert(iPhi >= 0 && iPhi < Grid().NumberOfLongitudes());
    auto i = Grid().NumberOfLongitudes() * iTheta + iPhi;
    return operator[](i);
  }

  template <typename Function>
  void Interpolate(Function&& f) {
    auto i = Int{0};
    for (auto theta : Grid().CoLatitudes()) {
      for (auto phi : Grid().Longitudes()) {
        operator[](i++) = f(theta, phi);
      }
    }
  }

  auto Integrate() { return Grid().Integrate(DataView()); }

  auto L2Norm() {
    auto f = *this * conj(*this);
    return std::sqrt(real(f).Integrate());
  }

 private:
  auto& GetDerived() { return static_cast<Derived&>(*this); }
  auto& GetDerived() const { return static_cast<const Derived&>(*this); }
};

//--------------------------------------------------------//
//      Canonical component class that owns its data      //
//--------------------------------------------------------//
template <typename GSHGrid, RealOrComplexFloatingPoint Scalar>
requires std::same_as < RemoveComplex<Scalar>,
typename GSHGrid::real_type > class CanonicalComponent
    : public CanonicalComponentBase<CanonicalComponent<GSHGrid, Scalar>> {
  using Int = std::ptrdiff_t;
  using Vector = FFTWpp::vector<Scalar>;

 public:
  using value_type = Scalar;

  CanonicalComponent(GSHGrid& grid)
      : _grid{grid}, _data{Vector(grid.ComponentSize())} {}

  CanonicalComponent(GSHGrid& grid, Scalar value)
      : _grid{grid}, _data{Vector(grid.ComponentSize(), value)} {}

  CanonicalComponent(CanonicalComponent&) = default;

  CanonicalComponent(CanonicalComponent&&) = default;

  template <typename Derived>
  CanonicalComponent(CanonicalComponentBase<Derived>& other)
      : _grid{other.Grid()}, _data{Vector(other.begin(), other.end())} {}

  template <typename Derived>
  CanonicalComponent(CanonicalComponentBase<Derived>&& other)
      : _grid{other.Grid()}, _data{Vector(other.begin(), other.end())} {}

  template <typename Derived>
  auto& operator=(CanonicalComponentBase<Derived>& other) {
    _grid = other.Grid();
    _data = Vector(other.begin(), other.end());
    return *this;
  }

  template <typename Derived>
  auto& operator=(CanonicalComponentBase<Derived>&& other) {
    _grid = other.Grid();
    _data = Vector(other.begin(), other.end());
    return *this;
  }

 private:
  Vector _data;
  GSHGrid& _grid;

  auto _DataView() { return MakeView(_data); }
  auto& _Grid() const { return _grid; }

  friend class CanonicalComponentBase<CanonicalComponent<GSHGrid, Scalar>>;
};

//-----------------------------------------------------------------//
//    Canonical component class that stores a view to its data     //
//-----------------------------------------------------------------//
template <typename GSHGrid, std::ranges::common_range View>
requires std::same_as < RemoveComplex<std::ranges::range_value_t<View>>,
typename GSHGrid::real_type > class CanonicalComponentView
    : public CanonicalComponentBase<CanonicalComponentView<GSHGrid, View>> {
  using Int = std::ptrdiff_t;
  using Scalar = std::ranges::range_value_t<View>;

 public:
  using value_type = Scalar;

  CanonicalComponentView(GSHGrid& grid, View view) : _grid{grid}, _view{view} {
    assert(_grid.ComponentSize() == view.size());
  }

  template <typename Derived>
  auto operator=(CanonicalComponentBase<Derived>& other) {
    assert(other.size() == _view.size());
    std::ranges::copy(other.DataView(), _view.begin());
    return *this;
  }

  template <typename Derived>
  auto operator=(CanonicalComponentBase<Derived>&& other) {
    assert(other.size() == _view.size());
    std::ranges::copy(other.DataView(), _view.begin());
    return *this;
  }

 private:
  GSHGrid& _grid;
  View _view;

  auto _DataView() { return _view; }
  auto& _Grid() const { return _grid; }

  friend class CanonicalComponentBase<CanonicalComponentView<GSHGrid, View>>;
};

//-------------------------------------------------------------//
//   Define functions to transform canonical component views   //
//-------------------------------------------------------------//

// Transform a view using unary operation.
template <typename Derived, typename Function>
auto CanonicalComponentUnaryOperation(CanonicalComponentBase<Derived>& view,
                                      Function f) {
  return CanonicalComponentView(
      view.Grid(), view.DataView() | boost::adaptors::transformed(f));
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
    return view;
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
    return view;
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
    return CanonicalComponentUnaryOperation(view, [](auto x) { return 0; });
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

// View formed from component-wise scalar operator
template <typename Derived, typename S, typename Function>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto CanonicalComponentScalarOperation(CanonicalComponentBase<Derived>& view,
                                       S a, Function f) {
  using T = typename Derived::value_type;
  auto g = std::function<T(T)>([a, f](T x) { return f(x, static_cast<T>(a)); });
  return CanonicalComponentUnaryOperation(view, g);
}

template <typename Derived, typename S, typename Function>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto CanonicalComponentScalarOperation(CanonicalComponentBase<Derived>&& view,
                                       S a, Function f) {
  return CanonicalComponentScalarOperation(view, a, f);
}

// Overloads for scalar multiplication.
template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator*(CanonicalComponentBase<Derived>& view, S a) {
  return CanonicalComponentScalarOperation(view, a, std::multiplies<>());
}

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator*(CanonicalComponentBase<Derived>&& view, S a) { return view * a; }

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator*(S a, CanonicalComponentBase<Derived>& view) { return view * a; }

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator*(S a, CanonicalComponentBase<Derived>&& view) { return view * a; }

// Overloads for scalar division.
template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator/(CanonicalComponentBase<Derived>& view, S a) {
  return CanonicalComponentScalarOperation(view, a, std::divides<>());
}

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator/(CanonicalComponentBase<Derived>&& view, S a) { return view / a; }

// Overloads for scalar addition.
template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator+(CanonicalComponentBase<Derived>& view, S a) {
  return CanonicalComponentScalarOperation(view, a, std::plus<>());
}

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator+(CanonicalComponentBase<Derived>&& view, S a) { return view + a; }

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator+(S a, CanonicalComponentBase<Derived>& view) { return view + a; }

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator+(S a, CanonicalComponentBase<Derived>&& view) { return view + a; }

// Overloads for scalar subtraction.
template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator-(CanonicalComponentBase<Derived>& view, S a) {
  return CanonicalComponentScalarOperation(view, a, std::minus<>());
}

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator-(CanonicalComponentBase<Derived>&& view, S a) { return view - a; }

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator-(S a, CanonicalComponentBase<Derived>& view) { return view - a; }

template <typename Derived, typename S>
requires requires {
  requires std::convertible_to<S, typename Derived::value_type>;
}
auto operator-(S a, CanonicalComponentBase<Derived>&& view) { return view - a; }

// Transform a pair of views using a binary operation.
template <typename Derived1, typename Derived2, typename Function>
auto CanonicalComponentBinaryOperation(CanonicalComponentBase<Derived1>& view1,
                                       CanonicalComponentBase<Derived2>& view2,
                                       Function f) {
  auto view = boost::combine(view1.DataView(), view2.DataView()) |
              boost::adaptors::transformed([f](auto pair) {
                auto x = boost::get<0>(pair);
                auto y = boost::get<1>(pair);
                return f(x, y);
              });
  return CanonicalComponentView(view1.Grid(), view);
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

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H
