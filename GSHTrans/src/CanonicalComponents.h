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

// Convert range into view within the Boost Framwork.
template <std::ranges::common_range Range>
auto MakeView(Range&& range) {
  return boost::sub_range<Range>(range);
}

//----------------------------------------------------------//
//             Forward declarations of the classes          //
//----------------------------------------------------------//

template <typename Derived>
class CanonicalComponentBase;

template <typename GSHGrid>
class RealCanonicalComponent;

template <typename GSHGrid>
class ComplexCanonicalComponent;

template <typename GSHGrid, std::ranges::common_range View>
requires std::same_as < RemoveComplex<std::ranges::range_value_t<View>>,
typename GSHGrid::real_type > class CanonicalComponentView;

//----------------------------------------------------------//
//              Forward declare some functions              //
//----------------------------------------------------------//

template <typename Derived>
auto conj(CanonicalComponentBase<Derived>& view);

template <typename Derived>
auto conj(CanonicalComponentBase<Derived>&& view);

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
//               Real canonical component class           //
//--------------------------------------------------------//

template <typename GSHGrid>
class RealCanonicalComponent
    : public CanonicalComponentBase<RealCanonicalComponent<GSHGrid>> {
  using Int = std::ptrdiff_t;
  using Real = typename GSHGrid::real_type;
  using Vector = FFTWpp::vector<Real>;

 public:
  using value_type = Real;

  RealCanonicalComponent(GSHGrid& grid)
      : _grid{grid}, _data{Vector(grid.ComponentSize())} {}

  RealCanonicalComponent(GSHGrid& grid, Real value)
      : _grid{grid}, _data{Vector(grid.ComponentSize(), value)} {}

  template <typename Function>
  requires ScalarFunction2D<Function, Real, Real> RealCanonicalComponent(
      GSHGrid& grid, Function&& f)
      : RealCanonicalComponent(grid) {
    this->Interpolate(f);
  }

  RealCanonicalComponent(RealCanonicalComponent&) = default;

  RealCanonicalComponent(RealCanonicalComponent&&) = default;

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, Real>
  RealCanonicalComponent(CanonicalComponentBase<Derived>& other)
      : _grid{other.Grid()}, _data{Vector(other.begin(), other.end())} {}

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, Real>
  RealCanonicalComponent(CanonicalComponentBase<Derived>&& other)
      : _grid{other.Grid()}, _data{Vector(other.begin(), other.end())} {}

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, Real>
  auto& operator=(CanonicalComponentBase<Derived>& other) {
    _grid = other.Grid();
    _data = Vector(other.begin(), other.end());
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, Real>
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

  friend class CanonicalComponentBase<RealCanonicalComponent<GSHGrid>>;
};

//-----------------------------------------------------------//
//               Complex canonical component class           //
//-----------------------------------------------------------//

template <typename GSHGrid>
class ComplexCanonicalComponent
    : public CanonicalComponentBase<ComplexCanonicalComponent<GSHGrid>> {
  using Int = std::ptrdiff_t;
  using Real = typename GSHGrid::real_type;
  using Complex = typename GSHGrid::complex_type;
  using Vector = FFTWpp::vector<Complex>;

 public:
  using value_type = Complex;

  ComplexCanonicalComponent(GSHGrid& grid)
      : _grid{grid}, _data{Vector(grid.ComponentSize())} {}

  ComplexCanonicalComponent(GSHGrid& grid, Complex value)
      : _grid{grid}, _data{Vector(grid.ComponentSize(), value)} {}

  template <typename Function>
  requires ScalarFunction2D<Function, Real, Complex> ComplexCanonicalComponent(
      GSHGrid& grid, Function&& f)
      : ComplexCanonicalComponent(grid) {
    this->Interpolate(f);
  }

  ComplexCanonicalComponent(ComplexCanonicalComponent&) = default;

  ComplexCanonicalComponent(ComplexCanonicalComponent&&) = default;

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, Complex>
  ComplexCanonicalComponent(CanonicalComponentBase<Derived>& other)
      : _grid{other.Grid()}, _data{Vector(other.begin(), other.end())} {}

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, Complex>
  ComplexCanonicalComponent(CanonicalComponentBase<Derived>&& other)
      : _grid{other.Grid()}, _data{Vector(other.begin(), other.end())} {}

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, Complex>
  auto& operator=(CanonicalComponentBase<Derived>& other) {
    _grid = other.Grid();
    _data = Vector(other.begin(), other.end());
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::value_type, Complex>
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

  friend class CanonicalComponentBase<ComplexCanonicalComponent<GSHGrid>>;
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
requires std::integral<S> or RealOrComplexFloatingPoint<S>
auto CanonicalComponentScalarOperation(CanonicalComponentBase<Derived>& view,
                                       S a, Function f) {
  using T = typename Derived::value_type;
  if constexpr (RealFloatingPoint<S> or std::integral<S>) {
    auto b = T(a);
    auto g = std::function<T(T)>([b, f](T x) { return f(x, b); });
    return CanonicalComponentUnaryOperation(view, g);
  } else {
    auto g = std::function<S(T)>([a, f](T x) { return f(S(x), a); });
    return CanonicalComponentUnaryOperation(view, g);
  }
}

template <typename Derived, typename S, typename Function>
auto CanonicalComponentScalarOperation(CanonicalComponentBase<Derived>&& view,
                                       S a, Function f) {
  return CanonicalComponentScalarOperation(view, a, f);
}

// Overloads for scalar multiplication.
template <typename Derived, typename S>
requires std::integral<S> or RealOrComplexFloatingPoint<S>
auto operator*(CanonicalComponentBase<Derived>& view, S a) {
  return CanonicalComponentScalarOperation(view, a, std::multiplies<>());
}

template <typename Derived, typename S>
auto operator*(CanonicalComponentBase<Derived>&& view, S a) {
  return view * a;
}

template <typename Derived, typename S>
auto operator*(S a, CanonicalComponentBase<Derived>& view) {
  return view * a;
}

template <typename Derived, typename S>
auto operator*(S a, CanonicalComponentBase<Derived>&& view) {
  return view * a;
}

// Overloads for scalar division.
template <typename Derived, typename S>
requires std::integral<S> or RealOrComplexFloatingPoint<S>
auto operator/(CanonicalComponentBase<Derived>& view, S a) {
  return CanonicalComponentScalarOperation(view, a, std::divides<>());
}

template <typename Derived, typename S>
auto operator/(CanonicalComponentBase<Derived>&& view, S a) {
  return view / a;
}

// Overloads for scalar addition.
template <typename Derived, typename S>
requires std::integral<S> or RealOrComplexFloatingPoint<S>
auto operator+(CanonicalComponentBase<Derived>& view, S a) {
  return CanonicalComponentScalarOperation(view, a, std::plus<>());
}

template <typename Derived, typename S>
auto operator+(CanonicalComponentBase<Derived>&& view, S a) {
  return view + a;
}

template <typename Derived, typename S>
auto operator+(S a, CanonicalComponentBase<Derived>& view) {
  return view + a;
}

template <typename Derived, typename S>
auto operator+(S a, CanonicalComponentBase<Derived>&& view) {
  return view + a;
}

// Overloads for scalar subtraction.
template <typename Derived, typename S>
requires std::integral<S> or RealOrComplexFloatingPoint<S>
auto operator-(CanonicalComponentBase<Derived>& view, S a) {
  return CanonicalComponentScalarOperation(view, a, std::minus<>());
}

template <typename Derived, typename S>
auto operator-(CanonicalComponentBase<Derived>&& view, S a) {
  return view - a;
}

template <typename Derived, typename S>
auto operator-(S a, CanonicalComponentBase<Derived>& view) {
  return view - a;
}

template <typename Derived, typename S>
auto operator-(S a, CanonicalComponentBase<Derived>&& view) {
  return view - a;
}

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
