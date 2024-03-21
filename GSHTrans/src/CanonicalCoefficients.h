#ifndef GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H
#define GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H

#include <FFTWpp/Core>
#include <algorithm>
#include <complex>
#include <iostream>
#include <memory>
#include <ranges>

#include "Concepts.h"
#include "GridBase.h"
#include "Indexing.h"

namespace GSHTrans {

//---------------------------------------------//
//             Define the base class           //
//---------------------------------------------//
template <typename Derived>
class CanonicalCoefficientBase {
  using Int = std::ptrdiff_t;

 public:
  // Data access methods.
  auto View() const { return _Derived().View(); }
  auto View() { return _Derived().View(); }
  auto begin() { return _Derived().View().begin(); }
  auto end() { return _Derived().View().end(); }

  auto& operator()(Int l, Int m) const {
    auto i = _Derived().Index(l, m);
    return _Derived()[i];
  }

  auto& operator()(Int l, Int m) {
    auto i = _Derived().Index(l, m);
    return _Derived()[i];
  }

  // Grid access methods.
  auto Grid() const { return _Derived().Grid(); }
  auto MaxDegree() const { return _Derived().MaxDegree(); }
  auto UpperIndex() const { return _Derived().UpperIndex(); }

 private:
  auto& _Derived() const { return static_cast<const Derived&>(*this); }
  auto& _Derived() { return static_cast<Derived&>(*this); }
};

//----------------------------------------------------//
//      Canonical coefficient storing its data        //
//----------------------------------------------------//

template <typename _Grid, RealOrComplexValued Type>
requires std::derived_from<_Grid, GridBase<_Grid>>
class CanonicalCoefficient
    : public CanonicalCoefficientBase<CanonicalCoefficient<_Grid, Type>>,
      public GSHIndices<std::conditional_t<std::same_as<Type, RealValued>,
                                           NonNegative, All>> {
  using Int = std::ptrdiff_t;
  using _GSHIndices = GSHIndices<
      std::conditional_t<std::same_as<Type, RealValued>, NonNegative, All>>;
  using Scalar = typename _Grid::complex_type;
  using Vector = FFTWpp::vector<Scalar>;

 public:
  using grid_type = _Grid;
  using real_or_complex_type = Type;
  using view_type = std::ranges::views::all_t<Vector>;

  using _GSHIndices::MaxDegree;
  using _GSHIndices::UpperIndex;

  // Methods required for CanonicalCoefficientBase.
  auto Grid() const { return _grid; }
  auto View() { return std::ranges::views::all(_data); }
  auto View() const { return std::ranges::views::all(_data); }
  auto operator[](Int i) const { return _data[i]; }
  auto& operator[](Int i) { return _data[i]; }

  // Constructors and assignement operators.
  CanonicalCoefficient() = default;

  CanonicalCoefficient(_Grid grid, Int n)
      : _GSHIndices(grid.MaxDegree(), grid.MaxDegree(), n),
        _grid{grid},
        _data{Vector(_GSHIndices::size())} {
    assert(
        std::ranges::contains(this->Grid().UpperIndices(), this->UpperIndex()));
  }

  CanonicalCoefficient(_Grid grid, Int n, Scalar s)
      : _GSHIndices(grid.MaxDegree(), grid.MaxDegree(), n),
        _grid{grid},
        _data{Vector(_GSHIndices::size(), s)} {
    assert(
        std::ranges::contains(this->Grid().UpperIndices(), this->UpperIndex()));
  }

  template <typename Function>
  requires std::regular_invocable<Function, Int, Int> &&
           std::convertible_to<std::invoke_result_t<Function, Int, Int>, Scalar>
  CanonicalCoefficient(_Grid grid, Int n, Function f)
      : CanonicalCoefficient(grid, n) {
    std::ranges::copy(this->Indices() | std::ranges::views::transform(
                                            [&f](auto index) -> Scalar {
                                              auto [l, m] = index;
                                              return f(l, m);
                                            }),
                      _data.begin());
  }

  CanonicalCoefficient(const CanonicalCoefficient&) = default;
  CanonicalCoefficient(CanonicalCoefficient&&) = default;

  template <typename Derived>
  CanonicalCoefficient(const CanonicalCoefficientBase<Derived>& other)
      : _GSHIndices(other.MaxDegree(), other.MaxDegree(), other.UpperIndex()),
        _grid{other.Grid()},
        _data{Vector(other.View().cbegin(), other.View().cend())} {}

  template <typename Derived>
  CanonicalCoefficient(CanonicalCoefficientBase<Derived>&& other)
      : CanonicalCoefficient(other) {}

  CanonicalCoefficient& operator=(const CanonicalCoefficient&) = default;
  CanonicalCoefficient& operator=(CanonicalCoefficient&&) = default;

  template <typename Derived>
  requires std::convertible_to<std::ranges::range_value_t<Derived>, Scalar>
  CanonicalCoefficient& operator=(
      const CanonicalCoefficientBase<Derived>& other) {
    assert(other.View().size() == this->size());
    std::ranges::copy(other.View() | std::ranges::views::transform(
                                         [](auto x) -> Scalar { return x; }),
                      _data.begin());
    return *this;
  }

  template <typename Derived>
  CanonicalCoefficient& operator=(CanonicalCoefficientBase<Derived>&& other) {
    *this = other;
    return *this;
  }

 private:
  _Grid _grid;
  Vector _data;
};

// Deduction guide for construction from another base type.
template <typename Derived>
CanonicalCoefficient(CanonicalCoefficientBase<Derived>& other)
    -> CanonicalCoefficient<typename Derived::grid_type,
                            typename Derived::real_or_complex_type>;

template <typename Derived>
CanonicalCoefficient(CanonicalCoefficientBase<Derived>&& other)
    -> CanonicalCoefficient<typename Derived::grid_type,
                            typename Derived::real_or_complex_type>;

// Type aliases for real and complex valued Coefficients.
template <typename _Grid>
using RealCanonicalCoefficient = CanonicalCoefficient<_Grid, RealValued>;

template <typename _Grid>
using ComplexCanonicalCoefficient = CanonicalCoefficient<_Grid, ComplexValued>;

//----------------------------------------------------------//
//      Canonical coefficient with a view to its data       //
//----------------------------------------------------------//
template <typename _Grid, RealOrComplexValued Type, std::ranges::view _View>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::ranges::input_range<_View>;
  // requires ComplexFloatingPoint<std::ranges::range_value_t<_View>>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::real_type>;
}
class CanonicalCoefficientView
    : public CanonicalCoefficientBase<
          CanonicalCoefficientView<_Grid, Type, _View>>,
      public std::ranges::view_interface<
          CanonicalCoefficientView<_Grid, Type, _View>>,
      public GSHIndices<std::conditional_t<std::same_as<Type, RealValued>,
                                           NonNegative, All>> {
  using Int = std::ptrdiff_t;
  using _GSHIndices = GSHIndices<
      std::conditional_t<std::same_as<Type, RealValued>, NonNegative, All>>;
  using Scalar = typename _Grid::complex_type;

 public:
  using grid_type = _Grid;
  using real_or_complex_type = Type;
  using view_type = _View;
  using std::ranges::view_interface<
      CanonicalCoefficientView<_Grid, Type, _View>>::size;

  using _GSHIndices::MaxDegree;
  using _GSHIndices::UpperIndex;

  // Methods required for CanonicalCoefficientBase.
  auto Grid() const { return _grid; }
  auto View() { return _view; }
  auto View() const { return _view; }
  auto operator[](Int i) const { return _view[i]; }
  auto& operator[](Int i) { return _view[i]; }

  CanonicalCoefficientView() = default;

  CanonicalCoefficientView(_Grid grid, Int n, _View view)
      : _GSHIndices(grid.MaxDegree(), grid.MaxDegree(), n),
        _grid{grid},
        _view{view} {
    assert(this->size() == _GSHIndices::size());
    assert(
        std::ranges::contains(this->Grid().UpperIndices(), this->UpperIndex()));
  }

  CanonicalCoefficientView(const CanonicalCoefficientView&) = default;
  CanonicalCoefficientView(CanonicalCoefficientView&&) = default;

  CanonicalCoefficientView& operator=(const CanonicalCoefficientView&) =
      default;
  CanonicalCoefficientView& operator=(CanonicalCoefficientView&&) = default;

  template <typename Derived>
  requires requires() {
    requires std::ranges::output_range<_View, Scalar>;
    requires std::convertible_to<std::ranges::range_value_t<Derived>, Scalar>;
  }
  CanonicalCoefficientView& operator=(
      const CanonicalCoefficientBase<Derived>& other) {
    assert(other.size() == this->size());
    assert(
        std::ranges::contains(this->Grid().UpperIndices(), other.UpperIndex()));
    std::ranges::copy(other.View() | std::ranges::views::transform(
                                         [](auto x) -> Scalar { return x; }),
                      _view.begin());

    return *this;
  }

  template <typename Derived>
  CanonicalCoefficientView& operator=(
      CanonicalCoefficientBase<Derived>&& other) {
    *this = other;
    return *this;
  }

 private:
  _Grid _grid;
  _View _view;
};

// Deduction guide to allow construction from ranges.
template <typename _Grid, RealOrComplexValued Type,
          std::ranges::viewable_range R>
CanonicalCoefficientView(_Grid, std::ptrdiff_t, R&&)
    -> CanonicalCoefficientView<_Grid, Type, std::ranges::views::all_t<R>>;

// Type aliases for real and complex valued Coefficient views.
template <typename _Grid, std::ranges::view _View>
using RealCanonicalCoefficientView =
    CanonicalCoefficientView<_Grid, RealValued, _View>;

template <typename _Grid, std::ranges::view _View>
using ComplexCanonicalCoefficientView =
    CanonicalCoefficientView<_Grid, ComplexValued, _View>;

// Range adaptors to form CanonicalCoefficientView from a view.
template <typename _Grid, RealOrComplexValued Type>
requires std::derived_from<_Grid, GridBase<_Grid>>
class FormCanonicalCoefficientView
    : public std::ranges::range_adaptor_closure<
          FormCanonicalCoefficientView<_Grid, Type>> {
  using Int = std::ptrdiff_t;

 public:
  FormCanonicalCoefficientView(_Grid grid, Int n)
      : _grid{std::move(grid)}, _n{n} {
    assert(std::ranges::contains(_grid.UpperIndices(), _n));
  }

  template <std::ranges::view View>
  auto operator()(View v) {
    return CanonicalCoefficientView<_Grid, Type, View>(_grid, _n, v);
  }

 private:
  _Grid _grid;
  Int _n;
};

// Type aliases for real and complex valued views.
template <typename _Grid>
using FormRealCanonicalCoefficientView =
    FormCanonicalCoefficientView<_Grid, RealValued>;

template <typename _Grid>
using FormComplexCanonicalCoefficientView =
    FormCanonicalCoefficientView<_Grid, ComplexValued>;

//--------------------------------------------------------//
//           Functions defined on the base class          //
//--------------------------------------------------------//

// Unary operations.

template <typename Derived, typename Function>
auto CanonicalCoefficientViewUnary(const CanonicalCoefficientBase<Derived>& v,
                                   Function f) {
  return v.View() | std::ranges::views::transform(f) |
         FormCanonicalCoefficientView<typename Derived::grid_type,
                                      typename Derived::real_or_complex_type>(
             v.Grid(), v.UpperIndex());
}

template <typename Derived>
auto operator-(CanonicalCoefficientBase<Derived>& v) {
  return CanonicalCoefficientViewUnary(v, [](auto x) { return -x; });
}

template <typename Derived>
auto operator-(CanonicalCoefficientBase<Derived>&& v) {
  return -v;
}

template <typename Derived>
auto real(CanonicalCoefficientBase<Derived>& v) {
  using S = std::ranges::range_value_t<CanonicalCoefficientBase<Derived>>;
  return CanonicalCoefficientViewUnary(v, [](auto x) { return std::real(x); });
}

template <typename Derived>
auto real(CanonicalCoefficientBase<Derived>&& v) {
  return real(v);
}

template <typename Derived>
auto imag(CanonicalCoefficientBase<Derived>& v) {
  using S = std::ranges::range_value_t<CanonicalCoefficientBase<Derived>>;
  return CanonicalCoefficientViewUnary(v, [](auto x) { return std::imag(x); });
}

template <typename Derived>
auto imag(CanonicalCoefficientBase<Derived>&& v) {
  return imag(v);
}

template <typename Derived>
auto conj(const CanonicalCoefficientBase<Derived>& v) {
  return CanonicalCoefficientViewUnary(v, [](auto x) { return std::conj(x); });
}

template <typename Derived>
auto conj(CanonicalCoefficientBase<Derived>&& v) {
  return conj(v);
}

template <typename Derived>
auto abs(CanonicalCoefficientBase<Derived>& v) {
  return CanonicalCoefficientViewUnary(v, [](auto x) { return std::abs(x); });
}

template <typename Derived>
auto abs(CanonicalCoefficientBase<Derived>&& v) {
  return abs(v);
}

// Unary operations involving a scalar.

template <typename Derived, Field S>
auto operator*(CanonicalCoefficientBase<Derived>& v, S s) {
  using T = std::ranges::range_value_t<CanonicalCoefficientBase<Derived>>;
  using R =
      std::conditional_t<ComplexFloatingPoint<T>, T,
                         std::conditional_t<ComplexFloatingPoint<S>, S, T>>;
  auto r = R(s);
  return CanonicalCoefficientViewUnary(v, [r](auto x) {
    std::cout << r << " " << x << std::endl;
    return r * x;
  });
}

template <typename Derived, Field S>
auto operator*(CanonicalCoefficientBase<Derived>&& v, S s) {
  return v * s;
}

template <typename Derived, Field S>
auto operator*(S s, CanonicalCoefficientBase<Derived>& v) {
  return v * s;
}

template <typename Derived, Field S>
auto operator*(S s, CanonicalCoefficientBase<Derived>&& v) {
  return v * s;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H
