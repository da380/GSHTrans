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
  auto begin() { return _Derived()._begin(); }
  auto end() { return _Derived()._end(); }

  auto GridPointer() const { return _Derived()._grid; }

  auto View() const { return _Derived()._View(); }
  auto View() { return _Derived()._View(); }

  auto NumberOfCoLatitudes() const {
    return GridPointer()->NumberOfCoLatitudes();
  }

  auto NumberOfLongitudes() const {
    return GridPointer()->NumberOfLongitudes();
  }

  auto Points() const { return GridPointer()->Points(); }

  auto operator()(Int iTheta, Int iPhi) const {
    auto i = iTheta * NumberOfLongitudes() + iPhi;
    return this->operator[](i);
  }

  auto& operator()(Int iTheta, Int iPhi) {
    auto i = iTheta * NumberOfLongitudes() + iPhi;
    return this->operator[](i);
  }

  template <typename Function>
  void Interpolate(Function&& f) {
    auto iter = Points().begin();
    for (auto& val : *this) {
      auto [theta, phi] = *iter++;
      val = f(theta, phi);
    }
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
  CanonicalComponent() = default;

  CanonicalComponent(std::shared_ptr<Grid> grid)
      : _grid{std::move(grid)}, _data{Vector(_grid->ComponentSize())} {}

  CanonicalComponent(const CanonicalComponent&) = default;
  CanonicalComponent(CanonicalComponent&&) = default;

  template <typename Derived>
  CanonicalComponent(CanonicalComponentBase<Derived>& other)
      : _grid{other.GridPointer()},
        _data{Vector(other.cbegin(), other.cend())} {}

  template <typename Derived>
  CanonicalComponent(CanonicalComponentBase<Derived>&& other)
      : CanonicalComponent(other) {}

  CanonicalComponent& operator=(const CanonicalComponent&) = default;
  CanonicalComponent& operator=(CanonicalComponent&&) = default;

  template <typename Derived>
  CanonicalComponent& operator=(CanonicalComponentBase<Derived>& other) {
    std::ranges::copy(other, _data.begin());
    return *this;
  }

  template <typename Derived>
  CanonicalComponent& operator=(CanonicalComponentBase<Derived>&& other) {
    *this = other;
    return *this;
  }

 private:
  std::shared_ptr<Grid> _grid;
  Vector _data;

  auto _begin() { return _data.begin(); }
  auto _end() { return _data.end(); }

  auto _View() { return std::ranges::views::all(_data); }
  auto _View() const {
    return std::ranges::views::all(_data) | std::ranges::views::as_const;
  }

  friend class CanonicalComponentBase<CanonicalComponent<Grid, Type>>;
};

//----------------------------------------------------------------//
//            Canonical component with view to its data           //
//----------------------------------------------------------------//

template <typename Grid, std::ranges::view View>
class CanonicalComponentView
    : public CanonicalComponentBase<CanonicalComponentView<Grid, View>> {
  using Int = std::ptrdiff_t;

 public:
  using view_type = View;

  CanonicalComponentView() = default;

  CanonicalComponentView(std::shared_ptr<Grid> grid, View view)
      : _grid{std::move(grid)}, _view{std::move(view)} {}

  CanonicalComponentView(const CanonicalComponentView&) = default;
  CanonicalComponentView(CanonicalComponentView&&) = default;

  CanonicalComponentView& operator=(const CanonicalComponentView&) = default;
  CanonicalComponentView& operator=(CanonicalComponentView&&) = default;

  template <typename Derived>
  CanonicalComponentView& operator=(CanonicalComponentBase<Derived>& other) {
    // std::ranges::copy(other, _view.begin());
    auto iter = _view.begin();
    for (auto val : other) *iter++ = val;
    return *this;
  }

  template <typename Derived>
  CanonicalComponentView& operator=(CanonicalComponentBase<Derived>&& other) {
    *this = other;
    return *this;
  }

 private:
  std::shared_ptr<Grid> _grid;
  View _view;

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

// Range adaptor to form CanonicalComponentView from a view.
template <typename Grid>
class PairWithGrid
    : public std::ranges::range_adaptor_closure<PairWithGrid<Grid>> {
 public:
  PairWithGrid(std::shared_ptr<Grid> grid) : _grid{std::move(grid)} {}

  template <std::ranges::view View>
  auto operator()(View v) {
    return CanonicalComponentView(_grid, v);
  }

 private:
  std::shared_ptr<Grid> _grid;
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
auto CanonicalComponentViewUnary(CanonicalComponentBase<Derived>& v,
                                 Function f) {
  return v.View() | std::ranges::views::transform(f) |
         PairWithGrid(v.GridPointer());
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
    return CanonicalComponentView(v.GridPointer(), v.View());
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
           PairWithGrid(v.GridPointer());
  } else {
    return CanonicalComponentViewUnary(v, [](auto x) { return std::imag(x); });
  }
}

template <typename Derived>
auto imag(CanonicalComponentBase<Derived>&& v) {
  return imag(v);
}

template <typename Derived>
auto conj(CanonicalComponentBase<Derived>& v) {
  using T = std::ranges::range_value_t<CanonicalComponentBase<Derived>>;
  if constexpr (RealFloatingPoint<T>) {
    return CanonicalComponentView(v.GridPointer(), v.View());
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
         PairWithGrid(v1.GridPointer());
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
