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

//---------------------------------------------------------------//
//                     View class definition                     //
//---------------------------------------------------------------//
template <std::ranges::view View, typename Grid>
requires requires() {
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<View>>,
                        typename Grid::real_type>;
}
class CanonicalComponentView
    : public std::ranges::view_interface<CanonicalComponentView<View, Grid>> {
  using Int = std::ptrdiff_t;

 public:
  using value_type = std::ranges::range_value_t<View>;
  using coefficient_type = std::conditional_t<RealFloatingPoint<value_type>,
                                              RealValued, ComplexValued>;

  CanonicalComponentView(View view, std::shared_ptr<Grid> grid)
      : _view{std::move(view)}, _grid{std::move(grid)} {
    assert(_grid->ComponentSize() == _view.size());
  }

  CanonicalComponentView(const CanonicalComponentView&) = default;
  CanonicalComponentView(CanonicalComponentView&&) = default;

  CanonicalComponentView& operator=(const CanonicalComponentView&) = default;
  CanonicalComponentView& operator=(CanonicalComponentView&&) = default;

  template <std::ranges::view OtherView>
  requires std::ranges::output_range<View,
                                     std::ranges::range_value_t<OtherView>>
  auto& operator=(CanonicalComponentView<OtherView, Grid> view) {
    assert(view.size() == _view.size());
    std::ranges::copy(view, _view.begin());
    return *this;
  }

  auto begin() { return _view.begin(); }
  auto end() { return _view.end(); }

  auto GridPointer() const { return _grid; }
  auto& GridReference() const { return _grid; }

  auto NumberOfCoLatitudes() const {
    return GridPointer()->NumberOfCoLatitudes();
  }

  auto NumberOfLongitudes() const {
    return GridPointer()->NumberOfLongitudes();
  }

 private:
  View _view;
  std::shared_ptr<Grid> _grid;
};

// Deduction guide to allow construction from ranges.
template <std::ranges::viewable_range R, typename Grid>
CanonicalComponentView(R&&, std::shared_ptr<Grid>)
    -> CanonicalComponentView<std::ranges::views::all_t<R>, Grid>;

namespace Views {
//------------------------------------------------------------//
//                    Range adaptor closure                   //
//------------------------------------------------------------//

template <typename Grid>
class CanonicalComponent
    : public std::ranges::range_adaptor_closure<CanonicalComponent<Grid>> {
 public:
  CanonicalComponent(std::shared_ptr<Grid> grid) : _grid{std::move(grid)} {}
  template <std::ranges::view View>
  auto operator()(View v) {
    return CanonicalComponentView(v, _grid);
  }

 private:
  std::shared_ptr<Grid> _grid;
};

}  // namespace Views

//--------------------------------------------------------------//
//          Operations on canonical compoenent views            //
//--------------------------------------------------------------//

template <Field S, Field T>
struct ReturnValueHelper {
  using type =
      std::conditional_t<ComplexFloatingPoint<T>, T,
                         std::conditional_t<ComplexFloatingPoint<S>, S, T>>;
};

template <Field S, Field T>
using ReturnValue = typename ReturnValueHelper<S, T>::type;

template <std::ranges::view View, typename Grid, typename Function>
auto CanonicalComponentViewUnary(CanonicalComponentView<View, Grid> v,
                                 Function f) {
  return v | std::ranges::views::transform(f) |
         Views::CanonicalComponent(v.GridPointer());
}

template <std::ranges::view View, typename Grid>
auto operator-(CanonicalComponentView<View, Grid> v) {
  return CanonicalComponentViewUnary(v, [](auto x) { return -x; });
}

template <std::ranges::view View, typename Grid>
auto real(CanonicalComponentView<View, Grid> v) {
  using T = std::ranges::range_value_t<View>;
  if constexpr (RealFloatingPoint<T>) {
    return v;
  } else {
    return CanonicalComponentViewUnary(v, [](auto x) { return std::real(x); });
  }
}

template <std::ranges::view View, typename Grid>
auto imag(CanonicalComponentView<View, Grid> v) {
  using T = std::ranges::range_value_t<View>;
  if constexpr (RealFloatingPoint<T>) {
    return std::ranges::views::repeat(T{0}, v.size()) |
           Views::CanonicalComponent(v.GridPointer());
  } else {
    return CanonicalComponentViewUnary(v, [](auto x) { return std::imag(x); });
  }
}

template <std::ranges::view View, typename Grid>
auto conj(CanonicalComponentView<View, Grid> v) {
  using T = std::ranges::range_value_t<View>;
  if constexpr (RealFloatingPoint<T>) {
    return v;
  } else {
    return CanonicalComponentViewUnary(v, [](auto x) { return std::conj(x); });
  }
}

template <std::ranges::view View, typename Grid>
auto abs(CanonicalComponentView<View, Grid> v) {
  return CanonicalComponentViewUnary(v, [](auto x) { return std::abs(x); });
}

template <std::ranges::view View, typename Grid, Field S>
auto operator*(CanonicalComponentView<View, Grid> v, S s) {
  using T = ReturnValue<S, std::ranges::range_value_t<View>>;
  auto t = T(s);
  return CanonicalComponentViewUnary(v, [t](auto x) { return t * x; });
}

template <std::ranges::view View, typename Grid, Field S>
auto operator*(S s, CanonicalComponentView<View, Grid> v) {
  return v * s;
}

template <std::ranges::view View, typename Grid, Field S>
auto operator+(CanonicalComponentView<View, Grid> v, S s) {
  using T = ReturnValue<S, std::ranges::range_value_t<View>>;
  auto t = T(s);
  return CanonicalComponentViewUnary(v, [t](auto x) { return t + x; });
}

template <std::ranges::view View, typename Grid, Field S>
auto operator+(S s, CanonicalComponentView<View, Grid> v) {
  return v + s;
}

template <std::ranges::view View, typename Grid, Field S>
auto operator-(CanonicalComponentView<View, Grid> v, S s) {
  using T = ReturnValue<S, std::ranges::range_value_t<View>>;
  auto t = T(s);
  return CanonicalComponentViewUnary(v, [t](auto x) { return x - t; });
}

template <std::ranges::view View, typename Grid, Field S>
auto operator/(CanonicalComponentView<View, Grid> v, S s) {
  using T = ReturnValue<S, std::ranges::range_value_t<View>>;
  auto t = T(s);
  return CanonicalComponentViewUnary(v, [t](auto x) { return x / t; });
}

template <std::ranges::view View, typename Grid, Field S>
auto pow(CanonicalComponentView<View, Grid> v, S s) {
  using T = ReturnValue<S, std::ranges::range_value_t<View>>;
  auto t = T(s);
  return CanonicalComponentViewUnary(v, [t](auto x) { return pow(x, t); });
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          typename Function>
auto CanonicalComponentViewBinary(CanonicalComponentView<View1, Grid> v1,
                                  CanonicalComponentView<View2, Grid> v2,
                                  Function f) {
  return std::ranges::zip_transform_view(f, v1, v2) |
         Views::CanonicalComponent(v1.GridPointer());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator+(CanonicalComponentView<View1, Grid> v1,
               CanonicalComponentView<View2, Grid> v2) {
  return CanonicalComponentViewBinary(v1, v2, std::plus<>());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator-(CanonicalComponentView<View1, Grid> v1,
               CanonicalComponentView<View2, Grid> v2) {
  return CanonicalComponentViewBinary(v1, v2, std::minus<>());
}

template <std::ranges::view View1, std::ranges::view View2, typename Grid>
auto operator*(CanonicalComponentView<View1, Grid> v1,
               CanonicalComponentView<View2, Grid> v2) {
  return CanonicalComponentViewBinary(v1, v2, std::multiplies<>());
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H
