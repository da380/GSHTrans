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

//----------------------------------------------------------------//
//                    Base class to set interface                 //
//----------------------------------------------------------------//
template <typename Derived>
class CanonicalComponentBase {
  using Int = std::ptrdiff_t;

 public:
  // Return grid information.
  auto GridPointer() const { return _Derived()._grid; }
  auto& GridReference() const { return *_Derived()._grid; }

  auto NumberOfCoLatitudes() const {
    return GridPointer()->NumberOfCoLatitudes();
  }

  auto NumberOfLongitudes() const {
    return GridPointer()->NumberOfLongitudes();
  }

 private:
  auto _Derived() { return static_cast<Derived&>(*this); }
  auto _Derived() const { return static_cast<const Derived&>(*this); }
};

//---------------------------------------------------------------//
//                     View class definition                     //
//---------------------------------------------------------------//
template <std::ranges::view View, typename Grid>
requires requires() {
  requires std::same_as<std::ranges::range_value_t<View>,
                        typename Grid::real_type>;
}
class CanonicalComponentView
    : public CanonicalComponentBase<CanonicalComponentView<View, Grid>>,
      public std::ranges::view_interface<CanonicalComponentView<View, Grid>> {
  using Int = std::ptrdiff_t;

 public:
  using coefficient_type =
      std::conditional_t<RealFloatingPoint<std::ranges::range_value_t<View>>,
                         RealValued, ComplexValued>;

  CanonicalComponentView(View view, std::shared_ptr<Grid> grid)
      : _view{std::move(view)}, _grid{std::move(grid)} {
    assert(_grid->ComponentSize() == _view.size());
  }

  auto begin() { return _view.begin(); }
  auto end() { return _view.end(); }

 private:
  View _view;
  std::shared_ptr<Grid> _grid;

  friend class CanonicalComponentBase<CanonicalComponentView<View, Grid>>;
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

template <std::ranges::view View, typename Grid, Field S>
auto operator*(CanonicalComponentView<View, Grid> v, S s) {
  return CanonicalComponentViewUnary(v, [s](auto x) { return s * x; });
}

template <std::ranges::view View, typename Grid, Field S>
auto operator*(S s, CanonicalComponentView<View, Grid> v) {
  return v * s;
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
