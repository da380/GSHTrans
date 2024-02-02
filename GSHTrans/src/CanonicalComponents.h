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

template <std::ranges::viewable_range R, typename Grid>
CanonicalComponentView(R&&, std::shared_ptr<Grid>)
    -> CanonicalComponentView<std::ranges::views::all_t<R>, Grid>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COMPONENTS_GUARD_H
