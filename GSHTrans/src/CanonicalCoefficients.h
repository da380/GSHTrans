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
    return this->operator[](i);
  }

  auto& operator()(Int l, Int m) {
    auto i = _Derived().Index(l, m);
    return this->operator[](i);
  }

  // Grid access methods.
  auto Grid() const { return _Derived().Grid(); }

 private:
  auto& _Derived() const { return static_cast<const Derived&>(*this); }
  auto& _Derived() { return static_cast<Derived&>(*this); }
};

//----------------------------------------------------//
//      Canonical coefficient storing its data        //
//----------------------------------------------------//
template <typename GSHGrid, RealOrComplexValued Type>
requires std::derived_from<GSHGrid, GridBase<GSHGrid>>
class CanonicalCoefficient
    : public CanonicalCoefficientBase<CanonicalCoefficient<GSHGrid, Type>>,
      public GSHIndices<std::conditional_t<std::same_as<Type, RealValued>,
                                           NonNegative, All>> {
  using Int = std::ptrdiff_t;
  using Indices = GSHIndices<
      std::conditional_t<std::same_as<Type, RealValued>, NonNegative, All>>;
  using Scalar = typename GSHGrid::complex_type;
  using Vector = FFTWpp::vector<Scalar>;

 public:
  using grid_type = GSHGrid;
  using view_type = std::ranges::views::all_t<Vector>;

  // Methods required for CanonicalCoefficientBase.
  auto Grid() const { return _grid; }
  auto View() { return std::ranges::views::all(_data); }
  auto View() const { return std::ranges::views::all(_data); }

  // Constructors and assignement operators.
  CanonicalCoefficient() = default;

  CanonicalCoefficient(GSHGrid grid, Int n)
      : Indices(grid.MaxDegree(), grid.MaxDegree(), n),
        _grid{grid},
        _data{Vector(Indices::size())} {}

 private:
  GSHGrid _grid;
  Vector _data;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H
