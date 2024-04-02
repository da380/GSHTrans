#ifndef GSH_TRANS_VECTOR_FIELD_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <algorithm>
#include <array>
#include <complex>
#include <ranges>
#include <utility>
#include <vector>

#include "Concepts.h"
#include "GridBase.h"

namespace GSHTrans {

template <std::ranges::view _View, typename _Grid, RealOrComplexValued _Value>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
  requires(std::same_as<_Value, ComplexValued> &&
           ComplexFloatingPoint<std::ranges::range_value_t<_View>>) ||
              (std::same_as<_Value, RealValued> &&
               RealFloatingPoint<std::ranges::range_value_t<_View>>);
  requires(std::same_as<_Value, ComplexValued> &&
           std::same_as<typename _Grid::MRange, All>) ||
              std::same_as<_Value, RealValued>;
}
class VectorField {
  using Int = std::ptrdiff_t;
  using _Field = Field<_View, _Grid, _Value>;
  using _Data = std::array<_Field, 3>;

 public:
  // Public type aliases.
  using Grid = _Grid;
  using View = _View;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;

  VectorField() = default;

  // Construct from field components.
  VectorField(_Field u0, _Field u1, _Field u2) : _data{_Data{u0, u1, u2}} {}

  // Construct from views to the data for each component.
  VectorField(_View v0, _View v1, _View v2, _Grid grid)
      : VectorField(_Field(v0, grid), _Field(v1, grid), _Field(v2, grid)) {}

  auto operator()(Int alpha) { return _data[alpha + 1]; }

 private:
  _Data _data;
};

// Deduction guides for construction from ranges.
template <std::ranges::viewable_range R, typename Grid,
          RealOrComplexValued Value>
VectorField(R&&, R&&, R&&, Grid)
    -> VectorField<std::ranges::views::all_t<R>, Grid, Value>;

template <std::ranges::view View, typename Grid, RealOrComplexValued Value>
VectorField(View, Grid)
    -> VectorField<std::ranges::subrange<std::ranges::iterator_t<View>>, Grid,
                   Value>;

// Type aliases for real and complex vector fields.
template <std::ranges::view View, typename Grid>
using RealVectorField = VectorField<View, Grid, RealValued>;

template <std::ranges::view View, typename Grid>
using ComplexVectorField = VectorField<View, Grid, ComplexValued>;

template <std::ranges::view View1, std::ranges::view View2, typename Grid,
          RealOrComplexValued Value>
auto operator+(VectorField<View1, Grid, Value> u1,
               VectorField<View2, Grid, Value> u2) {
  using Field = decltype(u1(-1) + u2(-1));
  return VectorField<typename Field::View, Grid, Value>(
      u1(-1) + u2(-1), u1(0) + u2(0), u1(1) + u2(1));
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_VECTOR_FIELD_GUARD_H