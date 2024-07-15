#ifndef GSH_TRANS_SCALAR_FIELD_COEFFICIENT_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_COEFFICIENT_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../CoefficientBase.h"
#include "../Concepts.h"
#include "../GridBase.h"
#include "ScalarFieldCoefficientBase.h"

namespace GSHTrans {

// Forward declare class.
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ScalarFieldCoefficient;

// Set traits.
namespace Internal {

template <typename _Grid, RealOrComplexValued _Value>
struct Traits<ScalarFieldCoefficient<_Grid, _Value>> {
  using Int = std::ptrdiff_t;
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Writeable = std::true_type;
};

}  // namespace Internal

template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ScalarFieldCoefficient
    : public ScalarFieldCoefficientBase<ScalarFieldCoefficient<_Grid, _Value>> {
 public:
  using Grid =
      typename Internal::Traits<ScalarFieldCoefficient<_Grid, _Value>>::Grid;
  using Value =
      typename Internal::Traits<ScalarFieldCoefficient<_Grid, _Value>>::Value;
  using Int =
      typename Internal::Traits<ScalarFieldCoefficient<_Grid, _Value>>::Int;
  using Real =
      typename Internal::Traits<ScalarFieldCoefficient<_Grid, _Value>>::Real;
  using Complex =
      typename Internal::Traits<ScalarFieldCoefficient<_Grid, _Value>>::Complex;
  using Writeable = typename Internal::Traits<
      ScalarFieldCoefficient<_Grid, _Value>>::Writeable;

  // Return the grid.
  auto GetGrid() const { return _grid; }

  // Read access to data.
  auto operator[](Int l, Int m) const { return _data[this->Index(l, m)]; }

  // Write access to data.
  auto& operator[](Int l, Int m)
  requires Writeable::value
  {
    return _data[this->Index(l, m)];
  }
  // Default constructor.
  ScalarFieldCoefficient() = default;

  // Construct from grid initialising values to zero.
  ScalarFieldCoefficient(_Grid grid)
      : _grid{grid}, _data{FFTWpp::vector<Complex>(this->CoefficientSize())} {}

  // Construction from grid initialising values with a function.
  template <typename Function>
  requires ScalarValuedFunction<Function, Real, Value>
  ScalarFieldCoefficient(_Grid grid, Function&& f)
      : ScalarFieldCoefficient(grid) {
    std::ranges::copy(_grid.InterpolateFunction(f), _data.begin());
  }

  // Construct from an element of the base class.
  template <typename Derived>
  requires std::convertible_to<typename Derived::Complex, Complex>
  ScalarFieldCoefficient(const ScalarFieldCoefficientBase<Derived>& other)
      : ScalarFieldCoefficient(other.GetGrid()) {
    for (auto [l, m] : this->Indices()) {
      operator[](l, m) = other[l, m];
    }
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Complex, Complex>
  ScalarFieldCoefficient(ScalarFieldCoefficientBase<Derived>&& other)
      : ScalarFieldCoefficient(other) {}

  // Default copy and move constructors.
  ScalarFieldCoefficient(const ScalarFieldCoefficient&) = default;
  ScalarFieldCoefficient(ScalarFieldCoefficient&&) = default;

  // Default copy and move assigment.
  ScalarFieldCoefficient& operator=(const ScalarFieldCoefficient&) = default;
  ScalarFieldCoefficient& operator=(ScalarFieldCoefficient&&) = default;

  // Use assignment defined in base class.
  using ScalarFieldCoefficientBase<
      ScalarFieldCoefficient<_Grid, _Value>>::operator=;

  // Return view to the data.
  auto Data() { return std::ranges::views::all(_data); }

 private:
  _Grid _grid;
  FFTWpp::vector<Complex> _data;
};

// Type aliases for real and complex fields.
template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealScalarFieldCoefficient = ScalarFieldCoefficient<Grid, RealValued>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexScalarFieldCoefficient =
    ScalarFieldCoefficient<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif