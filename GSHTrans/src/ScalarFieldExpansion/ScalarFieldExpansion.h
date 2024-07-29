#ifndef GSH_TRANS_SCALAR_FIELD_EXPANSION_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_EXPANSION_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../ExpansionBase.h"
#include "../GridBase.h"
#include "ScalarFieldExpansionBase.h"

namespace GSHTrans {

// Forward declare class.
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ScalarFieldExpansion;

// Set traits.
namespace Internal {

template <typename _Grid, RealOrComplexValued _Value>
struct Traits<ScalarFieldExpansion<_Grid, _Value>> {
  using Int = std::ptrdiff_t;
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<Value, RealValued>, Real, Complex>;
  using Writeable = std::true_type;
};

}  // namespace Internal

template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ScalarFieldExpansion
    : public ScalarFieldExpansionBase<ScalarFieldExpansion<_Grid, _Value>> {
 public:
  using Grid =
      typename Internal::Traits<ScalarFieldExpansion<_Grid, _Value>>::Grid;
  using Value =
      typename Internal::Traits<ScalarFieldExpansion<_Grid, _Value>>::Value;
  using Int =
      typename Internal::Traits<ScalarFieldExpansion<_Grid, _Value>>::Int;
  using Real =
      typename Internal::Traits<ScalarFieldExpansion<_Grid, _Value>>::Real;
  using Complex =
      typename Internal::Traits<ScalarFieldExpansion<_Grid, _Value>>::Complex;
  using Scalar =
      typename Internal::Traits<ScalarFieldExpansion<_Grid, _Value>>::Scalar;
  using Writeable =
      typename Internal::Traits<ScalarFieldExpansion<_Grid, _Value>>::Writeable;

  // Return the grid.
  auto& GetGrid() const { return _grid; }

  // Read access to data.
  auto operator[](Int l, Int m) const { return _data[this->Index(l, m)]; }

  // Write access to data.
  auto& operator[](Int l, Int m)
  requires Writeable::value
  {
    return _data[this->Index(l, m)];
  }
  // Default constructor.
  ScalarFieldExpansion() = default;

  // Construct from grid initialising values to zero.
  ScalarFieldExpansion(_Grid& grid)
      : _grid{grid}, _data{FFTWpp::vector<Complex>(this->ExpansionSize())} {}

  // Construction from grid using function to initialise values.
  template <typename Function>
  requires ScalarFunctionS2Exapansion<Function, Int, Complex>
  ScalarFieldExpansion(_Grid& grid, Function&& f) : ScalarFieldExpansion(grid) {
    std::ranges::copy(
        this->Indices() | std::ranges::views::transform([&f](auto index) {
          auto [l, m] = index;
          return f(l, m);
        }),
        _data.begin());
  }

  // Construct from an element of the base class.
  template <typename Derived>
  requires std::convertible_to<typename Derived::Complex, Complex>
  ScalarFieldExpansion(const ScalarFieldExpansionBase<Derived>& other)
      : ScalarFieldExpansion(other.GetGrid()) {
    for (auto [l, m] : this->Indices()) {
      operator[](l, m) = other[l, m];
    }
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Complex, Complex>
  ScalarFieldExpansion(ScalarFieldExpansionBase<Derived>&& other)
      : ScalarFieldExpansion(other) {}

  // Default copy and move constructors.
  ScalarFieldExpansion(const ScalarFieldExpansion&) = default;
  ScalarFieldExpansion(ScalarFieldExpansion&&) = default;

  // Default copy and move assigment.
  ScalarFieldExpansion& operator=(const ScalarFieldExpansion&) = default;
  ScalarFieldExpansion& operator=(ScalarFieldExpansion&&) = default;

  // Use assignment defined in base class.
  using ScalarFieldExpansionBase<ScalarFieldExpansion<_Grid, _Value>>::operator=
      ;

  // Return view to the data.
  auto Data() { return std::ranges::views::all(_data); }

 private:
  _Grid& _grid;
  FFTWpp::vector<Complex> _data;
};

// Type aliases for real and complex fields.
template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealScalarFieldExpansion = ScalarFieldExpansion<Grid, RealValued>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexScalarFieldExpansion = ScalarFieldExpansion<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif