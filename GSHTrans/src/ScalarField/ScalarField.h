#ifndef GSH_TRANS_SCALAR_FIELD_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "ScalarFieldBase.h"

namespace GSHTrans {

// Forward declare class.
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ScalarField;

// Set traits.
namespace Internal {

template <typename _Grid, RealOrComplexValued _Value>
struct Traits<ScalarField<_Grid, _Value>> {
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
class ScalarField : public ScalarFieldBase<ScalarField<_Grid, _Value>> {
 public:
  using Grid = typename Internal::Traits<ScalarField<_Grid, _Value>>::Grid;
  using Value = typename Internal::Traits<ScalarField<_Grid, _Value>>::Value;
  using Int = typename Internal::Traits<ScalarField<_Grid, _Value>>::Int;
  using Real = typename Internal::Traits<ScalarField<_Grid, _Value>>::Real;
  using Complex =
      typename Internal::Traits<ScalarField<_Grid, _Value>>::Complex;
  using Scalar = typename Internal::Traits<ScalarField<_Grid, _Value>>::Scalar;
  using Writeable =
      typename Internal::Traits<ScalarField<_Grid, _Value>>::Writeable;

  // Return the grid.
  auto& GetGrid() const { return _grid; }

  // Read access to data.
  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

  // Write access to data.
  auto& operator[](Int iTheta, Int iPhi) {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

  // Default constructor.
  ScalarField() = default;

  // Construct from grid initialising values to zero.
  ScalarField(_Grid& grid)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->FieldSize())} {}

  // Construction from grid initialising values with a function.
  template <typename Function>
  requires ScalarFunctionS2<Function, Real, Scalar>
  ScalarField(_Grid& grid, Function&& f) : ScalarField(grid) {
    std::ranges::copy(_grid.ProjectFunction(f), _data.begin());
  }

  // Construct from an element of the base class.
  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  ScalarField(const ScalarFieldBase<Derived>& other)
      : ScalarField(other.GetGrid()) {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) = other[iTheta, iPhi];
    }
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  ScalarField(ScalarFieldBase<Derived>&& other) : ScalarField(other) {}

  // Default copy and move constructors.
  ScalarField(const ScalarField&) = default;
  ScalarField(ScalarField&&) = default;

  // Default copy and move assigment.
  ScalarField& operator=(const ScalarField&) = default;
  ScalarField& operator=(ScalarField&&) = default;

  // Use assignment defined in base class.
  using ScalarFieldBase<ScalarField<_Grid, _Value>>::operator=;

  // Return view to the data.
  auto Data() { return std::ranges::views::all(_data); }

 private:
  _Grid& _grid;
  FFTWpp::vector<Scalar> _data;

  auto Index(Int iTheta, int iPhi) const {
    return iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

// Type aliases for real and complex fields.
template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealScalarField = ScalarField<Grid, RealValued>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexScalarField = ScalarField<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_GUARD_H
