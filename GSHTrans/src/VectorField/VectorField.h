#ifndef GSH_TRANS_VECTOR_FIELD_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../ScalarField/ScalarFieldBase.h"
#include "VectorFieldBase.h"
#include "VectorFieldComponent.h"
#include "VectorFieldConstComponent.h"

namespace GSHTrans {

// Forward declare class.
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class VectorField;

// Set traits.
namespace Internal {

template <typename _Grid, RealOrComplexValued _Value>
struct Traits<VectorField<_Grid, _Value>> {
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
class VectorField : public VectorFieldBase<VectorField<_Grid, _Value>> {
 public:
  using Grid = typename Internal::Traits<VectorField<_Grid, _Value>>::Grid;
  using Value = typename Internal::Traits<VectorField<_Grid, _Value>>::Value;
  using Int = typename Internal::Traits<VectorField<_Grid, _Value>>::Int;
  using Real = typename Internal::Traits<VectorField<_Grid, _Value>>::Real;
  using Complex =
      typename Internal::Traits<VectorField<_Grid, _Value>>::Complex;
  using Scalar = typename Internal::Traits<VectorField<_Grid, _Value>>::Scalar;
  using Writeable =
      typename Internal::Traits<VectorField<_Grid, _Value>>::Writeable;

  // Return the grid.
  auto GetGrid() const { return _grid; }

  // Read access to data.
  auto operator[](Int alpha, Int iTheta, Int iPhi) const {
    this->CheckCanonicalIndices(alpha);
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(alpha, iTheta, iPhi)];
  }

  // Write access to data.
  auto& operator[](Int alpha, Int iTheta, Int iPhi) {
    this->CheckCanonicalIndices(alpha);
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(alpha, iTheta, iPhi)];
  }

  // Return read-write access component
  auto operator[](Int alpha)
  requires Writeable::value
  {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponent(*this, alpha);
  }

  // Return read access component
  auto operator[](Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldConstComponent(*this, alpha);
  }

  // Default constructor.
  VectorField() = default;

  // Construct from grid initialising values to zero.
  VectorField(_Grid grid)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->FieldSize() * 3)} {}

  // Construction from grid initialising values with a function that
  // returns a GSHTrans::Vector
  template <typename Function>
  requires CanonicalVectorValuedFunction<Function, Real, Value>
  VectorField(_Grid grid, Function&& f) : VectorField(grid) {
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [point, index] :
           std::ranges::views::zip(this->Points(), this->PointIndices())) {
        auto [theta, phi] = point;
        auto [iTheta, iPhi] = index;
        operator[](alpha, iTheta, iPhi) = f(theta, phi)[alpha];
      }
    }
  }

  // Construct from an element of the base class.
  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  VectorField(const VectorFieldBase<Derived>& other)
      : VectorField(other.GetGrid()) {
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator[](alpha, iTheta, iPhi) = other[alpha, iTheta, iPhi];
      }
    }
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  VectorField(VectorFieldBase<Derived>&& other) : VectorField(other) {}

  // Default copy and move constructors.
  VectorField(const VectorField&) = default;
  VectorField(VectorField&&) = default;

  // Default copy and move assigment.
  VectorField& operator=(const VectorField&) = default;
  VectorField& operator=(VectorField&&) = default;

  // Return view to the data.
  auto Data() { return std::ranges::views::all(_data); }

 private:
  _Grid _grid;
  FFTWpp::vector<Scalar> _data;

  auto Offset(Int alpha) const { return (alpha + 1) * (this->FieldSize()); }

  auto Index(Int alpha, Int iTheta, int iPhi) const {
    return Offset(alpha) + iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

// Type aliases for real and complex fields.
template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealVectorField = VectorField<Grid, RealValued>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexVectorField = VectorField<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_VECTOR_FIELD_GUARD_H
