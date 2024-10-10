#ifndef GSH_TRANS_CANONICAL_COMPONENT_FIELD_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <cstddef>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "CanonicalComponentFieldBase.h"

namespace GSHTrans {

// Forward declare class.
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class CanonicalComponentField;

// Set traits.
namespace Internal {

template <typename _Grid, RealOrComplexValued _Value>
struct Traits<CanonicalComponentField<_Grid, _Value>> {
  using Int = std::ptrdiff_t;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<Value, RealValued>, Real, Complex>;
  using Writeable = std::true_type;
};

}  // namespace Internal

template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class CanonicalComponentField : public CanonicalComponentFieldBase<
                                    CanonicalComponentField<_Grid, _Value>> {
 public:
  using Value =
      typename Internal::Traits<CanonicalComponentField<_Grid, _Value>>::Value;
  using Int =
      typename Internal::Traits<CanonicalComponentField<_Grid, _Value>>::Int;
  using Real =
      typename Internal::Traits<CanonicalComponentField<_Grid, _Value>>::Real;
  using Complex = typename Internal::Traits<
      CanonicalComponentField<_Grid, _Value>>::Complex;
  using Scalar =
      typename Internal::Traits<CanonicalComponentField<_Grid, _Value>>::Scalar;
  using Writeable = typename Internal::Traits<
      CanonicalComponentField<_Grid, _Value>>::Writeable;

  // Return the upper index.
  auto& UpperIndex() const { return _N; }

  // Return the grid.
  auto& Grid() const { return _grid; }

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
  CanonicalComponentField() = default;

  // Construct from grid initialising values to zero.
  CanonicalComponentField(Int N, _Grid& grid)
      : _N{N}, _grid{grid}, _data{FFTWpp::vector<Scalar>(this->FieldSize())} {}

  // Construction from grid initialising values with a function.
  template <typename Function>
  requires ScalarFunctionS2<Function, Real, Scalar>
  CanonicalComponentField(Int N, _Grid& grid, Function&& f)
      : CanonicalComponentField(N, grid) {
    std::ranges::copy(_grid.ProjectFunction(f), _data.begin());
  }

  // Construct from an element of the base class.
  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  CanonicalComponentField(const CanonicalComponentFieldBase<Derived>& other)
      : CanonicalComponentField(other.UpperIndex(), other.Grid()) {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) = other[iTheta, iPhi];
    }
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  CanonicalComponentField(CanonicalComponentFieldBase<Derived>&& other)
      : CanonicalComponentField(other) {}

  // Default copy and move constructors.
  CanonicalComponentField(const CanonicalComponentField&) = default;
  CanonicalComponentField(CanonicalComponentField&&) = default;

  // Default copy and move assignment.
  CanonicalComponentField& operator=(const CanonicalComponentField&) = default;
  CanonicalComponentField& operator=(CanonicalComponentField&&) = default;

  // Use assignment defined in base class.
  using CanonicalComponentFieldBase<
      CanonicalComponentField<_Grid, _Value>>::operator=;

  // Return view to the data.
  auto Data() { return std::ranges::views::all(_data); }

 private:
  Int _N;
  _Grid& _grid;
  FFTWpp::vector<Scalar> _data;

  auto Index(Int iTheta, int iPhi) const {
    return iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

// Type aliases for real and complex fields.
template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealCanonicalComponentField = CanonicalComponentField<Grid, RealValued>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexCanonicalComponentField =
    CanonicalComponentField<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_GUARD_H
