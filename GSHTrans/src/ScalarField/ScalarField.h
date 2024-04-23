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

//-------------------------------------------------//
//       Scalar field that stores its data         //
//-------------------------------------------------//
template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class ScalarField : public ScalarFieldBase<ScalarField<_Grid, _Value>> {
 public:
  using Int = typename ScalarFieldBase<ScalarField<_Grid, _Value>>::Int;
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

  // Constructors.
  ScalarField() = default;

  ScalarField(_Grid grid)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->Size())} {}

  ScalarField(_Grid grid, Scalar s)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->Size(), s)} {}

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  ScalarField(const ScalarFieldBase<Derived>& other)
      : ScalarField(other.GetGrid()) {
    assert(this->Size() == other.Size());
    CopyValues(other);
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  ScalarField(ScalarFieldBase<Derived>&& other) : ScalarField(other) {}

  ScalarField(const ScalarField&) = default;
  ScalarField(ScalarField&&) = default;

  // Assignment.
  ScalarField& operator=(const ScalarField&) = default;
  ScalarField& operator=(ScalarField&&) = default;

  auto& operator=(Scalar s) {
    std::ranges::copy(std::ranges::views::repeat(s, this->Size()),
                      _data.begin());
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator=(const ScalarFieldBase<Derived>& other) {
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  auto& operator=(ScalarFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Iterators.
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  // Value assignement.
  auto& operator()(Int iTheta, Int iPhi) {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

  template <typename Function>
  requires requires() {
    requires std::invocable<Function, Real, Real>;
    requires std::convertible_to<std::invoke_result_t<Function, Real, Real>,
                                 Scalar>;
  }
  void Interpolate(Function f) {
    std::ranges::copy(_grid.InterpolateFunction(f), _data.begin());
  }

 private:
  _Grid _grid;
  FFTWpp::vector<Scalar> _data;

  template <typename Derived>
  void CopyValues(const ScalarFieldBase<Derived>& other) {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator()(iTheta, iPhi) = other(iTheta, iPhi);
    }
  }

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

#endif  // namespace GSHTrans
