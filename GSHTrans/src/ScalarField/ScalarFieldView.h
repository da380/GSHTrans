#ifndef GSH_TRANS_SCALAR_FIELD_VIEW_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_VIEW_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "ScalarField.h"

namespace GSHTrans {

//-------------------------------------------------//
//      Scalar field with a view to its data       //
//-------------------------------------------------//
template <typename _Grid, std::ranges::view _View>
requires requires() {
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>;
  requires std::ranges::random_access_range<_View>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
}
class ScalarFieldView : public ScalarFieldBase<ScalarFieldView<_Grid, _View>> {
 public:
  using Int = typename ScalarFieldBase<ScalarFieldView<_Grid, _View>>::Int;
  using Grid = _Grid;
  using Scalar = std::ranges::range_value_t<_View>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

  // Constructors.
  ScalarFieldView() = default;

  ScalarFieldView(_Grid grid, _View data) : _grid{grid}, _data{data} {
    assert(this->FieldSize() == _data.size());
  }

  ScalarFieldView(const ScalarFieldView&) = default;
  ScalarFieldView(ScalarFieldView&&) = default;

  // Assignment.
  ScalarFieldView& operator=(const ScalarFieldView&) = default;
  ScalarFieldView& operator=(ScalarFieldView&&) = default;

  auto& operator=(Scalar s) {
    std::ranges::copy(std::ranges::views::repeat(s, this->FieldSize()),
                      _data.begin());
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(const ScalarFieldBase<Derived>& other) {
    assert(this->FieldSize() == other.FieldSize());
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(ScalarFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

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
  _View _data;

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

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_VIEW_GUARD_H