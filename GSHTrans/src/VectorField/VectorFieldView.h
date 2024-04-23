#ifndef GSH_TRANS_VECTOR_FIELD_VIEW_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_VIEW_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "VectorFieldBase.h"

namespace GSHTrans {

//-------------------------------------------------//
//      Vector field with a view to its data       //
//-------------------------------------------------//
template <typename _Grid, std::ranges::view _View>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>;
  requires std::ranges::random_access_range<_View>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
}
class VectorFieldView : public VectorFieldBase<VectorFieldView<_Grid, _View>> {
 public:
  using Int = typename VectorFieldBase<VectorFieldView<_Grid, _View>>::Int;
  using Grid = _Grid;
  using Scalar = std::ranges::range_value_t<_View>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return _data[Index(alpha, iTheta, iPhi)];
  }
  auto operator()(Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldComponentView(*this, alpha);
  }

  // Constructors.
  VectorFieldView() = default;

  VectorFieldView(_Grid grid, _View data) : _grid{grid}, _data{data} {
    assert(3 * (this->ComponentSize()) == _data.size());
  }

  VectorFieldView(const VectorFieldView&) = default;
  VectorFieldView(VectorFieldView&&) = default;

  // Assignment.
  VectorFieldView& operator=(const VectorFieldView&) = default;
  VectorFieldView& operator=(VectorFieldView&&) = default;

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(const VectorFieldBase<Derived>& other) {
    assert(this->ComponentSize() == other.CompnentSize());
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(VectorFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Methods to make it a range.
  auto size() const { return this->Size(); }
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  // Value assignement.
  auto& operator()(Int alpha, Int iTheta, Int iPhi) {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    return _data[Index(alpha, iTheta, iPhi)];
  }
  auto operator()(Int alpha) {
    this->CheckCanonicalIndices(alpha);
    auto start = std::next(begin(), Offset(alpha));
    auto finish = std::next(start, this->ComponentSize());
    auto data = std::ranges::subrange(start, finish);
    return ScalarFieldView(_grid, data);
  }

 private:
  _Grid _grid;
  _View _data;

  template <typename Derived>
  void CopyValues(const VectorFieldBase<Derived>& other) {
    for (auto alpha : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator()(alpha, iTheta, iPhi) = other(alpha, iTheta, iPhi);
      }
    }
  }

  auto Offset(Int alpha) const { return (alpha + 1) * (this->ComponentSize()); }

  auto Index(Int alpha, Int iTheta, int iPhi) const {
    return Offset(alpha) + iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_VECTOR_FIELD_VIEW_GUARD_H