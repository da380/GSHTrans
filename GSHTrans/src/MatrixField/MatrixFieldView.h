#ifndef GSH_TRANS_MATRIX_FIELD_VIEW_GUARD_H
#define GSH_TRANS_MATRIX_FIELD_VIEW_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "MatrixFieldBase.h"

namespace GSHTrans {

//-------------------------------------------------//
//       Matrix field with a view to its data      //
//-------------------------------------------------//
template <typename _Grid, std::ranges::view _View>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::ranges::output_range<_View, std::ranges::range_value_t<_View>>;
  requires std::ranges::random_access_range<_View>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
}
class MatrixFieldView : public MatrixFieldBase<MatrixFieldView<_Grid, _View>> {
 public:
  using Int = typename MatrixFieldBase<MatrixFieldView<_Grid, _View>>::Int;
  using Grid = _Grid;
  using Scalar = std::ranges::range_value_t<_View>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;

  // Methods needed to inherit from MatrixField Base.
  auto GetGrid() const { return _grid; }
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    this->CheckCanonicalIndices(alpha, beta);
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(alpha, beta, iTheta, iPhi)];
  }
  auto operator()(Int alpha, Int beta) const {
    this->CheckCanonicalIndices(alpha, beta);
    return MatrixFieldComponentView(*this, alpha, beta);
  }

  // Constructors.
  MatrixFieldView() = default;

  MatrixFieldView(_Grid grid, _View data) : _grid{grid}, _data{data} {
    assert(9 * (this->FieldSize()) == _data.size());
  }

  // Assignment.
  MatrixFieldView& operator=(const MatrixFieldView&) = default;
  MatrixFieldView& operator=(MatrixFieldView&&) = default;

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(const MatrixFieldBase<Derived>& other) {
    assert(this->FieldSize() == other.FieldSize());
    CopyValues(other);
    return *this;
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::Value, Value>
  auto& operator=(MatrixFieldBase<Derived>&& other) {
    *this = other;
    return *this;
  }

  // Value assignment.
  auto& operator()(Int alpha, Int beta, Int iTheta, Int iPhi) {
    this->CheckCanonicalIndices(alpha, beta);
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(alpha, beta, iTheta, iPhi)];
  }
  auto operator()(Int alpha, Int beta) {
    this->CheckCanonicalIndices(alpha, beta);
    auto start = std::next(_data.begin(), Offset(alpha, beta));
    auto finish = std::next(start, this->FieldSize());
    auto data = std::ranges::subrange(start, finish);
    return ScalarFieldView(_grid, data);
  }

 private:
  _Grid _grid;
  _View _data;

  auto Offset(Int alpha, Int beta) const {
    return (3 * (alpha + 1) + (beta + 1)) * (this->FieldSize());
  }

  template <typename Derived>
  void CopyValues(const MatrixFieldBase<Derived>& other) {
    for (auto [alpha, beta] : this->CanonicalIndices()) {
      for (auto [iTheta, iPhi] : this->PointIndices()) {
        operator()(alpha, beta, iTheta, iPhi) =
            other(alpha, beta, iTheta, iPhi);
      }
    }
  }

  auto Index(Int alpha, Int beta, Int iTheta, int iPhi) const {
    return Offset(alpha, beta) + iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_MATRIX_FIELD_VIEW_GUARD_H