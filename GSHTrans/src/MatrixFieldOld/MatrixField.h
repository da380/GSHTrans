#ifndef GSH_TRANS_MATRIX_FIELD_GUARD_H
#define GSH_TRANS_MATRIX_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "MatrixFieldBase.h"

namespace GSHTrans {

template <typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class MatrixField : public MatrixFieldBase<MatrixField<_Grid, _Value>> {
 public:
  using Int = typename MatrixFieldBase<MatrixField<_Grid, _Value>>::Int;
  using Grid = _Grid;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;
  using SelfAdjoint = std::false_type;
  using Diagonal = std::false_type;
  using Isotropic = std::false_type;

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
  MatrixField() = default;

  MatrixField(_Grid grid)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(9 * (this->FieldSize()))} {}

  MatrixField(_Grid grid, std::array<Scalar, 9>&& A) : MatrixField(grid) {
    for (auto [i, index] :
         std::ranges::views::enumerate(this->CanonicalIndices())) {
      auto [alpha, beta] = index;
      this->operator()(alpha, beta) = A[i];
    }
  }

  // Assignment.
  MatrixField& operator=(const MatrixField&) = default;
  MatrixField& operator=(MatrixField&&) = default;

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
  FFTWpp::vector<Scalar> _data;

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

// Type aliases for real and complex Matrix fields
template <typename Grid>
using RealMatrixField = MatrixField<Grid, RealValued>;

template <typename Grid>
using ComplexMatrixField = MatrixField<Grid, ComplexValued>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_MATRIX_FIELD_GUARD_H