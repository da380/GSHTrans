#ifndef GSH_TRANS_MATRIX_FIELD_BASE_GUARD_H
#define GSH_TRANS_MATRIX_FIELD_BASE_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "../ScalarField/ScalarFieldBase.h"

namespace GSHTrans {

//-------------------------------------------------//
//               Define the base class             //
//-------------------------------------------------//
template <typename _Derived>
class MatrixFieldBase : public FieldBase<MatrixFieldBase<_Derived>> {
 public:
  using Int = typename FieldBase<MatrixFieldBase<_Derived>>::Int;

  // Methods related to the grid.
  auto GetGrid() const { return Derived().GetGrid(); }

  // Methods related to the data.
  auto ComponentSize() const { return GetGrid().FieldSize(); }

  auto CanonicalIndices() const {
    return std::ranges::views::cartesian_product(
        std::ranges::views::iota(-1, 2), std::ranges::views::iota(-1, 2));
  }

  void CheckCanonicalIndices(Int alpha, Int beta) const {
    assert(std::abs(alpha) <= 1);
    assert(std::abs(beta) <= 1);
  }

  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    return Derived().operator()(alpha, beta, iTheta, iPhi);
  }
  auto operator()(Int alpha, Int beta) const {
    return Derived().operator()(alpha, beta);
  }

  void Print() const {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      for (auto alpha = -1; alpha <= 1; alpha++) {
        for (auto beta = -1; beta <= 1; beta++) {
          std::cout << operator()(alpha, beta, iTheta, iPhi) << " ";
        }
        std::cout << std::endl;
      }
      std::cout << "-------------" << std::endl;
    }
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

//-------------------------------------------------//
//              Matrix component view              //
//-------------------------------------------------//
template <typename Derived>
class MatrixFieldComponentView
    : public ScalarFieldBase<MatrixFieldComponentView<Derived>> {
 public:
  using Int = typename ScalarFieldBase<MatrixFieldComponentView<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _A.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    return _A(_alpha, _beta, iTheta, iPhi);
  }

  // Constructors.
  MatrixFieldComponentView(const MatrixFieldBase<Derived>& A, Int alpha,
                           Int beta)
      : _A{A}, _alpha{alpha}, _beta{beta} {}

  MatrixFieldComponentView(const MatrixFieldComponentView&) = default;
  MatrixFieldComponentView(MatrixFieldComponentView&&) = default;

  // Assignment.
  MatrixFieldComponentView& operator=(const MatrixFieldComponentView&) =
      default;
  MatrixFieldComponentView& operator=(MatrixFieldComponentView&&) = default;

 private:
  const VectorFieldBase<Derived>& _A;
  Int _alpha;
  Int _beta;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_MATRIX_FIELD_BASE_GUARD_H