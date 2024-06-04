#ifndef GSH_TRANS_VECTOR_FIELD_BASE_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_BASE_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../ScalarField/ScalarFieldBase.h"

namespace GSHTrans {

//-------------------------------------------------//
//               Define the base class             //
//-------------------------------------------------//
template <typename _Derived>
class VectorFieldBase : public FieldBase<VectorFieldBase<_Derived>> {
 public:
  using Int = typename FieldBase<VectorFieldBase<_Derived>>::Int;

  // Methods related to the grid.
  auto GetGrid() const { return Derived().GetGrid(); }

  // Methods related to the data.
  auto CanonicalIndices() const { return std::ranges::views::iota(-1, 2); }

  void CheckCanonicalIndices(Int alpha) const { assert(std::abs(alpha) <= 1); }

  auto operator()(Int alpha, Int iTheta, Int iPhi) const {
    return Derived().operator()(alpha, iTheta, iPhi);
  }
  auto operator()(Int alpha) const { return Derived().operator()(alpha); }

  void Print() const {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      for (auto alpha : CanonicalIndices()) {
        std::cout << operator()(alpha, iTheta, iPhi) << " ";
      }
      std::cout << std::endl;
    }
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

//-------------------------------------------------//
//              Vector component view              //
//-------------------------------------------------//
template <typename Derived>
class VectorFieldComponentView
    : public ScalarFieldBase<VectorFieldComponentView<Derived>> {
 public:
  using Int = typename ScalarFieldBase<VectorFieldComponentView<Derived>>::Int;
  using Grid = typename Derived::Grid;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _u(_alpha, iTheta, iPhi);
  }

  // Constructors.
  VectorFieldComponentView(VectorFieldBase<Derived>& u, Int alpha)
      : _u{u}, _alpha{alpha} {
    this->CheckCanonicalIndices(alpha);
  }

  VectorFieldComponentView(const VectorFieldComponentView&) = default;
  VectorFieldComponentView(VectorFieldComponentView&&) = default;

  // Assignment.
  VectorFieldComponentView& operator=(const VectorFieldComponentView&) =
      default;
  VectorFieldComponentView& operator=(VectorFieldComponentView&&) = default;

 private:
  VectorFieldBase<Derived>& _u;
  Int _alpha;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_VECTOR_FIELD_BASE_GUARD_H