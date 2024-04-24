#ifndef GSH_TRANS_ISOTROPIC_MATRIX_FIELD_BASE_GUARD_H
#define GSH_TRANS_ISOTROPIC_MATRIX_FIELD_BASE_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <vector>

#include "../Concepts.h"
#include "../GridBase.h"
#include "../MatrixField/MatrixFieldBase.h"
#include "../ScalarField/ScalarFieldBase.h"

namespace GSHTrans {

template <typename _Derived>
class IsotropicMatrixFieldBase
    : public MatrixFieldBase<IsotropicMatrixFieldBase<_Derived>> {
 public:
  // Methods related to the grid.
  auto GetGrid() const { return Derived().GetGrid(); }

  // Methods related to the data.
  auto operator()(Int alpha, Int beta, Int iTheta, Int iPhi) const {
    return Derived().operator()(alpha, beta, iTheta, iPhi);
  }
  auto operator()(Int alpha, Int beta) const {
    return Derived().operator()(alpha, beta);
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_ISOTROPIC_MATRIX_FIELD_BASE_GUARD_H