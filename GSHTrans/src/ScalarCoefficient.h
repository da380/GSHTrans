#ifndef GSH_TRANS_SCALAR_COEFFICIENT_GUARD_H
#define GSH_TRANS_SCALAR_COEFFICIENT_GUARD_H

#include <algorithm>
#include <cassert>
#include <cmath>

namespace GHSTrans {

template <typename _Derived>
class ScalarCoefficientBase
    : public CoefficientBase<ScalarCoefficientBase<_Derived>> {
  using Int = typename ScalarCoefficientBase<_Derived>::Int;
};

}  // namespace GHSTrans

#endif  // GSH_TRANS_SCALAR_COEFFICIENT_GUARD_H