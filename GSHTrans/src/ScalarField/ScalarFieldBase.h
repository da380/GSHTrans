#ifndef GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"

namespace GSHTrans {

//-------------------------------------------------//
//               Define the base class             //
//-------------------------------------------------//
template <typename _Derived>
class ScalarFieldBase : public FieldBase<ScalarFieldBase<_Derived>> {
 public:
  using Int = typename FieldBase<ScalarFieldBase<_Derived>>::Int;

  // Methods related to the grid.
  auto GetGrid() const { return Derived().GetGrid(); }

  // Methods related to the data.
  auto Size() const { return GetGrid().FieldSize(); }
  auto operator()(Int iTheta, Int iPhi) const {
    return Derived().operator()(iTheta, iPhi);
  }

  void Print() const {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      std::cout << iTheta << " " << iPhi << " " << operator()(iTheta, iPhi)
                << std::endl;
    }
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_BASE_GUARD_H