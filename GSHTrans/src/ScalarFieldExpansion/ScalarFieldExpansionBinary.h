#ifndef GSH_TRANS_SCALAR_FIELD_EXPANSION_BINARY_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_EXPANSION_BINARY_GUARD_H

#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "ScalarFieldExpansionBase.h"

namespace GSHTrans {

// Forward declare class.
template <typename _Derived0, typename _Derived1, typename Function>
requires requires() {
  requires std::same_as<typename _Derived0::Value, typename _Derived1::Value>;
  requires std::same_as<typename _Derived0::Complex,
                        typename _Derived1::Complex>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename _Derived0::Complex,
                           typename _Derived1::Complex>,
      typename _Derived0::Complex>;
}
class ScalarFieldExpansionBinary;

// Set traits.
namespace Internal {

template <typename _Derived0, typename _Derived1, typename Function>
struct Traits<ScalarFieldExpansionBinary<_Derived0, _Derived1, Function>> {
  using Int = std::ptrdiff_t;
  using Real = typename _Derived0::Real;
  using Complex = typename _Derived0::Complex;
  using Value = typename _Derived0::Value;
  using Scalar = typename _Derived0::Scalar;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename _Derived0, typename _Derived1, typename Function>
requires requires() {
  requires std::same_as<typename _Derived0::Value, typename _Derived1::Value>;
  requires std::same_as<typename _Derived0::Complex,
                        typename _Derived1::Complex>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename _Derived0::Complex,
                           typename _Derived1::Complex>,
      typename _Derived0::Complex>;
}
class ScalarFieldExpansionBinary
    : public ScalarFieldExpansionBase<
          ScalarFieldExpansionBinary<_Derived0, _Derived1, Function>> {
 public:
  using Int = typename Internal::Traits<
      ScalarFieldExpansionBinary<_Derived0, _Derived1, Function>>::Int;
  using Value = typename Internal::Traits<
      ScalarFieldExpansionBinary<_Derived0, _Derived1, Function>>::Value;
  using Real = typename Internal::Traits<
      ScalarFieldExpansionBinary<_Derived0, _Derived1, Function>>::Real;
  using Complex = typename Internal::Traits<
      ScalarFieldExpansionBinary<_Derived0, _Derived1, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      ScalarFieldExpansionBinary<_Derived0, _Derived1, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      ScalarFieldExpansionBinary<_Derived0, _Derived1, Function>>::Writeable;

  // Methods needed to inherit from ScalarFieldExpansion Base.
  auto& Grid() const { return _u0.Grid(); }
  auto operator[](Int l, Int m) const { return _f(_u0[l, m], _u1[l, m]); }

  // Constructors.
  ScalarFieldExpansionBinary() = delete;
  ScalarFieldExpansionBinary(const ScalarFieldExpansionBase<_Derived0>& u0,
                             const ScalarFieldExpansionBase<_Derived0>& u1,
                             Function&& f)
      : _u0{u0}, _u1{u1}, _f{f} {
    assert(_u1.MaxDegree() == _u1.MaxDegree());
  }

  ScalarFieldExpansionBinary(const ScalarFieldExpansionBinary&) = default;
  ScalarFieldExpansionBinary(ScalarFieldExpansionBinary&&) = default;

  // Assignment.
  ScalarFieldExpansionBinary& operator=(ScalarFieldExpansionBinary&) = default;
  ScalarFieldExpansionBinary& operator=(ScalarFieldExpansionBinary&&) = default;

 private:
  const ScalarFieldExpansionBase<_Derived0>& _u0;
  const ScalarFieldExpansionBase<_Derived1>& _u1;
  Function& _f;
};

}  // namespace GSHTrans

#endif