#ifndef GSH_TRANS_SCALAR_FIELD_BINARY_WITH_SCALAR_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_BINARY_WITH_SCALAR_GUARD_H

#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "ScalarFieldBase.h"

namespace GSHTrans {

// Forward declare class.
template <typename _Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename _Derived::Scalar,
                          typename _Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename _Derived::Scalar,
                           typename _Derived::Scalar>,
      typename _Derived::Scalar>;
}
class ScalarFieldBinaryWithScalar;

// Set traits.
namespace Internal {

template <typename _Derived, typename Function>
struct Traits<ScalarFieldBinaryWithScalar<_Derived, Function>> {
  using Int = std::ptrdiff_t;
  using Real = typename _Derived::Real;
  using Complex = typename _Derived::Complex;
  using Scalar = std::invoke_result_t<Function, typename _Derived::Scalar,
                                      typename _Derived::Scalar>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename _Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename _Derived::Scalar,
                          typename _Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename _Derived::Scalar,
                           typename _Derived::Scalar>,
      typename _Derived::Scalar>;
}
class ScalarFieldBinaryWithScalar
    : public ScalarFieldBase<ScalarFieldBinaryWithScalar<_Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      ScalarFieldBinaryWithScalar<_Derived, Function>>::Int;
  using Value = typename Internal::Traits<
      ScalarFieldBinaryWithScalar<_Derived, Function>>::Value;
  using Real = typename Internal::Traits<
      ScalarFieldBinaryWithScalar<_Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      ScalarFieldBinaryWithScalar<_Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      ScalarFieldBinaryWithScalar<_Derived, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      ScalarFieldBinaryWithScalar<_Derived, Function>>::Writeable;

  // Methods needed to inherit from ScalarField Base.
  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u[iTheta, iPhi], _s);
  }

  // Constructors.
  ScalarFieldBinaryWithScalar() = delete;
  ScalarFieldBinaryWithScalar(const ScalarFieldBase<_Derived>& u, Scalar s,
                              Function&& f)
      : _u{u}, _f{f}, _s{s} {}

  ScalarFieldBinaryWithScalar(const ScalarFieldBinaryWithScalar&) = default;
  ScalarFieldBinaryWithScalar(ScalarFieldBinaryWithScalar&&) = default;

  // Assignment.
  ScalarFieldBinaryWithScalar& operator=(ScalarFieldBinaryWithScalar&) =
      default;
  ScalarFieldBinaryWithScalar& operator=(ScalarFieldBinaryWithScalar&&) =
      default;

 private:
  const ScalarFieldBase<_Derived>& _u;
  Function& _f;
  Scalar _s;
};

}  // namespace GSHTrans

#endif