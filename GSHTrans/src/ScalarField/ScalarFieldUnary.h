#ifndef GSH_TRANS_SCALAR_FIELD_UNARY_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_UNARY_GUARD_H

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
  requires std::invocable<Function, typename _Derived::Scalar>;
  requires std::convertible_to<
               std::invoke_result_t<Function, typename _Derived::Scalar>,
               typename _Derived::Scalar> ||
               (std::convertible_to<
                    std::invoke_result_t<Function, typename _Derived::Scalar>,
                    typename _Derived::Complex> &&
                std::same_as<typename _Derived::Grid::MRange, All> &&
                std::same_as<typename _Derived::Grid::MRange, All>);
}
class ScalarFieldUnary;

// Set traits.
namespace Internal {

template <typename _Derived, typename Function>
struct Traits<ScalarFieldUnary<_Derived, Function>> {
  using Int = std::ptrdiff_t;
  using Real = typename _Derived::Real;
  using Complex = typename _Derived::Complex;
  using Scalar = std::invoke_result_t<Function, typename _Derived::Scalar>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename _Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename _Derived::Scalar>;
  requires std::convertible_to<
               std::invoke_result_t<Function, typename _Derived::Scalar>,
               typename _Derived::Scalar> ||
               (std::convertible_to<
                    std::invoke_result_t<Function, typename _Derived::Scalar>,
                    typename _Derived::Complex> &&
                std::same_as<typename _Derived::Grid::MRange, All> &&
                std::same_as<typename _Derived::Grid::MRange, All>);
}
class ScalarFieldUnary
    : public ScalarFieldBase<ScalarFieldUnary<_Derived, Function>> {
 public:
  using Int =
      typename Internal::Traits<ScalarFieldUnary<_Derived, Function>>::Int;
  using Value =
      typename Internal::Traits<ScalarFieldUnary<_Derived, Function>>::Value;
  using Real =
      typename Internal::Traits<ScalarFieldUnary<_Derived, Function>>::Real;
  using Complex =
      typename Internal::Traits<ScalarFieldUnary<_Derived, Function>>::Complex;
  using Scalar =
      typename Internal::Traits<ScalarFieldUnary<_Derived, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      ScalarFieldUnary<_Derived, Function>>::Writeable;

  // Methods needed to inherit from ScalarField Base.
  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u[iTheta, iPhi]);
  }

  // Constructors.
  ScalarFieldUnary() = delete;
  ScalarFieldUnary(const ScalarFieldBase<_Derived>& u, Function&& f)
      : _u{u}, _f{f} {}

  ScalarFieldUnary(const ScalarFieldUnary&) = default;
  ScalarFieldUnary(ScalarFieldUnary&&) = default;

  // Assignment.
  ScalarFieldUnary& operator=(ScalarFieldUnary&) = default;
  ScalarFieldUnary& operator=(ScalarFieldUnary&&) = default;

 private:
  const ScalarFieldBase<_Derived>& _u;
  Function& _f;
};

}  // namespace GSHTrans

#endif