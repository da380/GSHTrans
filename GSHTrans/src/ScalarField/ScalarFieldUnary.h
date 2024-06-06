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
template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
               std::invoke_result_t<Function, typename Derived::Scalar>,
               typename Derived::Scalar> ||
               (std::convertible_to<
                    std::invoke_result_t<Function, typename Derived::Scalar>,
                    typename Derived::Complex> &&
                std::same_as<typename Derived::Grid::MRange, All> &&
                std::same_as<typename Derived::Grid::MRange, All>);
}
class ScalarFieldUnary;

// Set traits.
namespace Internal {

template <typename Derived, typename Function>
struct Traits<ScalarFieldUnary<Derived, Function>> {
  using Int = std::ptrdiff_t;
  using Grid = typename Derived::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar = std::invoke_result_t<Function, typename Derived::Scalar>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
               std::invoke_result_t<Function, typename Derived::Scalar>,
               typename Derived::Scalar> ||
               (std::convertible_to<
                    std::invoke_result_t<Function, typename Derived::Scalar>,
                    typename Derived::Complex> &&
                std::same_as<typename Derived::Grid::MRange, All> &&
                std::same_as<typename Derived::Grid::MRange, All>);
}
class ScalarFieldUnary
    : public ScalarFieldBase<ScalarFieldUnary<Derived, Function>> {
 public:
  using Int =
      typename Internal::Traits<ScalarFieldUnary<Derived, Function>>::Int;
  using Grid =
      typename Internal::Traits<ScalarFieldUnary<Derived, Function>>::Grid;
  using Value =
      typename Internal::Traits<ScalarFieldUnary<Derived, Function>>::Value;
  using Real =
      typename Internal::Traits<ScalarFieldUnary<Derived, Function>>::Real;
  using Complex =
      typename Internal::Traits<ScalarFieldUnary<Derived, Function>>::Complex;
  using Scalar =
      typename Internal::Traits<ScalarFieldUnary<Derived, Function>>::Scalar;
  using Writeable =
      typename Internal::Traits<ScalarFieldUnary<Derived, Function>>::Writeable;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u[iTheta, iPhi]);
  }

  // Constructors.
  ScalarFieldUnary() = delete;
  ScalarFieldUnary(const ScalarFieldBase<Derived>& u, Function&& f)
      : _u{u}, _f{f} {}

  ScalarFieldUnary(const ScalarFieldUnary&) = default;
  ScalarFieldUnary(ScalarFieldUnary&&) = default;

  // Assignment.
  ScalarFieldUnary& operator=(ScalarFieldUnary&) = default;
  ScalarFieldUnary& operator=(ScalarFieldUnary&&) = default;

 private:
  const ScalarFieldBase<Derived>& _u;
  Function& _f;
};

}  // namespace GSHTrans

#endif