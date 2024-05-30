#ifndef GSH_TRANS_SCALAR_FIELD_POINTWISE_UNARY_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_POINTWISE_UNARY_GUARD_H

#include <FFTWpp/Core>
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
class ScalarFieldPointwiseUnary;

// Set traits.
namespace Internal {

template <typename Derived, typename Function>
struct Traits<ScalarFieldPointwiseUnary<Derived, Function>> {
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
class ScalarFieldPointwiseUnary
    : public ScalarFieldBase<ScalarFieldPointwiseUnary<Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      ScalarFieldPointwiseUnary<Derived, Function>>::Int;
  using Grid = typename Internal::Traits<
      ScalarFieldPointwiseUnary<Derived, Function>>::Grid;
  using Value = typename Internal::Traits<
      ScalarFieldPointwiseUnary<Derived, Function>>::Value;
  using Real = typename Internal::Traits<
      ScalarFieldPointwiseUnary<Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      ScalarFieldPointwiseUnary<Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      ScalarFieldPointwiseUnary<Derived, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      ScalarFieldPointwiseUnary<Derived, Function>>::Writeable;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u(iTheta, iPhi));
  }

  // Constructors.
  ScalarFieldPointwiseUnary() = delete;
  ScalarFieldPointwiseUnary(const ScalarFieldBase<Derived>& u, Function f)
      : _u{u}, _f{f} {}

  ScalarFieldPointwiseUnary(const ScalarFieldPointwiseUnary&) = default;
  ScalarFieldPointwiseUnary(ScalarFieldPointwiseUnary&&) = default;

  // Assignment.
  ScalarFieldPointwiseUnary& operator=(ScalarFieldPointwiseUnary&) = default;
  ScalarFieldPointwiseUnary& operator=(ScalarFieldPointwiseUnary&&) = default;

 private:
  const ScalarFieldBase<Derived>& _u;
  Function _f;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SCALAR_FIELD_POINTWISE_UNARY_EXPRESSIONS_GUARD_H