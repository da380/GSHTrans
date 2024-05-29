#ifndef GSH_TRANS_POINTWISE_UNARY_SCALAR_FIELD_EXPRESSIONS_GUARD_H
#define GSH_TRANS_POINTWISE_UNARY_SCALAR_FIELD_EXPRESSIONS_GUARD_H

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
class PointwiseUnaryScalarField;

// Set traits.
namespace Internal {

template <typename Derived, typename Function>
struct Traits<PointwiseUnaryScalarField<Derived, Function>> {
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
class PointwiseUnaryScalarField
    : public ScalarFieldBase<PointwiseUnaryScalarField<Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      PointwiseUnaryScalarField<Derived, Function>>::Int;
  using Grid = typename Internal::Traits<
      PointwiseUnaryScalarField<Derived, Function>>::Grid;
  using Value = typename Internal::Traits<
      PointwiseUnaryScalarField<Derived, Function>>::Value;
  using Real = typename Internal::Traits<
      PointwiseUnaryScalarField<Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      PointwiseUnaryScalarField<Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      PointwiseUnaryScalarField<Derived, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      PointwiseUnaryScalarField<Derived, Function>>::Writeable;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u(iTheta, iPhi));
  }

  // Constructors.
  PointwiseUnaryScalarField() = delete;
  PointwiseUnaryScalarField(const ScalarFieldBase<Derived>& u, Function f)
      : _u{u}, _f{f} {}

  PointwiseUnaryScalarField(const PointwiseUnaryScalarField&) = default;
  PointwiseUnaryScalarField(PointwiseUnaryScalarField&&) = default;

  // Assignment.
  PointwiseUnaryScalarField& operator=(PointwiseUnaryScalarField&) = default;
  PointwiseUnaryScalarField& operator=(PointwiseUnaryScalarField&&) = default;

 private:
  const ScalarFieldBase<Derived>& _u;
  Function _f;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_POINTWISE_UNARY_SCALAR_FIELD_EXPRESSIONS_GUARD_H