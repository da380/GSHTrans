#ifndef GSH_TRANS_COMPLEXIFIED_SCALAR_FIELD_EXPRESSIONS_GUARD_H
#define GSH_TRANS_COMPLEXIFIED_SCALAR_FIELD_EXPRESSIONS_GUARD_H

#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "ScalarFieldBase.h"

namespace GSHTrans {

// Forward declare class.
template <typename Derived>
requires requires() {
  requires std::same_as<typename Derived::Value, RealValued>;
  requires std::same_as<typename Derived::Grid::MRange, All> &&
               std::same_as<typename Derived::Grid::NRange, All>;
}
class ComplexifiedScalarField;

// Set traits.
namespace Internal {

template <typename Derived>
struct Traits<ComplexifiedScalarField<Derived>> {
  using Int = std::ptrdiff_t;
  using Grid = typename Derived::Grid;
  using Value = ComplexValued;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar = Complex;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename Derived>
requires requires() {
  requires std::same_as<typename Derived::Value, RealValued>;
  requires std::same_as<typename Derived::Grid::MRange, All> &&
               std::same_as<typename Derived::Grid::NRange, All>;
}
class ComplexifiedScalarField
    : public ScalarFieldBase<ComplexifiedScalarField<Derived>> {
 public:
  using Int = typename Internal::Traits<ComplexifiedScalarField<Derived>>::Int;
  using Grid =
      typename Internal::Traits<ComplexifiedScalarField<Derived>>::Grid;
  using Value =
      typename Internal::Traits<ComplexifiedScalarField<Derived>>::Value;
  using Real =
      typename Internal::Traits<ComplexifiedScalarField<Derived>>::Real;
  using Complex =
      typename Internal::Traits<ComplexifiedScalarField<Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<ComplexifiedScalarField<Derived>>::Scalar;
  using Writeable =
      typename Internal::Traits<ComplexifiedScalarField<Derived>>::Writeable;

  // Methods needed to inherit from ScalarField Base.
  auto GetGrid() const { return _u.GetGrid(); }
  auto operator()(Int iTheta, Int iPhi) const -> Complex {
    this->CheckPointIndices(iTheta, iPhi);
    return _u(iTheta, iPhi);
  }

  // Constructors.
  ComplexifiedScalarField(const ScalarFieldBase<Derived>& u) : _u{u} {}

  ComplexifiedScalarField(const ComplexifiedScalarField&) = default;
  ComplexifiedScalarField(ComplexifiedScalarField&&) = default;

  // Assignment.
  ComplexifiedScalarField& operator=(const ComplexifiedScalarField&) = default;
  ComplexifiedScalarField& operator=(ComplexifiedScalarField&&) = default;

 private:
  const ScalarFieldBase<Derived>& _u;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_COMPLEXIFIED_SCALAR_FIELD_EXPRESSIONS_GUARD_H