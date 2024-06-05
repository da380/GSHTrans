#ifndef GSH_TRANS_COMPLEXIFIED_VECTOR_FIELD_GUARD_H
#define GSH_TRANS_COMPLEXIFIED_VECTOR_FIELD_GUARD_H

#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "VectorFieldBase.h"
#include "VectorFieldConstComponent.h"

namespace GSHTrans {

// Forward declare class.
template <typename Derived>
requires requires() {
  requires std::same_as<typename Derived::Value, RealValued>;
  requires std::same_as<typename Derived::Grid::MRange, All> &&
               std::same_as<typename Derived::Grid::NRange, All>;
}
class ComplexifiedVectorField;

// Set traits.
namespace Internal {

template <typename Derived>
struct Traits<ComplexifiedVectorField<Derived>> {
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
class ComplexifiedVectorField
    : public VectorFieldBase<ComplexifiedVectorField<Derived>> {
 public:
  using Int = typename Internal::Traits<ComplexifiedVectorField<Derived>>::Int;
  using Grid =
      typename Internal::Traits<ComplexifiedVectorField<Derived>>::Grid;
  using Value =
      typename Internal::Traits<ComplexifiedVectorField<Derived>>::Value;
  using Real =
      typename Internal::Traits<ComplexifiedVectorField<Derived>>::Real;
  using Complex =
      typename Internal::Traits<ComplexifiedVectorField<Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<ComplexifiedVectorField<Derived>>::Scalar;
  using Writeable =
      typename Internal::Traits<ComplexifiedVectorField<Derived>>::Writeable;

  // Return grid.
  auto GetGrid() const { return _u.GetGrid(); }

  // Read access to data.
  auto operator[](Int alpha, Int iTheta, Int iPhi) const -> Complex {
    this->CheckPointIndices(iTheta, iPhi);
    this->CheckCanonicalIndices(alpha);
    constexpr auto i = Complex{0, 1};
    return alpha == 0 ? _u[0, iTheta, iPhi]
                      : alpha * _u[1, iTheta, iPhi] + i * _u[-1, iTheta, iPhi];
  }

  // Read access to component.
  auto operator[](Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return VectorFieldConstComponent(*this, alpha);
  }

  // Constructors.
  ComplexifiedVectorField(const VectorFieldBase<Derived>& u) : _u{u} {}

  ComplexifiedVectorField(const ComplexifiedVectorField&) = default;
  ComplexifiedVectorField(ComplexifiedVectorField&&) = default;

  // Assignment.
  ComplexifiedVectorField& operator=(const ComplexifiedVectorField&) = default;
  ComplexifiedVectorField& operator=(ComplexifiedVectorField&&) = default;

 private:
  const VectorFieldBase<Derived>& _u;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_COMPLEXIFIED_VECTOR_FIELD_GUARD_H