#ifndef GSH_TRANS_VECTOR_FIELD_BINAR_TO_SCALAR_FIELD_GUARD_H
#define GSH_TRANS_VECTOR_FIELD_BINAR_TO_SCALAR_FIELD_GUARD_H

#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "VectorFieldBase.h"
#include "VectorFieldConstComponent.h"

namespace GSHTrans {

// Forward declare class.
template <typename Derived0, typename Derived1, typename Function>
requires requires() {
  requires std::same_as<typename Derived0::Scalar, typename Derived1::Scalar>;
  requires std::invocable<Function, typename Derived0::Scalar,
                          typename Derived0::Scalar, typename Derived0::Scalar,
                          typename Derived0::Scalar, typename Derived0::Scalar,
                          typename Derived0::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived0::Scalar,
                           typename Derived0::Scalar, typename Derived0::Scalar,
                           typename Derived0::Scalar, typename Derived0::Scalar,
                           typename Derived0::Scalar>,
      typename Derived0::Scalar>;
}
class VectorFieldBinaryToScalarField;

// Set traits.
namespace Internal {

template <typename Derived0, typename Derived1, typename Function>
struct Traits<VectorFieldBinaryToScalarField<Derived0, Derived1, Function>> {
  using Int = typename Derived0::Int;
  using Grid = typename Derived0::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar = typename Derived0::Scalar;
  using Value = typename Derived0::Value;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename Derived0, typename Derived1, typename Function>
requires requires() {
  requires std::same_as<typename Derived0::Scalar, typename Derived1::Scalar>;
  requires std::invocable<Function, typename Derived0::Scalar,
                          typename Derived0::Scalar, typename Derived0::Scalar,
                          typename Derived0::Scalar, typename Derived0::Scalar,
                          typename Derived0::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived0::Scalar,
                           typename Derived0::Scalar, typename Derived0::Scalar,
                           typename Derived0::Scalar, typename Derived0::Scalar,
                           typename Derived0::Scalar>,
      typename Derived0::Scalar>;
}
class VectorFieldBinaryToScalarField
    : public ScalarFieldBase<
          VectorFieldBinaryToScalarField<Derived0, Derived1, Function>> {
 public:
  using Int = typename Internal::Traits<
      VectorFieldBinaryToScalarField<Derived0, Derived1, Function>>::Int;
  using Grid = typename Internal::Traits<
      VectorFieldBinaryToScalarField<Derived0, Derived1, Function>>::Grid;
  using Value = typename Internal::Traits<
      VectorFieldBinaryToScalarField<Derived0, Derived1, Function>>::Value;
  using Real = typename Internal::Traits<
      VectorFieldBinaryToScalarField<Derived0, Derived1, Function>>::Real;
  using Complex = typename Internal::Traits<
      VectorFieldBinaryToScalarField<Derived0, Derived1, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      VectorFieldBinaryToScalarField<Derived0, Derived1, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      VectorFieldBinaryToScalarField<Derived0, Derived1, Function>>::Writeable;

  // Methods needed to inherit from VectorField Base.
  auto GetGrid() const { return _u0.GetGrid(); }

  // Read access to data.
  auto operator[](Int iTheta, Int iPhi) const -> Scalar {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u0[-1, iTheta, iPhi], _u0[0, iTheta, iPhi], _u0[1, iTheta, iPhi],
              _u1[-1, iTheta, iPhi], _u1[0, iTheta, iPhi],
              _u1[1, iTheta, iPhi]);
    return 0;
  }

  // Constructors.
  VectorFieldBinaryToScalarField() = delete;
  VectorFieldBinaryToScalarField(const VectorFieldBase<Derived0>& u0,
                                 const VectorFieldBase<Derived1>& u1,
                                 Function&& f)
      : _u0{u0}, _u1{u1}, _f{f} {}

  VectorFieldBinaryToScalarField(const VectorFieldBinaryToScalarField&) =
      default;
  VectorFieldBinaryToScalarField(VectorFieldBinaryToScalarField&&) = default;

  // Assignment.
  VectorFieldBinaryToScalarField& operator=(VectorFieldBinaryToScalarField&) =
      default;
  VectorFieldBinaryToScalarField& operator=(VectorFieldBinaryToScalarField&&) =
      default;

 private:
  const VectorFieldBase<Derived0>& _u0;
  const VectorFieldBase<Derived1>& _u1;
  Function& _f;
};

}  // namespace GSHTrans

#endif