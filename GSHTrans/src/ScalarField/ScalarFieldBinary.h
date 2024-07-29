#ifndef GSH_TRANS_SCALAR_FIELD_BINARY_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_BINARY_GUARD_H

#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "ScalarFieldBase.h"

namespace GSHTrans {

// Forward declare class.
template <typename Derived1, typename Derived2, typename Function>
requires requires() {
  requires std::same_as<typename Derived1::Real, typename Derived2::Real>;
  requires std::same_as<typename Derived1::Value, typename Derived2::Value>;
  requires std::invocable<Function, typename Derived1::Scalar,
                          typename Derived2::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived1::Scalar,
                           typename Derived2::Scalar>,
      typename Derived1::Scalar>;
}
class ScalarFieldBinary;

// Set traits.
namespace Internal {

template <typename Derived1, typename Derived2, typename Function>
struct Traits<ScalarFieldBinary<Derived1, Derived2, Function>> {
  using Int = std::ptrdiff_t;
  using Grid = typename Derived1::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Scalar = typename Derived1::Scalar;
  using Value = typename Derived1::Value;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename Derived1, typename Derived2, typename Function>
requires requires() {
  requires std::same_as<typename Derived1::Real, typename Derived2::Real>;
  requires std::same_as<typename Derived1::Value, typename Derived2::Value>;
  requires std::invocable<Function, typename Derived1::Scalar,
                          typename Derived2::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived1::Scalar,
                           typename Derived2::Scalar>,
      typename Derived1::Scalar>;
}
class ScalarFieldBinary
    : public ScalarFieldBase<ScalarFieldBinary<Derived1, Derived2, Function>> {
 public:
  using Int = typename Internal::Traits<
      ScalarFieldBinary<Derived1, Derived2, Function>>::Int;
  using Grid = typename Internal::Traits<
      ScalarFieldBinary<Derived1, Derived2, Function>>::Grid;
  using Value = typename Internal::Traits<
      ScalarFieldBinary<Derived1, Derived2, Function>>::Value;
  using Real = typename Internal::Traits<
      ScalarFieldBinary<Derived1, Derived2, Function>>::Real;
  using Complex = typename Internal::Traits<
      ScalarFieldBinary<Derived1, Derived2, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      ScalarFieldBinary<Derived1, Derived2, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      ScalarFieldBinary<Derived1, Derived2, Function>>::Writeable;

  // Methods needed to inherit from ScalarField Base.
  auto& GetGrid() const { return _u1.GetGrid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u1[iTheta, iPhi], _u2[iTheta, iPhi]);
  }

  // Constructors.
  ScalarFieldBinary() = delete;
  ScalarFieldBinary(const ScalarFieldBase<Derived1>& u1,
                    const ScalarFieldBase<Derived2>& u2, Function&& f)
      : _u1{u1}, _u2{u2}, _f{f} {
    assert(_u1.FieldSize() == _u2.FieldSize());
  }

  ScalarFieldBinary(const ScalarFieldBinary&) = default;
  ScalarFieldBinary(ScalarFieldBinary&&) = default;

  // Assignment.
  ScalarFieldBinary& operator=(ScalarFieldBinary&) = default;
  ScalarFieldBinary& operator=(ScalarFieldBinary&&) = default;

 private:
  const ScalarFieldBase<Derived1>& _u1;
  const ScalarFieldBase<Derived2>& _u2;
  Function& _f;
};

}  // namespace GSHTrans

#endif