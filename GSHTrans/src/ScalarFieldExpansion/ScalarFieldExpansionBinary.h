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
template <typename Derived0, typename Derived1, typename Function>
requires requires() {
  requires std::same_as<typename Derived0::Value, typename Derived1::Value>;
  requires std::same_as<typename Derived0::Complex, typename Derived1::Complex>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived0::Complex,
                           typename Derived1::Complex>,
      typename Derived0::Complex>;
}
class ScalarFieldExpansionBinary;

// Set traits.
namespace Internal {

template <typename Derived0, typename Derived1, typename Function>
struct Traits<ScalarFieldExpansionBinary<Derived0, Derived1, Function>> {
  using Int = std::ptrdiff_t;
  using Grid = typename Derived0::Grid;
  using Real = typename Grid::Real;
  using Complex = typename Grid::Complex;
  using Value = typename Derived0::Value;
  using Scalar = typename Derived0::Scalar;
  using Writeable = std::false_type;
};

}  // namespace Internal

template <typename Derived0, typename Derived1, typename Function>
requires requires() {
  requires std::same_as<typename Derived0::Value, typename Derived1::Value>;
  requires std::same_as<typename Derived0::Complex, typename Derived1::Complex>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived0::Complex,
                           typename Derived1::Complex>,
      typename Derived0::Complex>;
}
class ScalarFieldExpansionBinary
    : public ScalarFieldExpansionBase<
          ScalarFieldExpansionBinary<Derived0, Derived1, Function>> {
 public:
  using Int = typename Internal::Traits<
      ScalarFieldExpansionBinary<Derived0, Derived1, Function>>::Int;
  using Grid = typename Internal::Traits<
      ScalarFieldExpansionBinary<Derived0, Derived1, Function>>::Grid;
  using Value = typename Internal::Traits<
      ScalarFieldExpansionBinary<Derived0, Derived1, Function>>::Value;
  using Real = typename Internal::Traits<
      ScalarFieldExpansionBinary<Derived0, Derived1, Function>>::Real;
  using Complex = typename Internal::Traits<
      ScalarFieldExpansionBinary<Derived0, Derived1, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      ScalarFieldExpansionBinary<Derived0, Derived1, Function>>::Scalar;
  using Writeable = typename Internal::Traits<
      ScalarFieldExpansionBinary<Derived0, Derived1, Function>>::Writeable;

  // Methods needed to inherit from ScalarFieldExpansion Base.
  auto& GetGrid() const { return _u0.GetGrid(); }
  auto operator[](Int l, Int m) const { return _f(_u0[l, m], _u1[l, m]); }

  // Constructors.
  ScalarFieldExpansionBinary() = delete;
  ScalarFieldExpansionBinary(const ScalarFieldExpansionBase<Derived0>& u0,
                             const ScalarFieldExpansionBase<Derived0>& u1,
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
  const ScalarFieldExpansionBase<Derived0>& _u0;
  const ScalarFieldExpansionBase<Derived1>& _u1;
  Function& _f;
};

}  // namespace GSHTrans

#endif