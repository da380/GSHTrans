#ifndef GSH_TRANS_CANONICAL_COMPONENT_FIELD_BINARY_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_FIELD_BINARY_GUARD_H

#include <complex>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../Utility.h"
#include "CanonicalComponentFieldBase.h"

namespace GSHTrans {

// Forward declare the classes.

template <typename Derived0, typename Derived1, typename Function>
requires requires() {
  requires std::same_as<typename Derived0::Value, typename Derived1::Value>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived0::Scalar,
                           typename Derived1::Scalar>,
      typename Derived0::Scalar>;
}
class CanonicalComponentFieldBinary;

// Set up the traits

namespace Internal {

template <typename Derived0, typename Derived1, typename Function>
struct Traits<CanonicalComponentFieldBinary<Derived0, Derived1, Function>> {
  using Int = typename Derived0::Int;
  using Real = typename Derived0::Real;
  using Complex = typename Derived0::Complex;
  using Scalar = typename Derived0::Scalar;
  using Value = typename Derived0::Value;
  using Writeable = std::false_type;
};
}  // namespace Internal

template <typename Derived0, typename Derived1, typename Function>
requires requires() {
  requires std::same_as<typename Derived0::Value, typename Derived1::Value>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived0::Scalar,
                           typename Derived1::Scalar>,
      typename Derived0::Scalar>;
}
class CanonicalComponentFieldBinary
    : public CanonicalComponentFieldBase<
          CanonicalComponentFieldBinary<Derived0, Derived1, Function>> {
 public:
  using Int = Internal::Traits<
      CanonicalComponentFieldBinary<Derived0, Derived1, Function>>::Int;
  using Real = Internal::Traits<
      CanonicalComponentFieldBinary<Derived0, Derived1, Function>>::Real;
  using Complex = Internal::Traits<
      CanonicalComponentFieldBinary<Derived0, Derived1, Function>>::Complex;
  using Scalar = Internal::Traits<
      CanonicalComponentFieldBinary<Derived0, Derived1, Function>>::Scalar;
  using Value = Internal::Traits<
      CanonicalComponentFieldBinary<Derived0, Derived1, Function>>::Value;
  using Writeable = Internal::Traits<
      CanonicalComponentFieldBinary<Derived0, Derived1, Function>>::Writeable;

  auto UpperIndex() const { return _N; }
  auto& UpperIndex() { return _N; }
  auto& Grid() const { return _u0.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u0[iTheta, iPhi], _u1[iTheta, iPhi]);
  }

  // Constructors.
  CanonicalComponentFieldBinary() = delete;
  CanonicalComponentFieldBinary(const CanonicalComponentFieldBase<Derived0>& u0,
                                const CanonicalComponentFieldBase<Derived1>& u1,
                                Function&& f, Int N)
      : _u0{u0}, _u1{u1}, _f{f}, _N{N} {}

  CanonicalComponentFieldBinary(const CanonicalComponentFieldBinary&) = default;
  CanonicalComponentFieldBinary(CanonicalComponentFieldBinary&&) = default;

  // Assignment.
  CanonicalComponentFieldBinary& operator=(CanonicalComponentFieldBinary&) =
      default;
  CanonicalComponentFieldBinary& operator=(CanonicalComponentFieldBinary&&) =
      default;

 private:
  const CanonicalComponentFieldBase<Derived0>& _u0;
  const CanonicalComponentFieldBase<Derived1>& _u1;
  const Function& _f;
  Int _N;
};

}  // namespace GSHTrans

#endif