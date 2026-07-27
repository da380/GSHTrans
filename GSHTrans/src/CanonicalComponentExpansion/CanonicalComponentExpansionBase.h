#ifndef GSH_TRANS_CANONICAL_COMPONENT_EXPANSION_BASE_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_EXPANSION_BASE_GUARD_H

#include <cassert>
#include <concepts>
#include <cstddef>

#include "../Concepts.h"
#include "../GridBase.h"
#include "../Indexing.h"

namespace GSHTrans {

template <std::ptrdiff_t _N, typename _Derived>
class CanonicalComponentExpansionBase {
 public:
  using Int = typename Internal::Traits<_Derived>::Int;
  using Value = typename Internal::Traits<_Derived>::Value;
  using Real = typename Internal::Traits<_Derived>::Real;
  using Complex = typename Internal::Traits<_Derived>::Complex;
  using Scalar = typename Internal::Traits<_Derived>::Scalar;
  using MRange = typename Internal::Traits<_Derived>::MRange;
  using Writeable = typename Internal::Traits<_Derived>::Writeable;

  constexpr auto UpperIndex() const { return _N; }

  auto MinDegree() const { return Derived().MinDegree(); }
  auto MaxDegree() const { return Derived().MaxDegree(); }
  auto MaxOrder() const { return Derived().MaxOrder(); }
  auto Degrees() const { return Derived().Degrees(); }
  auto Indices() const { return Derived().Indices(); }
  auto Orders() const { return Derived().Orders(); }
  auto Size() const { return Derived().Size(); }
  auto Index(Int l, Int m) const { return Derived().Index(l, m); }

  auto& Grid() const { return Derived().Grid(); }

  auto operator[](Int l, Int m) const { return Derived()[l, m]; }

  auto& operator[](Int l, Int m)
  requires Writeable::value
  {
    return Derived()[l, m];
  }

  template <typename __Derived>
  requires Writeable::value &&
           std::same_as<typename __Derived::Scalar, Scalar> &&
           std::same_as<typename __Derived::MRange, MRange>
  auto& operator=(const CanonicalComponentExpansionBase<_N, __Derived>& other) {
    assert(other.MaxDegree() == MaxDegree());
    for (auto [l, m] : Indices()) {
      operator[](l, m) = other[l, m];
    }
    return Derived();
  }

  template <typename __Derived>
  requires Writeable::value &&
           std::same_as<typename __Derived::Scalar, Scalar> &&
           std::same_as<typename __Derived::MRange, MRange>
  auto& operator=(CanonicalComponentExpansionBase<_N, __Derived>&& other) {
    return operator=(other);
  }

  template <typename __Derived>
  requires Writeable::value &&
           std::same_as<typename __Derived::Scalar, Scalar> &&
           std::same_as<typename __Derived::MRange, MRange>
  auto& operator+=(
      const CanonicalComponentExpansionBase<_N, __Derived>& other) {
    assert(other.MaxDegree() == MaxDegree());
    for (auto [l, m] : Indices()) {
      operator[](l, m) += other[l, m];
    }
    return Derived();
  }

  template <typename __Derived>
  requires Writeable::value &&
           std::same_as<typename __Derived::Scalar, Scalar> &&
           std::same_as<typename __Derived::MRange, MRange>
  auto& operator+=(CanonicalComponentExpansionBase<_N, __Derived>&& other) {
    return operator+=(other);
  }

  template <typename __Derived>
  requires Writeable::value &&
           std::same_as<typename __Derived::Scalar, Scalar> &&
           std::same_as<typename __Derived::MRange, MRange>
  auto& operator-=(
      const CanonicalComponentExpansionBase<_N, __Derived>& other) {
    assert(other.MaxDegree() == MaxDegree());
    for (auto [l, m] : Indices()) {
      operator[](l, m) -= other[l, m];
    }
    return Derived();
  }

  template <typename __Derived>
  requires Writeable::value &&
           std::same_as<typename __Derived::Scalar, Scalar> &&
           std::same_as<typename __Derived::MRange, MRange>
  auto& operator-=(CanonicalComponentExpansionBase<_N, __Derived>&& other) {
    return operator-=(other);
  }

  auto& operator*=(Scalar scalar)
  requires Writeable::value
  {
    for (auto [l, m] : Indices()) {
      operator[](l, m) *= scalar;
    }
    return Derived();
  }

  auto& operator/=(Scalar scalar)
  requires Writeable::value
  {
    for (auto [l, m] : Indices()) {
      operator[](l, m) /= scalar;
    }
    return Derived();
  }

 private:
  constexpr auto& Derived() const {
    return static_cast<const _Derived&>(*this);
  }
  constexpr auto& Derived() { return static_cast<_Derived&>(*this); }
};

}  // namespace GSHTrans

#endif
