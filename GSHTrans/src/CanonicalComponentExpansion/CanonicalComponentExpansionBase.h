#ifndef GSH_TRANS_CANONICAL_COMPONENT_EXPANSION_BASE_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_EXPANSION_BASE_GUARD_H

#include <concepts>
#include <format>
#include <iostream>
#include <type_traits>

#include "../Concepts.h"
#include "../ExpansionBase.h"
#include "../GridBase.h"
#include "../Indexing.h"

namespace GSHTrans {

template <std::ptrdiff_t _N, typename _Derived>
class CanonicalComponentExpansionBase
    : public ExpansionBase<CanonicalComponentExpansionBase<_N, _Derived>> {
 public:
  using Int = typename Internal::Traits<_Derived>::Int;
  using Value = typename Internal::Traits<_Derived>::Value;
  using Real = typename Internal::Traits<_Derived>::Real;
  using Complex = typename Internal::Traits<_Derived>::Complex;
  using Scalar = typename Internal::Traits<_Derived>::Scalar;
  using MRange = typename Internal::Traits<_Derived>::MRange;
  using Writeable = typename Internal::Traits<_Derived>::Writeable;

  // Return the upper index.
  constexpr auto UpperIndex() const { return _N; }

  // Return the minimum degree.
  auto MinimumDegree() const { return std::abs(UpperIndex()); }

  // Return the grid.
  auto& Grid() const { return Derived().Grid(); }

  // Return spherical harmonic indices.
  auto Indices() const {
    return GSHIndices<MRange>(this->MaxDegree(), this->MaxDegree(), _N)
        .Indices();
  }

  // Return the index for the (l,m)th coefficient.
  auto Index(Int l, Int m) const {
    return GSHIndices<MRange>(this->MaxDegree(), this->MaxDegree(), _N)
        .Index(l, m);
  }

  // Return the spherical harmonic degree for each datum.
  auto Degrees() const { return Indices() | std::ranges::views::keys; }

  // Return the spherical harmonic order for each datum.
  auto Orders() const { return Indices() | std::ranges::views::values; }

  // Return the total number of coefficients.
  auto Size() const {
    return GSHIndices<MRange>(this->MaxDegree(), this->MaxDegree(), _N).Size();
  }

  // Read access to data.
  auto operator[](Int l, Int m) const { return Derived()[l, m]; }

  // Write access to data.
  auto& operator[](Int l, Int m)
  requires Writeable::value
  {
    return Derived()[l, m];
  }

  // Assign values from another Expansion.
  template <typename __Derived>
  requires Writeable::value && std::same_as<typename __Derived::Scalar, Scalar>
  auto& operator=(const CanonicalComponentExpansionBase<_N, __Derived>& other) {
    assert(other.MaxDegree() == this->MaxDegree());
    for (auto [l, m] : this->Indices()) {
      operator[](l, m) = other[l, m];
    }
    return Derived();
  }

  template <typename __Derived>
  requires Writeable::value && std::same_as<typename __Derived::Scalar, Scalar>
  auto& operator=(CanonicalComponentExpansionBase<_N, __Derived>&& other) {
    Derived() = other;
    return Derived();
  }

  // Compound plus assignment with other Expansion.
  template <typename __Derived>
  requires Writeable::value && std::same_as<typename __Derived::Scalar, Scalar>
  auto& operator+=(
      const CanonicalComponentExpansionBase<_N, __Derived>& other) {
    assert(other.MaxDegree() == this->MaxDegree());
    for (auto [l, m] : this->Indices()) {
      operator[](l, m) += other[l, m];
    }
    return Derived();
  }

  template <typename __Derived>
  requires Writeable::value && std::same_as<typename __Derived::Scalar, Scalar>
  auto& operator+=(CanonicalComponentExpansionBase<_N, __Derived>&& other) {
    Derived() += other;
    return Derived();
  }

  // Compound minus assignment with other Expansion.
  template <typename __Derived>
  requires Writeable::value && std::same_as<typename __Derived::Scalar, Scalar>
  auto& operator-=(
      const CanonicalComponentExpansionBase<_N, __Derived>& other) {
    assert(other.MaxDegree() == this->MaxDegree());
    for (auto [l, m] : this->Indices()) {
      operator[](l, m) -= other(l, m);
    }
    return Derived();
  }

  template <typename __Derived>
  requires Writeable::value && std::same_as<typename __Derived::Scalar, Scalar>
  auto& operator-=(CanonicalComponentExpansionBase<_N, __Derived>&& other) {
    Derived() -= other;
    return Derived();
  }

  // Compound multiply assignment with scalar.
  auto& operator*=(Scalar s)
  requires Writeable::value
  {
    for (auto [l, m] : this->Indices()) {
      operator[](l, m) *= s;
    }
    return Derived();
  }

  // Compound divide assignment with scalar.
  auto& operator/=(Scalar s)
  requires Writeable::value
  {
    for (auto [l, m] : this->Indices()) {
      operator[](l, m) /= s;
    }
    return Derived();
  }

  // Write values to ostream.
  friend std::ostream& operator<<(
      std::ostream& os,
      const CanonicalComponentExpansionBase<_N, _Derived>& u) {
    auto values = u.Indices() | std::ranges::views::transform([&u](auto pair) {
                    auto [l, m] = pair;
                    return u[l, m];
                  });
    auto range = std::ranges::views::zip(u.Degrees(), u.Orders(), values);
    auto n = u.Size();
    auto printValues = [&os](auto l, auto m, auto val, auto newLine) {
      os << std::format("{:+8}  {:+8}  ({:+.8e},  {:+.8e})", l, m, val.real(),
                        val.imag());
      if (newLine) os << '\n';
    };
    for (auto [l, m, val] : range | std::ranges::views::take(n - 1)) {
      printValues(l, m, val, true);
    }
    for (auto [l, m, val] : range | std::ranges::views::drop(n - 1)) {
      printValues(l, m, val, false);
    }
    return os;
  }

 private:
  constexpr auto& Derived() const {
    return static_cast<const _Derived&>(*this);
  }
  constexpr auto& Derived() { return static_cast<_Derived&>(*this); }

  // Return the GSHIndex.
  auto GSHIndex() const {
    return std::move(
        GSHIndices<MRange>(this->MaxDegree(), this->MaxDegree(), _N));
  }
};

}  // namespace GSHTrans

#endif