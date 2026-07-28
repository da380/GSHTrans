#ifndef GSH_TRANS_CANONICAL_COMPONENT_EXPANSION_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_EXPANSION_GUARD_H

#include <FFTWpp/Core>
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <stdexcept>
#include <type_traits>

#include "../Concepts.h"
#include "CanonicalComponentExpansionBase.h"

namespace GSHTrans {

template <std::ptrdiff_t _N, typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class CanonicalComponentExpansion;

namespace Internal {

template <std::ptrdiff_t _N, typename _Grid, RealOrComplexValued _Value>
struct Traits<CanonicalComponentExpansion<_N, _Grid, _Value>> {
  using Int = std::ptrdiff_t;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;

  // RealValued describes real grid samples represented by the reduced m >= 0
  // spectrum. Its stored coefficients remain complex; negative Fourier modes
  // are implicit and are related to the expansion at -n by GSH symmetry.
  using Scalar = Complex;
  using MRange =
      std::conditional_t<std::same_as<Value, RealValued>, NonNegative, All>;
  using Writeable = std::true_type;
};

}  // namespace Internal

template <std::ptrdiff_t _N, typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class CanonicalComponentExpansion
    : public CanonicalComponentExpansionBase<
          _N, CanonicalComponentExpansion<_N, _Grid, _Value>> {
 public:
  using Self = CanonicalComponentExpansion<_N, _Grid, _Value>;
  using Base = CanonicalComponentExpansionBase<_N, Self>;
  using Traits = Internal::Traits<Self>;
  using Value = typename Traits::Value;
  using Int = typename Traits::Int;
  using Real = typename Traits::Real;
  using Complex = typename Traits::Complex;
  using Scalar = typename Traits::Scalar;
  using MRange = typename Traits::MRange;
  using Writeable = typename Traits::Writeable;

  auto& Grid() const { return _grid; }

  auto MinDegree() const { return _indices.MinDegree(); }
  auto MaxDegree() const { return _indices.MaxDegree(); }
  auto MaxOrder() const { return _indices.MaxOrder(); }
  auto Degrees() const { return _indices.Degrees(); }
  auto Indices() const { return _indices.Indices(); }
  auto Orders() const { return Indices() | std::ranges::views::values; }
  auto Size() const { return _indices.Size(); }
  auto Index(Int l, Int m) const { return _indices.Index(l, m); }

  auto operator[](Int l, Int m) const { return _data[Index(l, m)]; }
  auto& operator[](Int l, Int m) { return _data[Index(l, m)]; }

  auto Data() { return std::ranges::views::all(_data); }
  auto Data() const { return std::ranges::views::all(_data); }

  CanonicalComponentExpansion() = delete;

  explicit CanonicalComponentExpansion(_Grid& grid)
      : _grid{grid},
        _indices{grid.MaxDegree(), grid.MaxDegree(), _N},
        _data{FFTWpp::vector<Complex>(_indices.Size())} {}

  template <typename Derived>
  requires std::same_as<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::MRange, MRange>
  CanonicalComponentExpansion(
      const CanonicalComponentExpansionBase<_N, Derived>& other)
      : CanonicalComponentExpansion(other.Grid()) {
    for (auto [l, m] : Indices()) {
      operator[](l, m) = other[l, m];
    }
  }

  template <typename Derived>
  requires std::same_as<typename Derived::Scalar, Scalar> &&
           std::same_as<typename Derived::MRange, MRange>
  CanonicalComponentExpansion(
      CanonicalComponentExpansionBase<_N, Derived>&& other)
      : CanonicalComponentExpansion(other) {}

  CanonicalComponentExpansion(const CanonicalComponentExpansion&) = default;
  CanonicalComponentExpansion(CanonicalComponentExpansion&&) = default;

  CanonicalComponentExpansion& operator=(
      const CanonicalComponentExpansion& other) {
    return AssignValues(other);
  }
  CanonicalComponentExpansion& operator=(CanonicalComponentExpansion&& other) {
    return AssignValues(other);
  }

  using Base::operator=;

 private:
  _Grid& _grid;
  GSHIndices<MRange> _indices;
  FFTWpp::vector<Complex> _data;

  auto& AssignValues(const CanonicalComponentExpansion& other) {
    if (other.Size() != Size()) {
      throw std::invalid_argument(
          "Cannot assign canonical component expansions with different sizes");
    }
    std::ranges::copy(other._data, _data.begin());
    return *this;
  }
};

template <std::ptrdiff_t N, typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealCanonicalComponentExpansion =
    CanonicalComponentExpansion<N, Grid, RealValued>;

template <std::ptrdiff_t N, typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexCanonicalComponentExpansion =
    CanonicalComponentExpansion<N, Grid, ComplexValued>;

template <typename Grid, RealOrComplexValued Value>
requires std::derived_from<Grid, GridBase<Grid>>
using ScalarExpansion = CanonicalComponentExpansion<0, Grid, Value>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealScalarExpansion =
    CanonicalComponentExpansion<0, Grid, RealValued>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexScalarExpansion =
    CanonicalComponentExpansion<0, Grid, ComplexValued>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COMPONENT_EXPANSION_GUARD_H
