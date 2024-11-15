#ifndef GSH_TRANS_CANONICAL_COMPONENT_EXPANSION_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_EXPANSION_GUARD_H

#include <concepts>
#include <format>
#include <type_traits>

#include "../Concepts.h"
#include "CanonicalComponentExpansionBase.h"

namespace GSHTrans {

// Forward declare class.
template <std::ptrdiff_t _N, typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class CanonicalComponentExpansion;

// Set traits.
namespace Internal {

template <std::ptrdiff_t _N, typename _Grid, RealOrComplexValued _Value>
struct Traits<CanonicalComponentExpansion<_N, _Grid, _Value>> {
  using Int = std::ptrdiff_t;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<Value, RealValued>, Real, Complex>;
  using MRange = typename std::conditional_t<std::same_as<Value, RealValued>,
                                             NonNegative, All>;
  using Writeable = std::true_type;
};

}  // namespace Internal

template <std::ptrdiff_t _N, typename _Grid, RealOrComplexValued _Value>
requires std::derived_from<_Grid, GridBase<_Grid>>
class CanonicalComponentExpansion
    : public CanonicalComponentExpansionBase<
          _N, CanonicalComponentExpansion<_N, _Grid, _Value>>,
      public GSHIndices<std::conditional_t<std::same_as<_Value, RealValued>,
                                           NonNegative, All>> {
 public:
  using Value = typename Internal::Traits<
      CanonicalComponentExpansion<_N, _Grid, _Value>>::Value;
  using Int = typename Internal::Traits<
      CanonicalComponentExpansion<_N, _Grid, _Value>>::Int;
  using Real = typename Internal::Traits<
      CanonicalComponentExpansion<_N, _Grid, _Value>>::Real;
  using Complex = typename Internal::Traits<
      CanonicalComponentExpansion<_N, _Grid, _Value>>::Complex;
  using Scalar = typename Internal::Traits<
      CanonicalComponentExpansion<_N, _Grid, _Value>>::Scalar;
  using MRange = typename Internal::Traits<
      CanonicalComponentExpansion<_N, _Grid, _Value>>::MRange;
  using Writeable = typename Internal::Traits<
      CanonicalComponentExpansion<_N, _Grid, _Value>>::Writeable;

  // Return the grid.
  auto& Grid() const { return _grid; }

  // Read access to data.
  auto operator[](Int l, Int m) const { return _data[this->Index(l, m)]; }

  // Write access to data.
  auto& operator[](Int l, Int m) { return _data[this->Index(l, m)]; }

  // Return a view to the data.
  auto Data() { return std::ranges::views::all(_data); }

  // Default constructor.
  CanonicalComponentExpansion() = default;

  // Construct from grid initialising values to zero.
  CanonicalComponentExpansion(_Grid& grid)
      : GSHIndices<MRange>(grid.MaxDegree(), grid.MaxDegree(), _N),
        _grid{grid},
        _data{FFTWpp::vector<Complex>(this->Size())} {}

  /*



  // Default constructor.
  CanonicalComponentExpansion() = default;

  // Construct from grid initialising values to zero.
  CanonicalComponentExpansion(_Grid& grid)
      : _grid{grid}, _data{FFTWpp::vector<Scalar>(this->ExpansionSize())} {}

  // Construction from grid initialising values with a function.
  template <typename Function>
  requires ScalarFunctionS2<Function, Real, Scalar>
  CanonicalComponentExpansion(_Grid& grid, Function&& f)
      : CanonicalComponentExpansion(grid) {
    std::ranges::copy(_grid.ProjectFunction(f), _data.begin());
  }

  // Construct from an element of the base class.
  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  CanonicalComponentExpansion(
      const CanonicalComponentExpansionBase<_N, Derived>& other)
      : CanonicalComponentExpansion(other.UpperIndex(), other.Grid()) {
    for (auto [iTheta, iPhi] : this->PointIndices()) {
      operator[](iTheta, iPhi) = other[iTheta, iPhi];
    }
  }

  template <typename Derived>
  requires std::convertible_to<typename Derived::Scalar, Scalar>
  CanonicalComponentExpansion(
      CanonicalComponentExpansionBase<_N, Derived>&& other)
      : CanonicalComponentExpansion(other) {}

  // Default copy and move constructors.
  CanonicalComponentExpansion(const CanonicalComponentExpansion&) = default;
  CanonicalComponentExpansion(CanonicalComponentExpansion&&) = default;

  // Default copy and move assignment.
  CanonicalComponentExpansion& operator=(const CanonicalComponentExpansion&)
  = default; CanonicalComponentExpansion&
  operator=(CanonicalComponentExpansion&&) = default;

  // Use assignment defined in base class.
  using CanonicalComponentExpansionBase<
      _N, CanonicalComponentExpansion<_N, _Grid, _Value>>::operator=;

*/

 private:
  _Grid& _grid;
  FFTWpp::vector<Complex> _data;
};

/*

// Type aliases for real and complex Expansions.
template <std::ptrdiff_t N, typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealCanonicalComponentExpansion =
    CanonicalComponentExpansion<N, Grid, RealValued>;

template <std::ptrdiff_t N, typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexCanonicalComponentExpansion =
    CanonicalComponentExpansion<N, Grid, ComplexValued>;

// Type aliases for scalar Expansions.
template <typename Grid, RealOrComplexValued Value>
requires std::derived_from<Grid, GridBase<Grid>>
using ScalarExpansion = CanonicalComponentExpansion<0, Grid, Value>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using RealScalarExpansion = CanonicalComponentExpansion<0, Grid, RealValued>;

template <typename Grid>
requires std::derived_from<Grid, GridBase<Grid>>
using ComplexScalarExpansion =
    CanonicalComponentExpansion<0, Grid, ComplexValued>;

*/

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COMPONENT_EXPANSION_GUARD_H