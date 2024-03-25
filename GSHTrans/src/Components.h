#ifndef GSH_TRANS_COMPONENTS_GUARD_HVector
#define GSH_TRANS_COMPONENTS_GUARD_H

#include <FFTWpp/Core>
#include <algorithm>
#include <complex>
#include <functional>
#include <iostream>
#include <ranges>

#include "Concepts.h"
#include "GridBase.h"

namespace GSHTrans {

namespace Testing {

template <typename S>
requires std::integral<S> || RealOrComplexFloatingPoint<S>
class ScalarMultiply {
 public:
  ScalarMultiply() = default;
  ScalarMultiply(S s) : _s{s} {}

  ScalarMultiply(const ScalarMultiply&) = default;
  ScalarMultiply(ScalarMultiply&&) = default;

  template <typename T>
  auto operator()(T t) const {
    return t * _s;
  }

 private:
  S _s;
};

struct Field {};
struct Coefficient {};
struct RealValued {};
struct ComplexValued {};

template <typename T>
concept ComponentType = std::same_as<T, Field> || std::same_as<T, Coefficient>;

template <typename T>
concept ComponentValue =
    std::same_as<T, RealValued> || std::same_as<T, ComplexValued>;

template <typename _Derived>
class ComponentBase {
 public:
  using Int = std::ptrdiff_t;

  // Methods related to the grid.
  auto GetGrid() const { return Derived().GetGrid(); }
  auto NumberOfLongitudes() const
  requires std::same_as<typename _Derived::Type, Field>
  {
    return GetGrid().NumberOfLongitudes();
  }
  auto NumberOfLatitudes() const
  requires std::same_as<typename _Derived::Type, Field>
  {
    return GetGrid().NumberOfLatitudes();
  }

  // Methods related to the canonical indices.
  auto CanonicalIndices() const { return Derived().CanonicalIndices(); }
  auto UpperIndex() const {
    return std::ranges::fold_left_first(CanonicalIndices(), std::plus<>())
        .value_or(Int(0));
  }

  // Methods related to the data.

  auto GetSize() const {
    using Type = typename _Derived::Type;
    using Value = typename _Derived::Value;
    if constexpr (std::same_as<Type, Field>) {
      return GetGrid().ComponentSize();
    } else {
      if (std::same_as<Value, RealValued> && UpperIndex() == 0) {
        return GetGrid().CoefficientSizeNonNegative(UpperIndex());
      } else {
        return GetGrid().CoefficientSize(UpperIndex());
      }
    }
  }

  // Application operators for fields.
  auto operator()(Int iTheta, Int iPhi) const
  requires std::same_as<typename _Derived::Type, Field>
  {
    auto i = iTheta * NumberOfLongitudes() + iPhi;
    return operator[](i);
  }
  auto& operator()(Int iTheta, Int iPhi) const
  requires std::same_as<typename _Derived::Type, Field> &&
           std::ranges::output_range<
               typename _Derived::View,
               std::ranges::range_value_t<typename _Derived::View>>
  {
    auto i = iTheta * NumberOfLongitudes() + iPhi;
    return operator[](i);
  }

 private:
  auto& Derived() const { return static_cast<const _Derived&>(*this); }
  auto& Derived() { return static_cast<_Derived&>(*this); }
};

// Component class that stores its own data.
template <typename _Grid, ComponentType _Type, ComponentValue _Value>
class Component : public ComponentBase<Component<_Grid, _Type, _Value>> {
 public:
  // Public type aliases.
  using Int = std::ptrdiff_t;
  using Grid = _Grid;
  using Type = _Type;
  using Value = _Value;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar = std::conditional_t<
      std::same_as<_Type, Field>,
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>,
      Complex>;
  using Vector = std::vector<Scalar>;
  using MRange = typename _Grid::MRange;
  using NRange = typename _Grid::NRange;
  using View = std::ranges::views::all_t<Vector>;

  // Methods required for the base class.
  auto GetGrid() const { return _grid; }
  auto CanonicalIndices() const { return _indices; }

  // Constructors.
  Component() = default;

  Component(_Grid grid, std::vector<Int>&& indices)
      : _grid{grid},
        _indices{indices},
        _data{std::vector<Scalar>(this->GetSize())} {}

  template <typename... Indices>
  requires(std::integral<Indices> && ...)
  Component(_Grid grid, Indices... indices)
      : _grid{grid},
        _indices{std::vector<Int>{indices...}},
        _data{std::vector<Scalar>(this->GetSize())} {}

 private:
  _Grid _grid;
  std::vector<Int> _indices;
  std::vector<Scalar> _data;
};

}  // namespace Testing

/*



template <typename Derived>
class ComponentBase {
  using Int = std::ptrdiff_t;

 public:
  auto Grid() const { return _Derived().Grid(); }
  auto CanonicalIndices() const { return _Derived().CanonicalIndices(); }
  auto UpperIndex() const {
    return std::ranges::fold_left(CanonicalIndices(), Int(0), std::plus<>());
  }

  auto operator[](Int i) const { return _Derived()[i]; }
  auto& operator[](Int i) { return _Derived()[i]; }

 private:
  auto& _Derived() const { return static_cast<const Derived&>(*this); }
  auto& _Derived() { return static_cast<Derived&>(*this); }
};

// Component options.
struct RealField {};
struct ComplexField {};
struct Coefficient {};

template <typename T>
concept ComponentType =
    std::same_as<T, RealField> or std::same_as<T, ComplexField> or
    std::same_as<T, Coefficient>;

template <typename _GridType, ComponentType _Type>
requires requires() {
  requires std::derived_from<_GridType, GridBase<_GridType>>;
}
class Component : public ComponentBase<Component<_GridType, _Type>> {
  using Int = std::ptrdiff_t;

 public:
  using GridType = _GridType;
  using Type = _Type;

  using Real = typename _GridType::real_type;
  using Complex = typename _GridType::complex_type;
  using Scalar =
      std::conditional_t<std::same_as<_Type, RealField>, Real, Complex>;

  Component() = default;

  template <typename... Indices>
  requires(std::integral<Indices> && ...)
  Component(_GridType grid, Indices... indices)
      : _grid{grid},
        _indices{std::vector<Int>{indices...}},
        _data{std::vector<Scalar>(_grid.ComponentSize())} {}

  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  auto CanonicalIndices() const { return std::ranges::views::all(_indices); }

 private:
  _GridType _grid;
  std::vector<Int> _indices;
  std::vector<Scalar> _data;
};

*/

}  // namespace GSHTrans

#endif  // GSH_TRANS_COMPONENTS_GUARD_H