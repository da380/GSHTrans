#ifndef GSH_TRANS_COMPONENTS_GUARD_HVector
#define GSH_TRANS_COMPONENTS_GUARD_H

#include <FFTWpp/Core>
#include <algorithm>
#include <complex>
#include <iostream>
#include <ranges>

#include "Concepts.h"
#include "GridBase.h"

namespace GSHTrans {

namespace Testing {

struct Field {};
struct Coefficient {};
struct RealValued {};
struct ComplexValued {};

template <typename T>
concept ComponentType = std::same_as<T, Field> || std::same_as<T, Coefficient>;

template <typename T>
concept ComponentValue =
    std::same_as<T, RealValued> || std::same_as<T, ComplexValued>;

template <typename Derived>
class ComponentBase {
  using Int = std::ptrdiff_t;

 public:
  auto CanonicalIndices() const { return _Derived().CanonicalIndices(); }
  auto UpperIndex() const {
    return std::ranges::fold_left_first(CanonicalIndices(), std::plus<>())
        .value_or(Int(0));
  }

  auto size() const { return _Derived().All().size(); }
  auto begin() { return _Derived().All().begin(); }
  auto end() { return _Derived().All().end(); }
  auto cbegin() const { return _Derived().All().cbegin(); }
  auto cend() const { return _Derived().All().cend(); }

  auto operator[](Int i) const { return _Derived().All()[i]; }
  auto& operator[](Int i) { return _Derived().All()[i]; }

 private:
  auto& _Derived() const { return static_cast<const Derived&>(*this); }
  auto& _Derived() { return static_cast<Derived&>(*this); }
};

template <typename _Grid, ComponentType _Type, ComponentValue _Value>
class Component : public ComponentBase<Component<_Grid, _Type, _Value>> {
  using Int = std::ptrdiff_t;

 public:
  using Grid = _Grid;
  using Type = _Type;
  using Value = _Value;

  using Real = typename _Grid::real_type;
  using Complex = typename _Grid::complex_type;

  using Scalar = std::conditional_t<
      std::same_as<_Type, Field>,
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>,
      Complex>;

  using MRange = typename _Grid::MRange_type;

  Component() = default;

  template <typename... Indices>
  requires(std::integral<Indices> && ...)
  Component(_Grid grid, Indices... indices)
      : _grid{grid},
        _indices{std::vector<Int>{indices...}},
        _data{std::vector<Scalar>(GetSize())} {
    std::cout << _data.size() << std::endl;
  }

  auto CanonicalIndices() const { return _indices; }

  auto All() { return std::ranges::views::all(_data); }
  auto All() const { return std::ranges::views::all(_data); }

 private:
  _Grid _grid;
  std::vector<Int> _indices;
  std::vector<Scalar> _data;

  auto GetSize() const {
    if constexpr (std::same_as<_Type, Field>) {
      return _grid.ComponentSize();
    } else {
      if constexpr (std::same_as<_Value, RealValued>) {
        return _grid.RealCoefficientSize(this->UpperIndex());
      } else {
        return _grid.ComplexCoefficientSize(this->UpperIndex());
      }
    }
  }
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