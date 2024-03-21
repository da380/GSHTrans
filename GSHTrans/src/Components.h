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

template <typename _Derived>
class ComponentBase {
  using Int = std::ptrdiff_t;

 public:
  auto GetGrid() const { return Derived().GetGrid(); }
  auto CanonicalIndices() const { return Derived().CanonicalIndices(); }
  auto UpperIndex() const {
    return std::ranges::fold_left_first(CanonicalIndices(), std::plus<>())
        .value_or(Int(0));
  }

  auto Data() { return Derived().Data(); }
  auto Data() const { return Derived().Data(); }
  auto size() const { return Derived().Data().size(); }
  auto begin() { return Derived().Data().begin(); }
  auto end() { return Derived().Data().end(); }
  auto cbegin() const { return Derived().Data().cbegin(); }
  auto cend() const { return Derived().Data().cend(); }

  auto operator[](Int i) const { return Derived().Data()[i]; }
  auto& operator[](Int i) { return Derived().Data()[i]; }

  auto GetSize() const {
    using Type = typename _Derived::Type;
    using Value = typename _Derived::Value;
    if constexpr (std::same_as<Type, Field>) {
      return GetGrid().ComponentSize();
    } else {
      if constexpr (std::same_as<Value, RealValued>) {
        return GetGrid().RealCoefficientSize(UpperIndex());
      } else {
        return GetGrid().ComplexCoefficientSize(UpperIndex());
      }
    }
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
  using MRange = typename _Grid::MRange;
  using NRange = typename _Grid::NRange;

  // Methods required for the base class.
  auto GetGrid() const { return _grid; }
  auto CanonicalIndices() const { return _indices; }
  auto Data() { return std::ranges::views::all(_data); }
  auto Data() const { return std::ranges::views::all(_data); }

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

// Component class with a view to its data.
template <typename _Grid, ComponentType _Type, ComponentValue _Value,
          std::ranges::view _View>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::ranges::input_range<_View>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
}
class ComponentView
    : public ComponentBase<ComponentView<_Grid, _Type, _Value, _View>> {
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
  using MRange = typename _Grid::MRange;
  using NRange = typename _Grid::NRange;

  // Methods required for the base class.
  auto GetGrid() const { return _grid; }
  auto CanonicalIndices() const { return _indices; }
  auto Data() { return _data; }
  auto Data() const { return _data; }

  // Constructors.
  ComponentView() = default;

  ComponentView(_Grid grid, _View data, std::vector<Int>&& indices)
      : _grid{grid}, _indices{indices}, _data{data} {
    assert(_data.size() == this->GetSize());
  }

  template <typename... Indices>
  requires(std::integral<Indices> && ...)
  ComponentView(_Grid grid, _View data, Indices... indices)
      : _grid{grid}, _indices{std::vector<Int>{indices...}}, _data{data} {
    assert(_data.size() == this->GetSize());
  }

 private:
  _Grid _grid;
  std::vector<Int> _indices;
  _View _data;
};

// Component class with data formed a unary transformation of a view.
template <typename _Grid, ComponentType _Type, ComponentValue _Value,
          std::ranges::view _View, typename _Function>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::ranges::input_range<_View>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                        typename _Grid::Real>;
}
class ComponentViewUnary
    : public ComponentBase<
          ComponentViewUnary<_Grid, _Type, _Value, _View, _Function>> {
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
  using MRange = typename _Grid::MRange;
  using NRange = typename _Grid::NRange;

  // Methods required for the base class.
  auto GetGrid() const { return _grid; }
  auto CanonicalIndices() const { return _indices; }
  // auto Data() { return _data | std::ranges::views::transform(_f); }
  auto Data() const { return _data | std::ranges::views::transform(_f); }

  // Constructors.
  ComponentViewUnary() = default;

  ComponentViewUnary(_Grid grid, _View data, _Function f,
                     std::vector<Int>&& indices)
      : _grid{grid}, _indices{indices}, _data{data}, _f{f} {
    assert(_data.size() == this->GetSize());
  }

 private:
  _Grid _grid;
  std::vector<Int> _indices;
  _View _data;
  _Function _f;
};

// Component class with data formed a binary transformation of two views.
template <typename _Grid, ComponentType _Type, ComponentValue _Value,
          std::ranges::view _View1, std::ranges::view _View2,
          typename _Function>
requires requires() {
  requires std::derived_from<_Grid, GridBase<_Grid>>;
  requires std::ranges::input_range<_View1>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View1>>,
                        typename _Grid::Real>;
  requires std::ranges::input_range<_View2>;
  requires std::same_as<RemoveComplex<std::ranges::range_value_t<_View2>>,
                        typename _Grid::Real>;
}
class ComponentViewBinary
    : public ComponentBase<ComponentViewBinary<_Grid, _Type, _Value, _View1,
                                               _View2, _Function>> {
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
  using MRange = typename _Grid::MRange;
  using NRange = typename _Grid::NRange;

  // Methods required for the base class.
  auto GetGrid() const { return _grid; }
  auto CanonicalIndices() const { return _indices; }
  auto Data() const {
    return std::ranges::views::zip_transform(_f, _data1, _data2);
  }

  // Constructors.
  ComponentViewBinary() = default;

  ComponentViewBinary(_Grid grid, _View1 data1, _View2 data2, _Function f,
                      std::vector<Int>&& indices)
      : _grid{grid}, _indices{indices}, _data1{data1}, _data2{data2}, _f{f} {
    assert(_data1.size() == this->GetSize());
    assert(_data2.size() == this->GetSize());
  }

 private:
  _Grid _grid;
  std::vector<Int> _indices;
  _View1 _data1;
  _View2 _data2;
  _Function _f;
};

template <typename Derived>
auto operator-(ComponentBase<Derived>& u) {
  auto data = u.Data();
  auto f = [](auto x) { return -x; };
  using View = decltype(data);
  using Function = decltype(f);
  return ComponentViewUnary<typename Derived::Grid, typename Derived::Type,
                            typename Derived::Value, View, Function>(
      u.GetGrid(), data, f, u.CanonicalIndices());
}

template <typename Derived1, typename Derived2>
auto operator+(ComponentBase<Derived1>& u, ComponentBase<Derived2>& v) {
  auto data1 = u.Data();
  auto data2 = v.Data();
  auto f = std::plus<>();
  return ComponentViewBinary<typename Derived1::Grid, typename Derived1::Type,
                             typename Derived1::Value, decltype(data1),
                             decltype(data2), decltype(f)>(
      u.GetGrid(), data1, data2, f, u.CanonicalIndices());
}

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