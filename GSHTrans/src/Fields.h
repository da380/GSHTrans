#ifndef GSH_TRANS_FIELDS_GUARD_H
#define GSH_TRANS_FIELDS_GUARD_H

#include <FFTWpp/All>
#include <algorithm>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>
#include <numeric>
#include <ranges>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

// Base class for canonical components.
template <typename Derived>
class CanonicalComponentBase {
  using Int = std::ptrdiff_t;

 public:
  auto& GridReference() const { return GetDerived()._GridReference(); }
  auto GridPointer() const { return GetDerived()._GridPointer(); }

  auto MaxDegree() const { return GridReference().MaxDegree(); }
  auto MaxUpperIndex() const { return GridReference().MaxUpperIndex(); }

  auto CoLatitudes() const { return GridReference().CoLatitudes(); }
  auto NumberOfCoLatitudes() const {
    return GridReference().NumberOfCoLatitudes();
  }
  auto CoLatitudeIndices() const { return GridReference().CoLatitudeIndices(); }

  auto LongitudeSpacing() const { return GridReference().LongitudeSpacing(); }
  auto LongitudeIndices() const { return GridReference().LongitudeIndices(); }
  auto Longitudes() const { return GridReference().Lognitudes(); }
  auto NumberOfLongitudes() const {
    return GridReference().NumberOfLongitudes();
  }

  auto View() { return GetDerived()._View(); }

  auto begin() { return View().begin(); }
  auto end() { return View().end(); }

  auto size() { return View().size(); }

  auto operator[](Int i) { return View()[i]; }

  auto operator()(Int iTheta, Int iPhi) {
    auto i = NumberOfLongitudes() * iTheta + iPhi;
    return operator[](i);
  }

  auto Integrate() { return GridReference().Integrate(View()); }

 private:
  auto& GetDerived() { return static_cast<Derived&>(*this); }
  const auto& GetDerived() const { return static_cast<const Derived&>(*this); }
};

// Canonical component class that owns its own data.
template <RealOrComplexFloatingPoint Scalar, typename Grid>
class CanonicalComponent
    : public CanonicalComponentBase<CanonicalComponent<Scalar, Grid>> {
  using Int = std::ptrdiff_t;
  using Vector = FFTWpp::vector<Scalar>;

 public:
  CanonicalComponent() = default;
  CanonicalComponent(const CanonicalComponent&) = default;
  CanonicalComponent(CanonicalComponent&&) = default;

  CanonicalComponent(std::shared_ptr<Grid> grid)
      : _grid{grid},
        _data{Vector(_grid->NumberOfCoLatitudes() * _grid->NumberOfLongitudes(),
                     1)} {}

  template <typename OtherDerived>
  CanonicalComponent(CanonicalComponentBase<OtherDerived>& other)
      : _grid{other.GridPointer()}, _data(other.begin(), other.end()) {}

  CanonicalComponent& operator=(CanonicalComponent&) = default;
  CanonicalComponent& operator=(CanonicalComponent&&) = default;

  template <typename OtherDerived>
  CanonicalComponent& operator=(CanonicalComponentBase<OtherDerived>& other) {
    this->_grid = other.GridPointer();
    this->_data = Vector(other.begin(), other.end());
    return *this;
  }

  auto& operator[](Int i) { return _data[i]; }

  auto& operator()(Int iTheta, Int iPhi) {
    auto i = this->NumberOfLongitudes() * iTheta + iPhi;
    return operator[](i);
  }

 private:
  std::shared_ptr<Grid> _grid;
  Vector _data;

  auto& _GridReference() const { return *_grid; }
  auto _GridPointer() const { return _grid; }
  auto _View() { return boost::sub_range<Vector>(_data); }

  friend class CanonicalComponentBase<CanonicalComponent<Scalar, Grid>>;
};

// Canonical component class that does not own its own data.
template <typename Grid, typename Range>
class CanonicalComponentView
    : public CanonicalComponentBase<CanonicalComponentView<Grid, Range>> {
  using Int = std::ptrdiff_t;

 public:
  CanonicalComponentView(std::shared_ptr<Grid> grid, Range& range)
      : _grid{grid}, _range{range} {}

  template <typename OtherDerived>
  CanonicalComponentView& operator=(
      CanonicalComponentBase<OtherDerived>& other) {
    this->_grid = other.GridPointer();
    std::ranges::copy(other, this->begin());
    return *this;
  }

  auto& operator[](Int i) { return _range[i]; }

  auto& operator()(Int iTheta, Int iPhi) {
    auto i = this->NumberOfLongitudes() * iTheta + iPhi;
    return operator[](i);
  }

 private:
  std::shared_ptr<Grid> _grid;
  Range& _range;

  auto& _GridReference() const { return *_grid; }
  auto _GridPointer() const { return _grid; }
  auto _View() { return boost::sub_range<Range>(_range); }

  friend class CanonicalComponentBase<CanonicalComponentView<Grid, Range>>;
};

// View to a canonical component produced by the action of a unary function.
template <typename Field, typename Function>
class CanonicalComponentUnaryFunction
    : public CanonicalComponentBase<
          CanonicalComponentUnaryFunction<Field, Function>> {
 public:
  CanonicalComponentUnaryFunction(Field& field, Function&& function)
      : _field{field}, _function{function} {}

 private:
  Field& _field;
  Function& _function;

  auto& _GridReference() const { return _field.GridReference(); }
  auto _GridPointer() const { return _field.GridPointer(); }
  auto _View() {
    return _field.View() | boost::adaptors::transformed(_function);
  }

  friend class CanonicalComponentBase<
      CanonicalComponentUnaryFunction<Field, Function>>;
};

// View to a canonical component produced by the action of a binary function.
template <typename Field1, typename Field2, typename Function>
class CanonicalComponentBinaryFunction
    : public CanonicalComponentBase<
          CanonicalComponentBinaryFunction<Field1, Field2, Function>> {
 public:
  CanonicalComponentBinaryFunction(Field1& field1, Field2& field2,
                                   Function&& function)
      : _field1{field1}, _field2{field2}, _function{function} {}

 private:
  Field1& _field1;
  Field2& _field2;
  Function& _function;

  auto& _GridReference() const { return _field1.GridReference(); }

  auto _GridPointer() const { return _field1.GridPointer(); }

  auto _View() {
    return boost::combine(_field1.View(), _field2.View()) |
           boost::adaptors::transformed([this](auto pair) {
             auto x = boost::get<0>(pair);
             auto y = boost::get<1>(pair);
             return _function(x, y);
           });
  }

  friend class CanonicalComponentBase<
      CanonicalComponentBinaryFunction<Field1, Field2, Function>>;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELDS_GUARD_H
