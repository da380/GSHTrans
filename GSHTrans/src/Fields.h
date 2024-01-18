#ifndef GSH_TRANS_FIELDS_GUARD_H
#define GSH_TRANS_FIELDS_GUARD_H

#include <FFTWpp/All>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

// Base class for scalar fields.
template <typename Derived>
class ScalarFieldBase {
  using Int = std::ptrdiff_t;

 public:
  auto& GridReference() const { return GetDerived()._GridReference(); }
  auto GridPointer() const { return GetDerived()._GridPointer(); }

  Int MaxDegree() const { return GridReference().MaxDegree(); }
  Int MaxUpperIndex() const { return GridReference().MaxUpperIndex(); }

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

// Scalar field class that owns its own data.
template <typename Grid>
class ScalarField : public ScalarFieldBase<ScalarField<Grid>> {
 public:
  using real_type = Grid::real_type;
  using scalar_type = Grid::scalar_type;

 private:
  using Int = std::ptrdiff_t;
  using Vector = FFTWpp::vector<scalar_type>;

 public:
  ScalarField() = default;
  ScalarField(const ScalarField&) = default;
  ScalarField(ScalarField&&) = default;

  ScalarField(std::shared_ptr<Grid> grid)
      : _grid{grid},
        _data{Vector(_grid->NumberOfCoLatitudes() * _grid->NumberOfLongitudes(),
                     1)} {}

  ScalarField(Int lMax, FFTWpp::PlanFlag flag = FFTWpp::Measure)
      : _grid{std::make_shared<Grid>(lMax, 0, flag)},
        _data{Vector(_grid->NumberOfCoLatitudes() * _grid->NumberOfLongitudes(),
                     1)} {}

  template <typename OtherDerived>
  ScalarField(ScalarFieldBase<OtherDerived>& other)
      : _grid{other.GridPointer()}, _data(other.begin(), other.end()) {}

  ScalarField& operator=(ScalarField&) = default;
  ScalarField& operator=(ScalarField&&) = default;

  template <typename OtherDerived>
  ScalarField& operator=(ScalarFieldBase<OtherDerived>& other) {
    this->_grid = other.GridPointer();
    this->_data = Vector(other.begin(), other.end());
    return *this;
  }

  auto& operator[](Int i) { return _data[i]; }

  auto& operator()(Int iTheta, Int iPhi) {
    auto i = _grid->NumberOfLongitudes() * iTheta + iPhi;
    return operator[](i);
  }

 private:
  std::shared_ptr<Grid> _grid;
  Vector _data;

  auto& _GridReference() const { return *_grid; }
  auto _GridPointer() const { return _grid; }
  auto _View() { return boost::sub_range<Vector>(_data); }

  friend class ScalarFieldBase<ScalarField<Grid>>;
};

// Scalar field class that does not own its own data.
template <typename Grid, typename Range>
class ScalarFieldView : public ScalarFieldBase<ScalarFieldView<Grid, Range>> {
  using Int = std::ptrdiff_t;

 public:
  ScalarFieldView(std::shared_ptr<Grid> grid, Range& range)
      : _grid{grid}, _range{range} {}

  auto& operator[](Int i) { return _range[i]; }

  auto& operator()(Int iTheta, Int iPhi) {
    auto i = _grid->NumberOfLongitudes() * iTheta + iPhi;
    return operator[](i);
  }

 private:
  std::shared_ptr<Grid> _grid;
  Range& _range;

  auto& _GridReference() const { return *_grid; }
  auto _GridPointer() const { return _grid; }
  auto _View() { return boost::sub_range<Range>(_range); }

  friend class ScalarFieldBase<ScalarFieldView<Grid, Range>>;
};

template <typename Field, typename Function>
class ScalarFieldUnaryFunction
    : public ScalarFieldBase<ScalarFieldUnaryFunction<Field, Function>> {
 public:
  ScalarFieldUnaryFunction(Field& field, Function&& function)
      : _field{field}, _function{function} {}

 private:
  Field& _field;
  Function& _function;

  auto& _GridReference() const { return _field.GridReference(); }
  auto _GridPointer() const { return _field.GridPointer(); }
  auto _View() {
    return _field.View() | boost::adaptors::transformed(_function);
  }

  friend class ScalarFieldBase<ScalarFieldUnaryFunction<Field, Function>>;
};

template <typename Field1, typename Field2, typename Function>
class ScalarFieldBinaryFunction
    : public ScalarFieldBase<
          ScalarFieldBinaryFunction<Field1, Field2, Function>> {
 public:
  ScalarFieldBinaryFunction(Field1& field1, Field2& field2, Function&& function)
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

  friend class ScalarFieldBase<
      ScalarFieldBinaryFunction<Field1, Field2, Function>>;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELDS_GUARD_H
