#ifndef GSH_TRANS_FIELDS_GUARD_H
#define GSH_TRANS_FIELDS_GUARD_H

#include <FFTWpp/All>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

template <typename Derived>
class ScalarFieldBase {
  using Int = std::ptrdiff_t;

 public:
  // Return grid information.
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

  // Return range.
  auto DataView() { return GetDerived()._DataView(); }

  // Return iterators.
  auto begin() { return DataView().begin(); }
  auto end() { return DataView().end(); }

  // Return size.
  auto size() { return DataView().size(); }

  // Index operators
  auto operator[](Int i) { return DataView()[i]; }

  // Return value of the field at the (iTheta,iPhi)th location.
  auto operator()(Int iTheta, Int iPhi) {
    auto i = NumberOfLongitudes() * iTheta + iPhi;
    return operator[](i);
  }

  // Return the integral of the field.
  auto Integrate() { return GridReference().Integrate(DataView()); }

 private:
  // Return references to the instance of the derived class.
  auto& GetDerived() { return static_cast<Derived&>(*this); }
  const auto& GetDerived() const { return static_cast<const Derived&>(*this); }
};

template <typename Grid>
class ScalarField : public ScalarFieldBase<ScalarField<Grid>> {
 public:
  using real_type = Grid::real_type;
  using scalar_type = Grid::scalar_type;

 private:
  using Int = std::ptrdiff_t;
  using Vector = FFTWpp::vector<scalar_type>;

 public:
  // Set default constructors.
  ScalarField() = default;
  ScalarField(const ScalarField&) = default;
  ScalarField(ScalarField&&) = default;

  // Constructor given a grid initialising values to zero.
  ScalarField(std::shared_ptr<Grid> grid)
      : _grid{grid},
        _data{Vector(_grid->NumberOfCoLatitudes() * _grid->NumberOfLongitudes(),
                     1)} {}

  // Constructor given grid parameters.
  ScalarField(Int lMax, FFTWpp::PlanFlag flag = FFTWpp::Measure)
      : _grid{std::make_shared<Grid>(lMax, 0, flag)},
        _data{Vector(_grid->NumberOfCoLatitudes() * _grid->NumberOfLongitudes(),
                     1)} {}

  // Copy constructor from other derived class.
  template <typename OtherDerived>
  ScalarField(ScalarFieldBase<OtherDerived>& other)
      : _grid{other.GridPointer()}, _data(other.begin(), other.end()) {}

  // Set default assigment operators.
  ScalarField& operator=(ScalarField&) = default;
  ScalarField& operator=(ScalarField&&) = default;

  // Copy assignment from other derived class.
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
  auto _DataView() { return std::ranges::views::all(_data); }

  friend class ScalarFieldBase<ScalarField<Grid>>;
};

template <typename View>
class ScalarFieldView : public ScalarFieldBase<ScalarFieldView<View>> {
  using Int = std::ptrdiff_t;

 public:
  ScalarFieldView(View& view) : _view{view} {}

  auto& operator[](Int i) { return _view[i]; }

  auto& operator()(Int iTheta, Int iPhi) {
    auto i = _GridReference().NumberOfLongitudes() * iTheta + iPhi;
    return operator[](i);
  }

 private:
  View& _view;

  auto& _GridReference() const { return _view.GridReference(); }
  auto _GridPointer() const { return _view.GridPointer(); }
  auto _DataView() { return std::ranges::views::all(_view); }

  friend class ScalarFieldBase<ScalarFieldView<View>>;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELDS_GUARD_H
