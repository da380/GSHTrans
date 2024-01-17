#ifndef GSH_TRANS_FIELDS_GUARD_H
#define GSH_TRANS_FIELDS_GUARD_H

#include <FFTWpp/All>
#include <iostream>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

class UnaryPlus {
 public:
  auto operator()(auto&& x) { return x; }
  auto operator()(auto&& x) const { return x; }
};

class UnaryMinus {
 public:
  auto operator()(auto&& x) { return -x; }
  auto operator()(auto&& x) const { return -x; }
};

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

  // Return iterators.
  auto begin() { return GetDerived()._Range().begin(); }
  auto end() { return GetDerived()._Range().end(); }

  auto cbegin() const { return GetDerived()._cbegin(); }
  auto cend() const { return GetDerived()._cend(); }

  // Return size.
  auto size() const { return GetDerived()._Range().size(); }

  // Index operators
  auto operator[](Int i) { return GetDerived()._Range()[i]; }
  auto operator[](Int i) const { return GetDerived()._Index(i); }

  // Return value of the field at the (iTheta,iPhi)th location.
  auto operator()(Int iTheta, Int iPhi) const {
    auto i = NumberOfLongitudes() * iTheta + iPhi;
    return operator[](i);
  }

  // Return the integral of the field.
  auto Integrate() const {
    auto start = GetDerived()._cbegin();
    auto finish = GetDerived()._cend();
    auto range = std::ranges::subrange(start, finish);
    return GridReference().Integrate(range);
  }

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
  ScalarField(const ScalarFieldBase<OtherDerived>& other)
      : _grid{other.GridPointer()}, _data(other.cbegin(), other.cend()) {}

  // Set default assigment operators.
  ScalarField& operator=(const ScalarField&) = default;
  ScalarField& operator=(ScalarField&&) = default;

  // Copy assignment from other derived class.
  template <typename OtherDerived>
  ScalarField& operator=(const ScalarFieldBase<OtherDerived>& other) {
    this->_grid = other.GridPointer();
    this->_data = Vector(other.cbegin(), other.cend());
    return *this;
  }

  // Define the required public methods.
  auto& _GridReference() const { return *_grid; }
  auto _GridPointer() const { return _grid; }

  auto _Range() { return std::ranges::views::all(_data); }

  auto _cbegin() const { return _data.cbegin(); }
  auto _cend() const { return _data.cend(); }

  // Define additional public methods.

 private:
  std::shared_ptr<Grid> _grid;
  Vector _data;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_FIELDS_GUARD_H
