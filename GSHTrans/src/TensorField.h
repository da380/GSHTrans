#ifndef GSH_TRANS_TENSOR_FIELD_GUARD_H
#define GSH_TRANS_TENSOR_FIELD_GUARD_H

#include <FFTWpp/Memory>
#include <iostream>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

template <typename Derived>
class ScalarFieldBase {
  using Int = std::ptrdiff_t;

 public:
  auto operator()(Int iTheta, Int iPhi) const {
    const auto& derived = GetDerived();
    auto i = derived.NumberOfLongitudes() * iTheta + iPhi;
    return derived[i];
  }

 private:
  auto& GetDerived() { return static_cast<Derived&>(*this); }
  const auto& GetDerived() const { return static_cast<const Derived&>(*this); }
};

template <typename GSHGrid>
class ScalarField : public ScalarFieldBase<ScalarField<GSHGrid>> {
  using Real = GSHGrid::real_type;
  using Scalar = GSHGrid::scalar_type;
  using Int = std::ptrdiff_t;
  using Vector = FFTWpp::vector<Scalar>;

 public:
  using real_type = GSHGrid::real_type;
  using scalar_type = GSHGrid::scalar_type;

  ScalarField(std::shared_ptr<GSHGrid> grid)
      : _grid{grid},
        _data{Vector(_grid->NumberOfCoLatitudes() * _grid->NumberOfLongitudes(),
                     1)} {}

  auto operator[](Int i) { return _data[i]; }

  auto operator[](Int i) const { return _data[i]; }

  auto NumberOfLongitudes() const { return _grid->NumberOfLongitudes(); }

  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

 private:
  std::shared_ptr<GSHGrid> _grid;
  Vector _data;
};

template <typename Field>
class ScalarFieldView : public ScalarFieldBase<ScalarFieldView<Field>> {
  using Int = std::ptrdiff_t;

 public:
  using scalar_type = Field::scalar_type;

  ScalarFieldView(Field& field, scalar_type scale = 1)
      : _field{field}, _scale{scale} {}

  auto operator[](Int i) { return _field[i] * _scale; }
  auto operator[](Int i) const { return _field[i] * _scale; }

  auto NumberOfLongitudes() const { return _field.NumberOfLongitudes(); }

 private:
  scalar_type _scale;
  Field& _field;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_TENSOR_FIELD_GUARD_H
