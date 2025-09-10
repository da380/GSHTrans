#ifndef GSH_TRANS_SCALAR_FIELD_VIEW_GUARD_H
#define GSH_TRANS_SCALAR_FIELD_VIEW_GUARD_H

#include <FFTWpp/Core>
#include <concepts>
#include <iostream>
#include <limits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "ScalarFieldBase.h"

namespace GSHTrans {

// Forward declare class.
template <typename _Grid, std::ranges::view _View>
requires std::derived_from<_Grid, GridBase<_Grid>> &&
         std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                      typename _Grid::Real>
class ScalarFieldView;

// Set traits.
namespace Internal {

template <typename _Grid, std::ranges::view _View>
struct Traits<ScalarFieldView<_Grid, _View>> {
  using Int = std::ptrdiff_t;
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar = std::ranges::range_value_t<_View>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;
  using Writeable = std::conditional_t<
      std::ranges::output_range<_View, std::ranges::range_value_t<_View>>,
      std::true_type, std::false_type>;
};

}  // namespace Internal

template <typename _Grid, std::ranges::view _View>
requires std::derived_from<_Grid, GridBase<_Grid>> &&
         std::same_as<RemoveComplex<std::ranges::range_value_t<_View>>,
                      typename _Grid::Real>
class ScalarFieldView : public ScalarFieldBase<ScalarFieldView<_Grid, _View>> {
 public:
  using Value = typename Internal::Traits<ScalarFieldView<_Grid, _View>>::Value;
  using Int = typename Internal::Traits<ScalarFieldView<_Grid, _View>>::Int;
  using Real = typename Internal::Traits<ScalarFieldView<_Grid, _View>>::Real;
  using Complex =
      typename Internal::Traits<ScalarFieldView<_Grid, _View>>::Complex;
  using Scalar =
      typename Internal::Traits<ScalarFieldView<_Grid, _View>>::Scalar;
  using Writeable =
      typename Internal::Traits<ScalarFieldView<_Grid, _View>>::Writeable;

  // Return the grid.
  auto& Grid() const { return _grid; }

  // Read access to the data.
  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

  // Write access to the data.
  auto& operator[](Int iTheta, Int iPhi)
  requires Writeable::value
  {
    this->CheckPointIndices(iTheta, iPhi);
    return _data[Index(iTheta, iPhi)];
  }

  // Constructors.
  ScalarFieldView() = default;

  // Construct from grid initialising values to zero.
  ScalarFieldView(_Grid& grid, _View data) : _grid{grid}, _data{data} {}

  // Default copy and move constructors.
  ScalarFieldView(const ScalarFieldView&) = default;
  ScalarFieldView(ScalarFieldView&&) = default;

  // Default copy and move assignment.
  ScalarFieldView& operator=(const ScalarFieldView&) = default;
  ScalarFieldView& operator=(ScalarFieldView&&) = default;

  // Use assignment defined in base class.
  using ScalarFieldBase<ScalarFieldView<_Grid, _View>>::operator=;

  // Return view to the data.
  auto Data() { return _data; }

 private:
  _Grid& _grid;
  _View _data;

  auto Index(Int iTheta, int iPhi) const {
    return iTheta * this->NumberOfLongitudes() + iPhi;
  }
};

// Deduction guide to allow construction from a range.
template <typename Grid, std::ranges::viewable_range R>
ScalarFieldView(Grid&,
                R&&) -> ScalarFieldView<Grid, std::ranges::views::all_t<R>>;

}  // namespace GSHTrans

#endif  // #ifndef GSH_TRANS_SCALAR_FIELD_GUARD_H
