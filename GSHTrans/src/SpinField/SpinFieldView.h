#ifndef GSH_TRANS_SPIN_FIELD_VIEW_GUARD_H
#define GSH_TRANS_SPIN_FIELD_VIEW_GUARD_H

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "../Concepts.h"
#include "SpinWeighted.h"

namespace GSHTrans {

// A non-owning field: someone else's storage, plus a grid handle.
//
// This is what makes the tensor layer possible. A TensorField owns one
// contiguous buffer and hands out components as views at the correct upper
// index, rather than holding a tuple of separately allocated fields; the
// radial slices of a layered field are the same idea. Both need a view to be
// admissible wherever an owning field is, which is why operator[] returns by
// value on every node and why the grid is a handle rather than a reference.
//
// Views are terminals: they name storage, so an expression may hold an lvalue
// one by reference.
template <std::ptrdiff_t _N, AngularGrid _Grid,
          RealOrComplexValued _Value = ComplexValued,
          typename _Element = ScalarFor<typename _Grid::Real, _Value>>
class SpinFieldView {
 public:
  using Int = std::ptrdiff_t;

  static constexpr Int UpperIndex = _N;
  using Value = _Value;
  using GridType = _Grid;
  using Real = typename _Grid::Real;
  using Complex = std::complex<Real>;
  using Scalar = std::remove_const_t<_Element>;

  static_assert(std::same_as<Scalar, ScalarFor<Real, Value>>);
  static_assert(std::same_as<Value, ComplexValued> or UpperIndex == 0,
                "A spin-weighted field can be real-valued only at upper index "
                "zero");
  static_assert(
      not std::same_as<typename _Grid::NRange, NonNegative> or UpperIndex >= 0,
      "This grid stores only non-negative upper indices");

  SpinFieldView() = delete;

  SpinFieldView(GridType grid, std::span<_Element> data)
      : _grid{std::move(grid)}, _data{data} {
    if (!std::ranges::contains(_grid.UpperIndices(), UpperIndex)) {
      throw std::invalid_argument("This grid does not carry upper index " +
                                  std::to_string(UpperIndex));
    }
    if (_data.size() != static_cast<std::size_t>(_grid.FieldSize())) {
      throw std::invalid_argument(
          "A view over " + std::to_string(_data.size()) +
          " values does not cover this grid's " +
          std::to_string(_grid.FieldSize()) + " points");
    }
  }

  const GridType& Grid() const { return _grid; }

  Scalar operator[](Int iTheta, Int iPhi) const {
    return _data[FlatIndex(iTheta, iPhi)];
  }

  // Present only on a view over mutable storage.
  _Element& operator[](Int iTheta, Int iPhi)
  requires(not std::is_const_v<_Element>)
  {
    return _data[FlatIndex(iTheta, iPhi)];
  }

  template <typename S>
  requires std::convertible_to<Scalar, S>
  void EvaluateInto(std::span<S> target) const {
    if (target.size() != _data.size()) {
      throw std::invalid_argument(
          "Evaluation target has size " + std::to_string(target.size()) +
          ", but this field has " + std::to_string(_data.size()) + " points");
    }
    std::ranges::copy(_data, target.begin());
  }

  auto Size() const { return static_cast<Int>(_data.size()); }
  auto Data() const { return _data; }

 private:
  GridType _grid;
  std::span<_Element> _data;

  Int FlatIndex(Int iTheta, Int iPhi) const {
    const auto nPhi = static_cast<Int>(_grid.NumberOfLongitudes());
    assert(iTheta >= 0 &&
           iTheta < static_cast<Int>(_grid.NumberOfCoLatitudes()));
    assert(iPhi >= 0 && iPhi < nPhi);
    return iTheta * nPhi + iPhi;
  }
};

// A view over storage nobody may write through.
template <std::ptrdiff_t N, AngularGrid Grid,
          RealOrComplexValued Value = ComplexValued>
using ConstSpinFieldView =
    SpinFieldView<N, Grid, Value, const ScalarFor<typename Grid::Real, Value>>;

template <std::ptrdiff_t N, AngularGrid Grid, RealOrComplexValued Value,
          typename Element>
struct IsTerminalTrait<SpinFieldView<N, Grid, Value, Element>>
    : std::true_type {};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SPIN_FIELD_VIEW_GUARD_H
