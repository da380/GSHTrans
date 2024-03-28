#ifndef GSH_TRANS_TENSOR_FIELD_GUARD_H
#define GSH_TRANS_TENSOR_FIELD_GUARD_H

#include <array>
#include <map>
#include <ranges>
#include <tuple>

#include "TensorIndices.h"

namespace GSHTrans {

template <std::size_t _Rank, RealOrComplexValued _Value, typename _Grid,
          std::ranges::view _View>
class TensorField : public std::ranges::view_interface<
                        TensorField<_Rank, _Value, _Grid, _View>> {
  using Int = std::ptrdiff_t;

 public:
  using Real = typename _Grid::Real;
  using Complex = typename _Grid::Complex;
  using Scalar =
      std::conditional_t<std::same_as<_Value, RealValued>, Real, Complex>;
  using Grid = _Grid;

  TensorField() = default;

  TensorField(_Grid grid, _View view) : _grid{grid}, _view{view} {
    assert(_view.size() == TensorIndices<_Rank>().size() * _grid.FieldSize());
    Int count = 0;
    for (auto index : TensorIndex<_Rank>()) {
      _map[index] = count++;
    }
  }

  auto operator[](TensorIndex<_Rank>& index) { return 0; }

  auto operator[](TensorIndex<_Rank>&& index) { operator[](index); }

  template <typename... I>
  requires(sizeof...(I) == _Rank) && (std::convertible_to<I, Int> && ...)
  auto operator[](I... alpha) {
    return operator[](TensorIndex(alpha...));
  }

 private:
  _Grid _grid;
  _View _view;
  std::map<TensorIndex<_Rank>, Int> _map;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_TENSOR_FIELD_GUARD_H