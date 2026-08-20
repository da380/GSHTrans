#ifndef GSH_TRANS_TENSOR_EXPANSION_GUARD_H
#define GSH_TRANS_TENSOR_EXPANSION_GUARD_H

#include <FFTWpp/Core>

#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

#include "../Concepts.h"
#include "../Tensor/MultiIndex.h"
#include "../Tensor/Orbits.h"
#include "../Tensor/TensorField.h"
#include "SpinExpansion.h"

namespace GSHTrans {

// A tensor field in the spectral domain: one buffer, handing out its
// components as spin expansions.
//
// The mirror of TensorField, and deliberately so -- same component addressing,
// same stored set, same grouping by upper index. What it does *not* mirror is
// phase 4's second buffer. A pinned component is a real field, but its
// coefficients are complex numbers in the reduced m >= 0 storage, so the
// spectral side is one complex buffer throughout. What varies is the block
// length: a component's block is sized by its own upper index, and a pinned
// one by the reduced storage, which is why the total is computed rather than
// being a product.
//
// The saving carries across. A real rank-2 tensor stores five components
// rather than nine here too, and each pinned one costs about half what a
// complex one of the same degree would.
template <std::ptrdiff_t _Rank, TensorSymmetry<_Rank> _Symmetry,
          TensorReality _Reality, AngularGrid _Grid>
class TensorExpansion {
 public:
  using Int = std::ptrdiff_t;

  static constexpr Int Rank = _Rank;
  using Symmetry = _Symmetry;
  using Reality = _Reality;
  using GridType = _Grid;
  using Real = typename _Grid::Real;
  using Complex = std::complex<Real>;

  using FieldType = TensorField<Rank, Symmetry, Reality, GridType>;

  static constexpr auto& Orbits = FieldType::Orbits;
  static constexpr auto& ComponentLayout = FieldType::ComponentLayout;
  static constexpr Int Components = FieldType::Components;
  static constexpr Int StoredComponents = FieldType::StoredComponents;
  static constexpr Int ComplexComponents = FieldType::ComplexComponents;
  static constexpr Int RealComponents = FieldType::RealComponents;

  template <Int... Alphas>
  static constexpr bool Represents = FieldType::template Represents<Alphas...>;

  template <Int... Alphas>
  static constexpr bool Writable = FieldType::template Writable<Alphas...>;

  TensorExpansion() = delete;

  TensorExpansion(GridType grid, Int lMax)
      : _grid{std::move(grid)}, _lMax{lMax}, _data(BlockTotal(_grid, lMax)) {
    if (lMax < Rank) {
      throw std::invalid_argument(
          "A rank-" + std::to_string(Rank) +
          " tensor has components at upper index " + std::to_string(Rank) +
          ", so its expansion needs at least that degree");
    }
  }

  const GridType& Grid() const { return _grid; }
  auto MaxDegree() const { return _lMax; }
  auto Size() const { return static_cast<Int>(_data.size()); }
  auto Data() { return std::span<Complex>(_data); }
  auto Data() const { return std::span<const Complex>(_data); }

  // The component with this multi-index, as a spin expansion over the block
  // holding it.
  //
  // Only the *stored* components have blocks. A derived one is determined by
  // its representative, and deriving it in the spectral domain means applying
  // eq:complevel, T^{-N}_{l,-m} = (-1)^m conj(T^N_{lm}) -- an index reversal
  // rather than the pointwise relation the spatial side uses. That is not a
  // view over anything, so it is not offered here: derive in the spatial
  // domain, or ask the representative and apply the relation.
  template <Int... Alphas>
  requires Writable<Alphas...>
  auto Component() {
    constexpr auto flat = FieldType::template FlatOf<Alphas...>;
    constexpr auto slot = FieldType::SlotOfFlat(Orbits.representative[flat]);
    constexpr auto n = ComponentLayout.upperIndexOfSlot[slot];
    constexpr auto real = ComponentLayout.realOfSlot[slot];
    using Value = std::conditional_t<real, RealValued, ComplexValued>;
    return SpinExpansionView<n, GridType, Value>(_grid, _lMax,
                                                 BlockOf(slot));
  }

  template <Int... Alphas>
  requires Writable<Alphas...>
  auto Component() const {
    constexpr auto flat = FieldType::template FlatOf<Alphas...>;
    constexpr auto slot = FieldType::SlotOfFlat(Orbits.representative[flat]);
    constexpr auto n = ComponentLayout.upperIndexOfSlot[slot];
    constexpr auto real = ComponentLayout.realOfSlot[slot];
    using Value = std::conditional_t<real, RealValued, ComplexValued>;
    return ConstSpinExpansionView<n, GridType, Value>(_grid, _lMax,
                                                      BlockOf(slot));
  }

  // The block a stored component occupies, by slot.
  std::span<Complex> BlockOf(Int slot) {
    const auto [offset, size] = Block(_grid, _lMax, slot);
    return Data().subspan(offset, size);
  }

  std::span<const Complex> BlockOf(Int slot) const {
    const auto [offset, size] = Block(_grid, _lMax, slot);
    return Data().subspan(offset, size);
  }

 private:
  GridType _grid;
  Int _lMax;
  FFTWpp::vector<Complex> _data;

  // Where a slot's block starts and how long it is. The blocks follow the
  // component order, which is the order the transform writes them in.
  static std::pair<std::size_t, std::size_t> Block(const GridType& grid,
                                                   Int lMax, Int slot) {
    auto offset = std::size_t{0};
    for (auto i = Int{0}; i < slot; i++) offset += BlockSize(grid, lMax, i);
    return {offset, BlockSize(grid, lMax, slot)};
  }

  static std::size_t BlockSize(const GridType& grid, Int lMax, Int slot) {
    if (ComponentLayout.realOfSlot[slot]) {
      return static_cast<std::size_t>(grid.RealCoefficientSize(lMax));
    }
    return static_cast<std::size_t>(
        grid.CoefficientSize(lMax, ComponentLayout.upperIndexOfSlot[slot]));
  }

  static std::size_t BlockTotal(const GridType& grid, Int lMax) {
    auto total = std::size_t{0};
    for (auto slot = Int{0}; slot < StoredComponents; slot++) {
      total += BlockSize(grid, lMax, slot);
    }
    return total;
  }
};

// Between the two domains, arranging the batched transform the tensor field
// already knows how to make.
template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid, TensorLayout Layout>
auto Expand(const TensorField<Rank, Symmetry, Reality, Grid, Layout>& tensor,
            std::ptrdiff_t lMax,
            Execution policy = Execution::Sequential()) {
  auto expansion =
      TensorExpansion<Rank, Symmetry, Reality, Grid>(tensor.Grid(), lMax);
  tensor.ForwardTransformation(lMax, expansion.Data(), policy);
  return expansion;
}

template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid,
          TensorLayout Layout = ComponentMajor>
auto Evaluate(const TensorExpansion<Rank, Symmetry, Reality, Grid>& expansion,
              Execution policy = Execution::Sequential()) {
  auto tensor =
      TensorField<Rank, Symmetry, Reality, Grid, Layout>(expansion.Grid());
  tensor.InverseTransformation(expansion.MaxDegree(), expansion.Data(), policy);
  return tensor;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_TENSOR_EXPANSION_GUARD_H
