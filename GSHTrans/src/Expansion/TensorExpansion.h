#ifndef GSH_TRANS_TENSOR_EXPANSION_GUARD_H
#define GSH_TRANS_TENSOR_EXPANSION_GUARD_H

#include <FFTWpp/Core>

#include <cmath>
#include <complex>
#include <type_traits>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

#include "../Concepts.h"
#include "../Utility.h"
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

  // The coefficient of *any* representable component at (l, m).
  //
  // The block accessors above reach only stored components, because a block
  // is what a view can be taken over. A derived component has no block: on the
  // spectral side deriving one is eq:complevel,
  //
  //   T^{-N}_{l,-m} = (-1)^m conj(T^{N}_{lm}),
  //
  // which reverses the order index rather than acting at fixed (l, m), so it
  // is a computation and not a view. This is that computation, and with it
  // every component of the tensor is readable in either domain.
  //
  // Four things have to be resolved, and the orbit table says which applies:
  //
  //   stored                  read it
  //   permutation relative    the representative's, with a sign
  //   reality relative        the representative's at -m, conjugated
  //   vanishing orbit         zero
  //
  // and two more come from how the block itself is stored: a pinned component
  // is a real field, so its block holds only m >= 0 and the negative orders
  // follow from f_{l,-m} = (-1)^m conj(f_{lm}); and one pinned as imaginary
  // holds the real field whose i-multiple the component is.
  //
  // Degrees below the component's own |N| return zero rather than reading off
  // the end: a component at upper index N has no content there, which is not
  // a missing value but an absent one. So does a component whose orbit
  // vanishes -- unlike the spatial accessor, which refuses those because a
  // node must have a type and there is nothing to give it. Here the answer is
  // a value, and zero is the right one; the surface gradient reads shifted
  // components that may vanish and would otherwise have to special-case them.
  template <Int... Alphas>
  requires(sizeof...(Alphas) == Rank)
  Complex Coefficient(Int l, Int m) const {
    constexpr auto flat = FieldType::template FlatOf<Alphas...>;
    constexpr auto n = FieldType::template UpperIndexOf<Alphas...>;
    constexpr auto constraint = Orbits.constraint[flat];

    if constexpr (constraint == ComponentConstraint::Zero) {
      return Complex{};
    } else {
      if (l < (n < 0 ? -n : n) || l > _lMax || m < -l || m > l) {
        return Complex{};
      }

      constexpr auto rep = Orbits.representative[flat];
      constexpr auto slot = FieldType::SlotOfFlat(rep);
      constexpr auto repN = ComponentLayout.upperIndexOfSlot[slot];
      constexpr auto sign = static_cast<Real>(Orbits.sign[flat]);
      constexpr auto conjugated = Orbits.conjugate[flat];
      constexpr auto real = ComponentLayout.realOfSlot[slot];

      // The representative's own coefficient, at whichever order this term
      // needs, allowing for a real block's reduced storage.
      const auto stored = [&](Int order) {
        auto block = ConstSpinExpansionView<repN, GridType,
                                            std::conditional_t<real, RealValued,
                                                               ComplexValued>>(
            _grid, _lMax, BlockOf(slot));
        if constexpr (real) {
          if (order < 0) {
            return static_cast<Real>(MinusOneToPower(order)) *
                   std::conj(block[l, -order]);
          }
        }
        return Complex{block[l, order]};
      };

      // A pinned-imaginary component is i times the real field stored for it.
      constexpr auto turn =
          constraint == ComponentConstraint::Imaginary ? Complex{0, 1}
                                                       : Complex{1, 0};

      if constexpr (conjugated) {
        return sign * turn * static_cast<Real>(MinusOneToPower(m + repN)) *
               std::conj(stored(-m));
      } else {
        return sign * turn * stored(m);
      }
    }
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
