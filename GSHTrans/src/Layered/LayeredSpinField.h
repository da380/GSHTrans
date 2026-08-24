#ifndef GSH_TRANS_LAYERED_SPIN_FIELD_GUARD_H
#define GSH_TRANS_LAYERED_SPIN_FIELD_GUARD_H

#include <FFTWpp/Core>

#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

#include "../Concepts.h"
#include "../Policies.h"
#include "../Indexing.h"
#include "../SpinField/SpinField.h"
#include "../SpinField/SpinFieldOverloads.h"
#include "../SpinField/SpinFieldView.h"
#include "../SpinField/SpinWeighted.h"
#include "RadialGrid.h"

namespace GSHTrans {

/// A field on the product of a radial grid with an angular one: nR angular
/// fields laid end to end, radius outermost, each slice contiguous in the
/// canonical angular order.
///
/// Two-dimensional is the primitive and three-dimensional is a stack. The
/// angular field is never wrapped and a slice is not a new kind of object: it
/// is a SpinFieldView, an ordinary spin-weighted node, so the index algebra, the
/// evaluation and the aliasing theorem all lift unchanged and there is no
/// second expression system to keep consistent with the first. That is what
/// views are admissible everywhere an owning field is.
///
/// The layout is radius-major because that is what the applications this exists
/// for already use, and because it makes the radial axis a batch: the stack of
/// one field's slices is `nR` blocks of `FieldSize`, which the transform's
/// (count, stride, dist) descriptor reads directly.
///
/// A one-slice stack is *not* an angular field. The idiom stays available to
/// application code, but the two are distinct types here, so that a function
/// taking one cannot silently be handed the other.
template <std::ptrdiff_t _N, AngularGrid _Grid,
          RealOrComplexValued _Value = ComplexValued>
class LayeredSpinField {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.

  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = _N;
  using Value = _Value;  ///< Whether the samples are real-valued or complex.
  using GridType = _Grid;  ///< The angular grid this is defined on.
  using Real = typename _Grid::Real;  ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  using Scalar = ScalarFor<Real, Value>;  ///< The value type: Real when real-valued, Complex otherwise.
  using RadialGridType = RadialGrid<Real>;
  using SliceType = SpinFieldView<_N, _Grid, _Value>;
  using ConstSliceType = ConstSpinFieldView<_N, _Grid, _Value>;

  static_assert(std::same_as<Value, ComplexValued> or UpperIndex == 0,
                "A spin-weighted field can be real-valued only at upper index "
                "zero");

  LayeredSpinField() = delete;

  LayeredSpinField(RadialGridType radialGrid, GridType grid)
      : _radialGrid{std::move(radialGrid)},
        _grid{std::move(grid)},
        _data(static_cast<std::size_t>(_radialGrid.NumberOfRadii()) *
              static_cast<std::size_t>(_grid.FieldSize())) {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _grid; }
  /** @brief The radial grid this is defined on. */
  const RadialGridType& Radial() const { return _radialGrid; }

  /** @brief How many radii the stack holds. */
  auto NumberOfRadii() const { return _radialGrid.NumberOfRadii(); }
  /** @brief Indices of the stored radii. */
  auto RadiusIndices() const { return _radialGrid.RadiusIndices(); }
  /** @brief How many samples one angular field holds. */
  auto FieldSize() const { return static_cast<Int>(_grid.FieldSize()); }
  /** @brief How many elements are stored. */
  auto Size() const { return static_cast<Int>(_data.size()); }

  // The uniform names a radial operator sees. A field's slice is a set of
  // angular points and an expansion's is a set of coefficients, but the radial
  // axis does not care which: it runs over `NumberOfRadii()` values `SliceSize()`
  // apart, and that is all `ApplyRadially` needs to know about either.
  /** @brief How many elements one radial slice holds. */
  auto SliceSize() const { return FieldSize(); }
  auto SameShape() const { return LayeredSpinField(_radialGrid, _grid); }

  /// The same field on a different set of radii, which is what resampling
  /// needs and what SameShape cannot give: a radial operator maps a stack to
  /// one of the same shape, and changing nR is precisely not that.
  auto SameShapeOn(RadialGridType radialGrid) const {
    return LayeredSpinField(std::move(radialGrid), _grid);
  }

  /** @brief The underlying buffer. */
  auto Data() { return std::span<Scalar>(_data); }
  /** @brief The underlying buffer. */
  auto Data() const { return std::span<const Scalar>(_data); }

  /// One angular field, as a view over this stack's storage. Writing through it
  /// writes the stack.
  auto Slice(Int i) {
    return SliceType(_grid, Data().subspan(Offset(i), SliceExtent()));
  }

  auto Slice(Int i) const {
    return ConstSliceType(_grid, Data().subspan(Offset(i), SliceExtent()));
  }

  /// The whole stack as the transform's batch: nR fields, each contiguous,
  /// FieldSize apart.
  auto Batch() const {
    return GSHTrans::Batch::Contiguous(NumberOfRadii(), FieldSize());
  }

 private:
  RadialGridType _radialGrid;
  GridType _grid;
  FFTWpp::vector<Scalar> _data;

  std::size_t SliceExtent() const {
    return static_cast<std::size_t>(_grid.FieldSize());
  }

  std::size_t Offset(Int i) const {
    if (i < 0 || i >= NumberOfRadii()) {
      throw std::invalid_argument(
          "Radius index " + std::to_string(i) + " is outside the " +
          std::to_string(NumberOfRadii()) + " radii of this grid");
    }
    return static_cast<std::size_t>(i) * SliceExtent();
  }
};

//--------------------------------------------------------------------------//
//                       Between the two representations                     //
//--------------------------------------------------------------------------//

/// A stack of coefficients: one expansion per radius, radius-major, so that the
/// layout matches the field it came from and the angular transform sees a
/// contiguous batch on both sides.
///
/// This is the [r][(l,m)] of the two layouts the plan names. The other,
/// [(l,m)][r], is what a radial solve at fixed degree and order wants, and the
/// repack between them is a separate step.
template <std::ptrdiff_t _N, AngularGrid _Grid,
          RealOrComplexValued _Value = ComplexValued>
class LayeredSpinExpansion {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.

  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = _N;
  using Value = _Value;  ///< Whether the samples are real-valued or complex.
  using GridType = _Grid;  ///< The angular grid this is defined on.
  using Real = typename _Grid::Real;  ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  using RadialGridType = RadialGrid<Real>;
  using MRange =
      std::conditional_t<std::same_as<_Value, RealValued>, NonNegative, All>;

  LayeredSpinExpansion() = delete;

  LayeredSpinExpansion(RadialGridType radialGrid, GridType grid, Int lMax)
      : _radialGrid{std::move(radialGrid)},
        _grid{std::move(grid)},
        _indices{Checked(lMax), lMax, UpperIndex},
        _data(static_cast<std::size_t>(_radialGrid.NumberOfRadii()) *
              static_cast<std::size_t>(_indices.Size())) {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _grid; }
  /** @brief The radial grid this is defined on. */
  const RadialGridType& Radial() const { return _radialGrid; }

  /** @brief How many radii the stack holds. */
  auto NumberOfRadii() const { return _radialGrid.NumberOfRadii(); }
  /** @brief Indices of the stored radii. */
  auto RadiusIndices() const { return _radialGrid.RadiusIndices(); }
  /** @brief The largest degree stored. */
  auto MaxDegree() const { return _indices.MaxDegree(); }
  /** @brief The smallest degree stored. */
  auto MinDegree() const { return _indices.MinDegree(); }
  /** @brief Every degree stored. */
  auto Degrees() const { return _indices.Degrees(); }
  /** @brief Every order stored at degree @p l. */
  auto Orders(Int l) const {
    return GSHSubIndices<MRange>(l, MaxDegree()).Orders();
  }
  /** @brief How many coefficients one block holds. */
  auto CoefficientSize() const { return static_cast<Int>(_indices.Size()); }
  /** @brief How many elements are stored. */
  auto Size() const { return static_cast<Int>(_data.size()); }

  /** @brief How many elements one radial slice holds. */
  auto SliceSize() const { return CoefficientSize(); }

  /// Where a degree and order sit within one radius's block. The radial-major
  /// repack needs this: once the layout is [(l, m)][r] a line is addressed by
  /// its position in the block and no longer by (l, m) directly.
  auto CoefficientIndex(Int l, Int m) const {
    return static_cast<Int>(_indices.Index(l, m));
  }

  auto SameShape() const {
    return LayeredSpinExpansion(_radialGrid, _grid, MaxDegree());
  }

  auto SameShapeOn(RadialGridType radialGrid) const {
    return LayeredSpinExpansion(std::move(radialGrid), _grid, MaxDegree());
  }

  /** @brief The underlying buffer. */
  auto Data() { return std::span<Complex>(_data); }
  /** @brief The underlying buffer. */
  auto Data() const { return std::span<const Complex>(_data); }

  /// The coefficient at one radius.
  Complex operator[](Int i, Int l, Int m) const {
    return _data[Offset(i) + static_cast<std::size_t>(_indices.Index(l, m))];
  }
  Complex& operator[](Int i, Int l, Int m) {
    return _data[Offset(i) + static_cast<std::size_t>(_indices.Index(l, m))];
  }

  auto Batch() const {
    return GSHTrans::Batch::Contiguous(NumberOfRadii(), CoefficientSize());
  }

 private:
  RadialGridType _radialGrid;
  GridType _grid;
  GSHIndices<MRange> _indices;
  FFTWpp::vector<Complex> _data;

  static Int Checked(Int lMax) {
    if (lMax < (UpperIndex < 0 ? -UpperIndex : UpperIndex)) {
      throw std::invalid_argument(
          "An expansion's degree cannot be below its upper index");
    }
    return lMax;
  }

  std::size_t Offset(Int i) const {
    return static_cast<std::size_t>(i) *
           static_cast<std::size_t>(_indices.Size());
  }
};

// Transform every radius at once.
//
// The whole stack goes to the grid as one batch, so the Wigner block for
// this upper index is streamed once for all the radii rather than once per
// radius, and the threading policy is the caller's as everywhere else.
//
// A batch shares grid, degree and upper index. A stack of one field's slices
// shares all three by construction, which is why the radial axis is the batch
// axis and no other axis of a three-dimensional problem is.
template <std::ptrdiff_t N, AngularGrid Grid, RealOrComplexValued Value>
auto Expand(const LayeredSpinField<N, Grid, Value>& field, std::ptrdiff_t lMax,
            Execution policy = Execution::Sequential()) {
  auto expansion =
      LayeredSpinExpansion<N, Grid, Value>(field.Radial(), field.Grid(), lMax);
  auto out = expansion.Data();
  field.Grid().ForwardTransformation(lMax, N, field.Data(), field.Batch(), out,
                                     expansion.Batch(), policy);
  return expansion;
}

template <std::ptrdiff_t N, AngularGrid Grid, RealOrComplexValued Value>
auto Evaluate(const LayeredSpinExpansion<N, Grid, Value>& expansion,
              Execution policy = Execution::Sequential()) {
  auto field =
      LayeredSpinField<N, Grid, Value>(expansion.Radial(), expansion.Grid());
  auto out = field.Data();
  expansion.Grid().InverseTransformation(expansion.MaxDegree(), N,
                                         expansion.Data(), expansion.Batch(),
                                         out, field.Batch(), policy);
  return field;
}

// Evaluate a two-dimensional expression at every radius.
//
// The explicit way up. There is no implicit conversion in either direction: a
// function that wants a stack must be given one, and a function that wants an
// angular field cannot be handed a one-slice stack by accident.
template <typename Expr, typename RadialGridType>
requires SpinFieldExpr<Expr>
auto Broadcast(const RadialGridType& radialGrid, const Expr& expr) {
  using E = Node<Expr>;
  auto stack = LayeredSpinField<E::UpperIndex, typename E::GridType,
                                typename E::Value>(radialGrid, expr.Grid());
  for (auto i : stack.RadiusIndices()) {
    auto slice = stack.Slice(i);
    expr.EvaluateInto(slice.Data());
  }
  return stack;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_LAYERED_SPIN_FIELD_GUARD_H
