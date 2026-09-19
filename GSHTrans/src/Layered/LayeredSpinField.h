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
#include "../Indexing.h"
#include "../Policies.h"
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
/// is a SpinFieldView, an ordinary spin-weighted node, so the index algebra,
/// the evaluation and the aliasing theorem all lift unchanged and there is no
/// second expression system to keep consistent with the first. That is what
/// a view being admissible everywhere an owning field is was for.
///
/// The layout is radius-major because that is what the applications this exists
/// for already use, and because it makes the radial axis a batch: the stack of
/// one field's slices is `nR` blocks of `FieldSize`, which the transform's
/// (count, stride, dist) descriptor reads directly.
///
/// A one-slice stack is *not* an angular field. The idiom stays available to
/// application code, but the two are distinct types here, so that a function
/// taking one cannot silently be handed the other.
template <std::ptrdiff_t N_, AngularGrid Grid_,
          RealOrComplexValued Value_ = ComplexValued>
class LayeredSpinField {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.

  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = N_;
  using Value = Value_;    ///< Whether the samples are real-valued or complex.
  using GridType = Grid_;  ///< The angular grid this is defined on.
  using Real = typename Grid_::Real;   ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  /// The value type: Real when real-valued, Complex otherwise.
  using Scalar = ScalarFor<Real, Value>;
  using RadialGridType = RadialGrid<Real>;  ///< The radial grid type.
  using SliceType = SpinFieldView<N_, Grid_, Value_>;  ///< One angular field.
  using ConstSliceType =
      ConstSpinFieldView<N_, Grid_, Value_>;  ///< One angular field, read-only.

  static_assert(std::same_as<Value, ComplexValued> or UpperIndex == 0,
                "A spin-weighted field can be real-valued only at upper index "
                "zero");

  LayeredSpinField() = delete;

  /**
   * @brief A zero stack on the product of the two grids.
   * @param radialGrid The radii.
   * @param grid The angular grid every slice is on.
   */
  LayeredSpinField(RadialGridType radialGrid, GridType grid)
      : radialGrid_{std::move(radialGrid)},
        grid_{std::move(grid)},
        data_(static_cast<std::size_t>(radialGrid_.NumberOfRadii()) *
              static_cast<std::size_t>(grid_.FieldSize())) {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return grid_; }
  /** @brief The radial grid this is defined on. */
  const RadialGridType& Radial() const { return radialGrid_; }

  /** @brief How many radii the stack holds. */
  auto NumberOfRadii() const { return radialGrid_.NumberOfRadii(); }
  /** @brief Indices of the stored radii. */
  auto RadiusIndices() const { return radialGrid_.RadiusIndices(); }
  /** @brief How many samples one angular field holds. */
  auto FieldSize() const { return grid_.FieldSize(); }
  /** @brief How many elements are stored. */
  auto Size() const { return static_cast<Int>(data_.size()); }

  // The uniform names a radial operator sees. A field's slice is a set of
  // angular points and an expansion's is a set of coefficients, but the radial
  // axis does not care which: it runs over `NumberOfRadii()` values
  // `SliceSize()` apart, and that is all `ApplyRadially` needs to know about
  // either.
  /** @brief How many elements one radial slice holds. */
  auto SliceSize() const { return FieldSize(); }
  /** @brief A zero stack of the same shape, which is what an operator needs
   * to write into. */
  auto SameShape() const { return LayeredSpinField(radialGrid_, grid_); }

  /// The same field on a different set of radii, which is what resampling
  /// needs and what SameShape cannot give: a radial operator maps a stack to
  /// one of the same shape, and changing nR is precisely not that.
  auto SameShapeOn(RadialGridType radialGrid) const {
    return LayeredSpinField(std::move(radialGrid), grid_);
  }

  /** @brief The underlying buffer. */
  auto Data() { return std::span<Scalar>(data_); }
  /** @brief The underlying buffer. */
  auto Data() const { return std::span<const Scalar>(data_); }

  /// One angular field, as a view over this stack's storage. Writing through it
  /// writes the stack.
  auto Slice(Int i) {
    return SliceType(grid_, Data().subspan(Offset(i), SliceExtent()));
  }

  /// The same, read-only.
  auto Slice(Int i) const {
    return ConstSliceType(grid_, Data().subspan(Offset(i), SliceExtent()));
  }

  /// The whole stack as the transform's batch: nR fields, each contiguous,
  /// FieldSize apart.
  auto Batch() const {
    return GSHTrans::Batch::Contiguous(NumberOfRadii(), FieldSize());
  }

 private:
  RadialGridType radialGrid_;
  GridType grid_;
  FFTWpp::vector<Scalar> data_;

  std::size_t SliceExtent() const {
    return static_cast<std::size_t>(grid_.FieldSize());
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
/// This is the [r][(l,m)] of the two layouts a layered code needs. The other,
/// [(l,m)][r], is what a radial solve at fixed degree and order wants, and the
/// repack between them is a separate step.
template <std::ptrdiff_t N_, AngularGrid Grid_,
          RealOrComplexValued Value_ = ComplexValued>
class LayeredSpinExpansion {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.

  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = N_;
  using Value = Value_;    ///< Whether the samples are real-valued or complex.
  using GridType = Grid_;  ///< The angular grid this is defined on.
  using Real = typename Grid_::Real;   ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  using RadialGridType = RadialGrid<Real>;  ///< The radial grid type.
  using MRange =
      std::conditional_t<std::same_as<Value_, RealValued>, NonNegative,
                         All>;  ///< Which orders are stored.

  LayeredSpinExpansion() = delete;

  /**
   * @brief A zero expansion on the product of the two grids.
   * @param radialGrid The radii.
   * @param grid The angular grid the coefficients belong to.
   * @param lMax The largest degree stored.
   */
  LayeredSpinExpansion(RadialGridType radialGrid, GridType grid, Int lMax)
      : radialGrid_{std::move(radialGrid)},
        grid_{std::move(grid)},
        indices_{Checked(lMax), lMax, UpperIndex},
        data_(static_cast<std::size_t>(radialGrid_.NumberOfRadii()) *
              static_cast<std::size_t>(indices_.Size())) {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return grid_; }
  /** @brief The radial grid this is defined on. */
  const RadialGridType& Radial() const { return radialGrid_; }

  /** @brief How many radii the stack holds. */
  auto NumberOfRadii() const { return radialGrid_.NumberOfRadii(); }
  /** @brief Indices of the stored radii. */
  auto RadiusIndices() const { return radialGrid_.RadiusIndices(); }
  /** @brief The largest degree stored. */
  auto MaxDegree() const { return indices_.MaxDegree(); }
  /** @brief The smallest degree stored. */
  auto MinDegree() const { return indices_.MinDegree(); }
  /** @brief Every degree stored. */
  auto Degrees() const { return indices_.Degrees(); }
  /** @brief Every order stored at degree @p l. */
  auto Orders(Int l) const {
    return GSHSubIndices<MRange>(l, MaxDegree()).Orders();
  }
  /** @brief How many coefficients one block holds. */
  auto CoefficientSize() const { return static_cast<Int>(indices_.Size()); }
  /** @brief How many elements are stored. */
  auto Size() const { return static_cast<Int>(data_.size()); }

  /** @brief How many elements one radial slice holds. */
  auto SliceSize() const { return CoefficientSize(); }

  /// Where a degree and order sit within one radius's block. The radial-major
  /// repack needs this: once the layout is [(l, m)][r] a line is addressed by
  /// its position in the block and no longer by (l, m) directly.
  auto CoefficientIndex(Int l, Int m) const {
    return static_cast<Int>(indices_.Index(l, m));
  }

  /** @brief A zero expansion of the same shape, which is what an operator
   * needs to write into. */
  auto SameShape() const {
    return LayeredSpinExpansion(radialGrid_, grid_, MaxDegree());
  }

  /// The same expansion on a different set of radii, which is what resampling
  /// needs and what SameShape cannot give.
  auto SameShapeOn(RadialGridType radialGrid) const {
    return LayeredSpinExpansion(std::move(radialGrid), grid_, MaxDegree());
  }

  /** @brief The underlying buffer. */
  auto Data() { return std::span<Complex>(data_); }
  /** @brief The underlying buffer. */
  auto Data() const { return std::span<const Complex>(data_); }

  /// The coefficient at one radius.
  Complex operator[](Int i, Int l, Int m) const {
    return data_[Offset(i) + static_cast<std::size_t>(indices_.Index(l, m))];
  }
  /// The same, writable.
  Complex& operator[](Int i, Int l, Int m) {
    return data_[Offset(i) + static_cast<std::size_t>(indices_.Index(l, m))];
  }

  /// The whole stack as the transform's batch: nR blocks, each contiguous,
  /// CoefficientSize() apart.
  auto Batch() const {
    return GSHTrans::Batch::Contiguous(NumberOfRadii(), CoefficientSize());
  }

 private:
  RadialGridType radialGrid_;
  GridType grid_;
  GSHIndices<MRange> indices_;
  FFTWpp::vector<Complex> data_;

  static Int Checked(Int lMax) {
    if (lMax < (UpperIndex < 0 ? -UpperIndex : UpperIndex)) {
      throw std::invalid_argument(
          "An expansion's degree cannot be below its upper index");
    }
    return lMax;
  }

  std::size_t Offset(Int i) const {
    if (i < 0 || i >= NumberOfRadii()) {
      throw std::invalid_argument(
          "Radius index " + std::to_string(i) + " is outside the " +
          std::to_string(NumberOfRadii()) + " radii of this grid");
    }
    return static_cast<std::size_t>(i) *
           static_cast<std::size_t>(indices_.Size());
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
  auto stack =
      LayeredSpinField<E::UpperIndex, typename E::GridType, typename E::Value>(
          radialGrid, expr.Grid());
  for (auto i : stack.RadiusIndices()) {
    auto slice = stack.Slice(i);
    expr.EvaluateInto(slice.Data());
  }
  return stack;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_LAYERED_SPIN_FIELD_GUARD_H
