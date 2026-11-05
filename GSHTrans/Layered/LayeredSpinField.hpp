#pragma once

#include <FFTWpp/Core>
#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "../Concepts.hpp"
#include "../Expansion/SpinExpansion.hpp"
#include "../Indexing.hpp"
#include "../Parallel.hpp"
#include "../Policies.hpp"
#include "../SpinField/SpinField.hpp"
#include "../SpinField/SpinFieldOverloads.hpp"
#include "../SpinField/SpinFieldView.hpp"
#include "../SpinField/SpinWeighted.hpp"
#include "RadialGrid.hpp"

namespace GSHTrans {

/// Specialised by every type that can stand in a layered expression
/// (LayeredExpressions.hpp). Declared here so that the owning stack can be
/// assigned one.
template <typename T>
struct IsLayeredNodeTrait : std::false_type {};

/// A field on the product of a radial grid with an angular one: nR angular
/// fields laid end to end, radius outermost, each slice contiguous in the
/// canonical angular order.
///
/// Two-dimensional is the primitive and three-dimensional is a stack. The
/// angular field is never wrapped and a slice is not a new kind of object: it
/// is a SpinFieldView, an ordinary spin-weighted node, so the index algebra,
/// the evaluation and the aliasing theorem all lift unchanged and there is no
/// second expression system to keep consistent with the first. That follows
/// from a view being admissible everywhere an owning field is.
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

  LayeredSpinField(const LayeredSpinField&) = default;
  LayeredSpinField(LayeredSpinField&&) = default;
  LayeredSpinField& operator=(const LayeredSpinField&) = default;
  LayeredSpinField& operator=(LayeredSpinField&&) = default;

  /// Evaluates a layered expression into this stack, sequentially; `Assign`
  /// is the threaded form. The stack may appear in the expression.
  template <typename Expr>
  requires IsLayeredNodeTrait<std::remove_cvref_t<Expr>>::value
  LayeredSpinField& operator=(const Expr& expr) {
    Assign(*this, expr);
    return *this;
  }

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

//--------------------------------------------------------------------------//
//                        Views over a caller's storage                      //
//--------------------------------------------------------------------------//

/// A layered field over storage the caller owns: the same radius-major layout
/// as LayeredSpinField, slice `i` at offset `i * FieldSize()`, over a span.
///
/// What a code that keeps its own arrays needs in order to use the layered
/// transforms and expressions without copying into the library's types: a
/// view is read and written in place. Like the two-dimensional SpinFieldView
/// it is a handle -- two grid handles and a span -- and is held by value in
/// an expression. The storage must outlive it.
template <std::ptrdiff_t N_, AngularGrid Grid_,
          RealOrComplexValued Value_ = ComplexValued,
          typename Element_ = ScalarFor<typename Grid_::Real, Value_>>
class LayeredSpinFieldView {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.

  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = N_;
  using Value = Value_;    ///< Whether the samples are real-valued or complex.
  using GridType = Grid_;  ///< The angular grid this is defined on.
  using Real = typename Grid_::Real;   ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  /// The value type: Real when real-valued, Complex otherwise.
  using Scalar = std::remove_const_t<Element_>;
  using RadialGridType = RadialGrid<Real>;  ///< The radial grid type.
  /// One angular field: a two-dimensional view over this view's storage.
  using SliceType = SpinFieldView<N_, Grid_, Value_, Element_>;

  static_assert(std::same_as<Scalar, ScalarFor<Real, Value>>);

  LayeredSpinFieldView() = delete;

  /**
   * @brief A view of @p data as a layered field.
   * @param radialGrid The radii.
   * @param grid The angular grid every slice is on.
   * @param data Exactly `nR * FieldSize()` values, radius-major.
   * @throws std::invalid_argument if @p data is not that size, or the grid
   * does not carry this upper index.
   */
  LayeredSpinFieldView(RadialGridType radialGrid, GridType grid,
                       std::span<Element_> data)
      : radialGrid_{std::move(radialGrid)},
        grid_{std::move(grid)},
        data_{data} {
    if (!std::ranges::contains(grid_.UpperIndices(), UpperIndex)) {
      throw std::invalid_argument("This grid does not carry upper index " +
                                  std::to_string(UpperIndex));
    }
    const auto needed = static_cast<std::size_t>(NumberOfRadii()) *
                        static_cast<std::size_t>(FieldSize());
    if (data_.size() != needed) {
      throw std::invalid_argument(
          "A layered view over " + std::to_string(data_.size()) +
          " values does not hold the " + std::to_string(needed) +
          " that these radii and this grid need");
    }
  }

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
  /** @brief How many elements are viewed. */
  auto Size() const { return static_cast<Int>(data_.size()); }
  /** @brief How many elements one radial slice holds. */
  auto SliceSize() const { return FieldSize(); }

  /** @brief The storage viewed. */
  auto Data() const { return data_; }

  /// One angular field, as a two-dimensional view over the same storage.
  auto Slice(Int i) const {
    if (i < 0 || i >= NumberOfRadii()) {
      throw std::invalid_argument(
          "Radius index " + std::to_string(i) + " is outside the " +
          std::to_string(NumberOfRadii()) + " radii of this grid");
    }
    return SliceType(grid_,
                     data_.subspan(static_cast<std::size_t>(i) *
                                       static_cast<std::size_t>(FieldSize()),
                                   static_cast<std::size_t>(FieldSize())));
  }

  /// The whole stack as the transform's batch: nR fields, each contiguous,
  /// FieldSize apart.
  auto Batch() const {
    return GSHTrans::Batch::Contiguous(NumberOfRadii(), FieldSize());
  }

 private:
  RadialGridType radialGrid_;
  GridType grid_;
  std::span<Element_> data_;
};

/// A layered field view nobody may write through.
template <std::ptrdiff_t N, AngularGrid Grid,
          RealOrComplexValued Value = ComplexValued>
using ConstLayeredSpinFieldView =
    LayeredSpinFieldView<N, Grid, Value,
                         const ScalarFor<typename Grid::Real, Value>>;

/// A layered expansion over storage the caller owns: block `i` at offset
/// `i * CoefficientSize()`, as in LayeredSpinExpansion.
template <std::ptrdiff_t N_, AngularGrid Grid_,
          RealOrComplexValued Value_ = ComplexValued,
          typename Element_ = std::complex<typename Grid_::Real>>
class LayeredSpinExpansionView {
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
  /// One radius's coefficients: a two-dimensional expansion view.
  using SliceType = SpinExpansionBase<N_, Grid_, Value_, Element_>;

  static_assert(std::same_as<std::remove_const_t<Element_>, Complex>);

  LayeredSpinExpansionView() = delete;

  /**
   * @brief A view of @p data as a layered expansion to degree @p lMax.
   * @param radialGrid The radii.
   * @param grid The angular grid the coefficients belong to.
   * @param lMax The largest degree stored.
   * @param data Exactly `nR * CoefficientSize()` values, radius-major.
   * @throws std::invalid_argument if @p data is not that size, or @p lMax is
   * below the upper index.
   */
  LayeredSpinExpansionView(RadialGridType radialGrid, GridType grid, Int lMax,
                           std::span<Element_> data)
      : radialGrid_{std::move(radialGrid)},
        grid_{std::move(grid)},
        indices_{Checked(lMax), lMax, UpperIndex},
        data_{data} {
    const auto needed = static_cast<std::size_t>(NumberOfRadii()) *
                        static_cast<std::size_t>(indices_.Size());
    if (data_.size() != needed) {
      throw std::invalid_argument(
          "A layered expansion view over " + std::to_string(data_.size()) +
          " coefficients does not hold the " + std::to_string(needed) +
          " that these radii and this degree need");
    }
  }

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
  /** @brief How many elements are viewed. */
  auto Size() const { return static_cast<Int>(data_.size()); }
  /** @brief How many elements one radial slice holds. */
  auto SliceSize() const { return CoefficientSize(); }
  /** @brief Where a degree and order sit within one radius's block. */
  auto CoefficientIndex(Int l, Int m) const {
    return static_cast<Int>(indices_.Index(l, m));
  }

  /** @brief The storage viewed. */
  auto Data() const { return data_; }

  /// One radius's coefficients, as a two-dimensional expansion view.
  auto Slice(Int i) const {
    return SliceType(
        grid_, MaxDegree(),
        data_.subspan(Offset(i), static_cast<std::size_t>(CoefficientSize())));
  }

  /// The coefficient at one radius.
  Complex operator[](Int i, Int l, Int m) const {
    return data_[Offset(i) + static_cast<std::size_t>(indices_.Index(l, m))];
  }
  /// The same, writable. Present only on a view over mutable storage.
  Element_& operator[](Int i, Int l, Int m)
  requires(not std::is_const_v<Element_>)
  {
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
  std::span<Element_> data_;

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

/// A layered expansion view nobody may write through.
template <std::ptrdiff_t N, AngularGrid Grid,
          RealOrComplexValued Value = ComplexValued>
using ConstLayeredSpinExpansionView =
    LayeredSpinExpansionView<N, Grid, Value,
                             const std::complex<typename Grid::Real>>;

//--------------------------------------------------------------------------//
//                     What the layered transforms accept                    //
//--------------------------------------------------------------------------//

/// Specialised by the types that hold a layered field's samples, owning or
/// viewing: what a layered transform reads from and writes to on the grid
/// side. A trait and not a structural concept, for the reason IsTerminalTrait
/// is one: the interface does not say what a type is.
template <typename T>
struct IsLayeredFieldStorageTrait : std::false_type {};
template <std::ptrdiff_t N, AngularGrid G, RealOrComplexValued V>
struct IsLayeredFieldStorageTrait<LayeredSpinField<N, G, V>> : std::true_type {
};
template <std::ptrdiff_t N, AngularGrid G, RealOrComplexValued V, typename E>
struct IsLayeredFieldStorageTrait<LayeredSpinFieldView<N, G, V, E>>
    : std::true_type {};

/// The same for the coefficient side.
template <typename T>
struct IsLayeredExpansionStorageTrait : std::false_type {};
template <std::ptrdiff_t N, AngularGrid G, RealOrComplexValued V>
struct IsLayeredExpansionStorageTrait<LayeredSpinExpansion<N, G, V>>
    : std::true_type {};
template <std::ptrdiff_t N, AngularGrid G, RealOrComplexValued V, typename E>
struct IsLayeredExpansionStorageTrait<LayeredSpinExpansionView<N, G, V, E>>
    : std::true_type {};

/** @brief A layered field's samples, owning or viewed. */
template <typename T>
concept LayeredFieldStorage =
    IsLayeredFieldStorageTrait<std::remove_cvref_t<T>>::value;

/** @brief A layered expansion's coefficients, owning or viewed. */
template <typename T>
concept LayeredExpansionStorage =
    IsLayeredExpansionStorageTrait<std::remove_cvref_t<T>>::value;

namespace LayeredDetails {

// Storage of the given kind that can be written, with the same upper index,
// grid type and value kind as `Like`.
template <typename T, typename Like>
concept WritableLike =
    (std::remove_cvref_t<T>::UpperIndex == Like::UpperIndex) &&
    std::same_as<typename std::remove_cvref_t<T>::GridType,
                 typename Like::GridType> &&
    std::same_as<typename std::remove_cvref_t<T>::Value,
                 typename Like::Value> &&
    requires(std::remove_cvref_t<T>& t) {
      requires not std::is_const_v<typename decltype(t.Data())::element_type>;
    };

// The two stacks of a transform must be on the same grids: a stack on other
// radii would claim values at radii they were not computed for, and one on
// another angular grid of the same size would pair samples with the wrong
// points. Equal sizes would let either through in silence.
template <typename A, typename B>
void CheckSameGrids(const A& a, const B& b) {
  if (a.Radial().Identity() != b.Radial().Identity() ||
      a.Grid().Identity() != b.Grid().Identity()) {
    throw std::invalid_argument(
        "A layered field and the expansion it is transformed to or from must "
        "be on the same radial and angular grids");
  }
}

}  // namespace LayeredDetails

//--------------------------------------------------------------------------//
//                         Transforming every radius                         //
//--------------------------------------------------------------------------//

// The whole stack goes to the grid as one batch, so the Wigner block for
// this upper index is streamed once for all the radii rather than once per
// radius, and the threading policy is the caller's as everywhere else.
//
// A batch shares grid, degree and upper index. A stack of one field's slices
// shares all three by construction, which is why the radial axis is the batch
// axis and no other axis of a three-dimensional problem is.
//
// Each comes in two forms, as the two-dimensional Expand and Evaluate do: one
// writing into a stack the caller already holds, which allocates nothing and
// takes its degree from that stack, and one returning a new stack, which is
// the first applied to a new result. Either side may be an owning stack or a
// view over the caller's storage.

/**
 * @brief Expands every radius of @p field into @p out, to @p out's degree.
 * @param field A layered field or view.
 * @param out A layered expansion or writable view, on the same grids.
 * @param policy Whether the transform may thread.
 * @throws std::invalid_argument if the two are not on the same grids.
 */
template <LayeredFieldStorage Field, LayeredExpansionStorage Out>
requires LayeredDetails::WritableLike<Out, std::remove_cvref_t<Field>>
void Expand(const Field& field, Out&& out,
            Execution policy = Execution::Sequential()) {
  constexpr auto N = std::remove_cvref_t<Field>::UpperIndex;
  LayeredDetails::CheckSameGrids(field, out);
  auto target = out.Data();
  field.Grid().ForwardTransformation(out.MaxDegree(), N, field.Data(),
                                     field.Batch(), target, out.Batch(),
                                     policy);
}

/** @brief Expands every radius of @p field to degree @p lMax. */
template <LayeredFieldStorage Field>
auto Expand(const Field& field, std::ptrdiff_t lMax,
            Execution policy = Execution::Sequential()) {
  using F = std::remove_cvref_t<Field>;
  auto expansion = LayeredSpinExpansion<F::UpperIndex, typename F::GridType,
                                        typename F::Value>(field.Radial(),
                                                           field.Grid(), lMax);
  Expand(field, expansion, policy);
  return expansion;
}

/**
 * @brief Evaluates every radius of @p expansion into @p out.
 * @param expansion A layered expansion or view.
 * @param out A layered field or writable view, on the same grids.
 * @param policy Whether the transform may thread.
 * @throws std::invalid_argument if the two are not on the same grids.
 */
template <LayeredExpansionStorage Expansion, LayeredFieldStorage Out>
requires LayeredDetails::WritableLike<Out, std::remove_cvref_t<Expansion>>
void Evaluate(const Expansion& expansion, Out&& out,
              Execution policy = Execution::Sequential()) {
  constexpr auto N = std::remove_cvref_t<Expansion>::UpperIndex;
  LayeredDetails::CheckSameGrids(out, expansion);
  auto target = out.Data();
  expansion.Grid().InverseTransformation(expansion.MaxDegree(), N,
                                         expansion.Data(), expansion.Batch(),
                                         target, out.Batch(), policy);
}

/** @brief Evaluates every radius of @p expansion. */
template <LayeredExpansionStorage Expansion>
auto Evaluate(const Expansion& expansion,
              Execution policy = Execution::Sequential()) {
  using E = std::remove_cvref_t<Expansion>;
  auto field =
      LayeredSpinField<E::UpperIndex, typename E::GridType, typename E::Value>(
          expansion.Radial(), expansion.Grid());
  Evaluate(expansion, field, policy);
  return field;
}

// Evaluate a two-dimensional expression at every radius.
//
// The explicit way up. There is no implicit conversion in either direction: a
// function that wants a stack must be given one, and a function that wants an
// angular field cannot be handed a one-slice stack by accident.
//
// Under a parallel policy the slices are evaluated concurrently, which the
// expression layer allows: evaluation is const and stateless, and each slice
// is a disjoint target.
template <typename Expr, typename RadialGridType>
requires SpinFieldExpr<Expr>
auto Broadcast(const RadialGridType& radialGrid, const Expr& expr,
               Execution policy = Execution::Sequential()) {
  using E = Node<Expr>;
  auto stack =
      LayeredSpinField<E::UpperIndex, typename E::GridType, typename E::Value>(
          radialGrid, expr.Grid());
  ForEach(policy, stack.NumberOfRadii(), [&](std::ptrdiff_t i) {
    auto slice = stack.Slice(i);
    expr.EvaluateInto(slice.Data());
  });
  return stack;
}

}  // namespace GSHTrans
