#pragma once

#include <algorithm>
#include <cassert>
#include <complex>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "../Concepts.hpp"
#include "SpinWeighted.hpp"

namespace GSHTrans {

/// A non-owning field: someone else's storage, plus a grid handle.
///
/// This is what makes the tensor layer possible. A TensorField owns one
/// contiguous buffer and hands out components as views at the correct upper
/// index, rather than holding a tuple of separately allocated fields; the
/// radial slices of a layered field are the same idea. Both need a view to be
/// admissible wherever an owning field is, which is why operator[] returns by
/// value on every node and why the grid is a handle rather than a reference.
///
/// A view is a handle and not storage: a grid handle, a span and a stride. An
/// expression therefore holds one **by value**, whatever value category it
/// arrived with, exactly as it holds another expression. It used to hold an
/// lvalue view by reference, as it does an owning field, on the argument that
/// a view "names storage" -- but the storage a view names is not the view, and
/// a view is what one naturally binds to a local, so
///
///     auto u = t.Component<1>();  auto w = t.Component<-1>();  return u * w;
///
/// returned an expression referring to two dead locals. The copy that avoids
/// it is a shared_ptr and a span, once per expression and never per element.
template <std::ptrdiff_t N_, AngularGrid Grid_,
          RealOrComplexValued Value_ = ComplexValued,
          typename Element_ = ScalarFor<typename Grid_::Real, Value_>>
class SpinFieldView {
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

  static_assert(std::same_as<Scalar, ScalarFor<Real, Value>>);
  static_assert(std::same_as<Value, ComplexValued> or UpperIndex == 0,
                "A spin-weighted field can be real-valued only at upper index "
                "zero");
  static_assert(not std::same_as<typename Grid_::NRange, NonNegative> or
                    UpperIndex >= 0,
                "This grid stores only non-negative upper indices");

  SpinFieldView() = delete;

  /// The stride is how far apart successive samples are, and it defaults to
  /// one because most views are over contiguous storage.
  ///
  /// It is not one when a tensor field is laid out point by point: there each
  /// component's samples are separated by the number of components, and a view
  /// is the only way to hand that component to the field algebra. The
  /// transform needs no repack either, since a batch is described by (count,
  /// stride, dist) -- so the layout stays a choice rather than becoming a
  /// precondition.
  SpinFieldView(GridType grid, std::span<Element_> data, Int stride = 1)
      : grid_{std::move(grid)}, data_{data}, stride_{stride} {
    if (!std::ranges::contains(grid_.UpperIndices(), UpperIndex)) {
      throw std::invalid_argument("This grid does not carry upper index " +
                                  std::to_string(UpperIndex));
    }
    if (stride_ < 1) {
      throw std::invalid_argument("A view's stride must be positive");
    }
    const auto span =
        static_cast<std::size_t>((grid_.FieldSize() - 1) * stride_ + 1);
    if (data_.size() < span) {
      throw std::invalid_argument(
          "A view over " + std::to_string(data_.size()) + " values at stride " +
          std::to_string(stride_) + " does not cover this grid's " +
          std::to_string(grid_.FieldSize()) + " points");
    }
  }

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return grid_; }

  /** @brief The sample at the grid point @p iTheta, @p iPhi. */
  Scalar operator[](Int iTheta, Int iPhi) const {
    return data_[FlatIndex(iTheta, iPhi)];
  }

  /// Present only on a view over mutable storage.
  Element_& operator[](Int iTheta, Int iPhi)
  requires(not std::is_const_v<Element_>)
  {
    return data_[FlatIndex(iTheta, iPhi)];
  }

  /**
   * @brief Writes every sample into @p target, unstriding as it goes.
   * @throws std::invalid_argument if @p target is not the field's size.
   */
  template <typename S>
  requires std::convertible_to<Scalar, S>
  void EvaluateInto(std::span<S> target) const {
    const auto size = static_cast<std::size_t>(Size());
    if (target.size() != size) {
      throw std::invalid_argument(
          "Evaluation target has size " + std::to_string(target.size()) +
          ", but this field has " + std::to_string(size) + " points");
    }
    if (stride_ == 1) {
      std::ranges::copy(data_.first(size), target.begin());
      return;
    }
    for (auto i = std::size_t{0}; i < size; i++) {
      target[i] = static_cast<S>(data_[i * static_cast<std::size_t>(stride_)]);
    }
  }

  // The number of samples, which is the grid's point count and not the extent
  // of the storage those samples are spread over.
  /** @brief How many elements are stored. */
  auto Size() const { return grid_.FieldSize(); }
  /** @brief The underlying buffer. */
  auto Data() const { return data_; }
  /** @brief How many elements separate successive samples. */
  auto Stride() const { return stride_; }

 private:
  GridType grid_;
  std::span<Element_> data_;
  Int stride_;

  Int FlatIndex(Int iTheta, Int iPhi) const {
    const auto nPhi = grid_.NumberOfLongitudes();
    assert(iTheta >= 0 && iTheta < grid_.NumberOfCoLatitudes());
    assert(iPhi >= 0 && iPhi < nPhi);
    return (iTheta * nPhi + iPhi) * stride_;
  }
};

// A view over storage nobody may write through.
template <std::ptrdiff_t N, AngularGrid Grid,
          RealOrComplexValued Value = ComplexValued>
using ConstSpinFieldView =
    SpinFieldView<N, Grid, Value, const ScalarFor<typename Grid::Real, Value>>;

}  // namespace GSHTrans
