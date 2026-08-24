#ifndef GSH_TRANS_RADIAL_MAJOR_GUARD_H
#define GSH_TRANS_RADIAL_MAJOR_GUARD_H

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "../Concepts.h"
#include "../Policies.h"
#include "RadialOperator.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                              The other layout                             //
//--------------------------------------------------------------------------//

// The same data with the two axes exchanged: `[(l, m)][r]` rather than
// `[r][(l, m)]`, so that a radial line is contiguous.
//
// Both layouts are wanted and neither is wanted always. The angular transform
// wants radius-major, because a batch is `nR` contiguous angular blocks. A
// radial solve at fixed degree and order wants radial-major, because the line
// it solves along is then a contiguous vector it can hand to a band solver.
// The preconditioner of a matrix-free 3-D operator is a sequence of exactly
// such solves, one per `(l, m)`, so both appear in the same inner loop.
//
// `ApplyRadially` already covers the case where a line is touched once: it
// gathers, applies and scatters, which is this repack done one line at a time
// and thrown away. What it cannot amortise is a line touched *many* times --
// an iterative solve, a factorisation applied repeatedly, a sweep over several
// operators. That is what this is for, and the choice between the two is a
// measurement rather than a principle.
//
// This is a buffer and a shape, not a field type. It has no grid, no upper
// index and no algebra, because nothing angular is meaningful once the angular
// axis has been shredded into lines. Repack, work on lines, repack back.

/// @tparam Stack The radius-major stack this is repacked from.
template <typename Stack>
class RadialMajor {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  using Scalar = typename std::remove_cvref_t<
      decltype(std::declval<const Stack&>()
                   .Data())>::value_type;  ///< The stack's value type.

  RadialMajor() = delete;

  /**
   * @brief Repacks a stack into radial lines.
   * @param stack The radius-major stack to transpose.
   * @param policy Whether the transpose may thread.
   */
  explicit RadialMajor(const Stack& stack,
                       Execution policy = Execution::Sequential())
      : _nR{stack.NumberOfRadii()},
        _lines{stack.SliceSize()},
        _data(static_cast<std::size_t>(_nR) *
              static_cast<std::size_t>(_lines)) {
    Transpose(stack.Data().data(), _data.data(), _nR, _lines, policy);
  }

  /** @brief How many radii the stack holds. */
  auto NumberOfRadii() const { return _nR; }
  /** @brief How many radial lines there are, one per element of a slice. */
  auto NumberOfLines() const { return _lines; }

  /// A buffer of the same shape, not transposed from anything. Scratch for an
  /// operator to write into: transposing a buffer whose contents are about to
  /// be overwritten is pure waste, and without this there is no way to avoid
  /// it.
  auto SameShape() const { return RadialMajor(_nR, _lines); }
  /** @brief How many elements are stored. */
  auto Size() const { return static_cast<Int>(_data.size()); }

  /** @brief The underlying buffer. */
  auto Data() { return std::span<Scalar>(_data); }
  /** @brief The underlying buffer. */
  auto Data() const { return std::span<const Scalar>(_data); }

  /// One radial line, contiguous. `j` is the position within a slice: for an
  /// expansion, what `CoefficientIndex(l, m)` returns.
  auto Line(Int j) {
    return Data().subspan(Offset(j), static_cast<std::size_t>(_nR));
  }

  /// The same, read-only.
  auto Line(Int j) const {
    return Data().subspan(Offset(j), static_cast<std::size_t>(_nR));
  }

  /** @brief Every line index. */
  auto LineIndices() const { return std::ranges::views::iota(Int{0}, _lines); }

  /// Refill from a stack of the shape this already has, without allocating.
  ///
  /// The loop this exists for repacks on every operator application, and
  /// allocating a buffer the size of the whole field each time -- with the page
  /// faults that first touching it brings -- would cost more than the transpose
  /// it is there to serve.
  void CopyFrom(const Stack& stack,
                Execution policy = Execution::Sequential()) {
    if (stack.NumberOfRadii() != _nR || stack.SliceSize() != _lines) {
      throw std::invalid_argument(
          "A radial-major buffer can only be refilled from a stack of the "
          "shape it already has");
    }
    Transpose(stack.Data().data(), _data.data(), _nR, _lines, policy);
  }

  /// Back to radius-major, into a stack of the right shape.
  void CopyInto(Stack& stack,
                Execution policy = Execution::Sequential()) const {
    if (stack.NumberOfRadii() != _nR || stack.SliceSize() != _lines) {
      throw std::invalid_argument(
          "A radial-major buffer can only be copied back into a stack of the "
          "shape it came from");
    }
    auto target = stack.Data();
    Transpose(_data.data(), target.data(), _lines, _nR, policy);
  }

 private:
  Int _nR;
  Int _lines;
  std::vector<Scalar> _data;

  RadialMajor(Int nR, Int lines)
      : _nR{nR},
        _lines{lines},
        _data(static_cast<std::size_t>(nR) * static_cast<std::size_t>(lines)) {}

  std::size_t Offset(Int j) const {
    if (j < 0 || j >= _lines) {
      throw std::invalid_argument("Line index outside the slice");
    }
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(_nR);
  }

  // The tile side, in elements. A pair of tiles has to sit in L1 alongside
  // whatever else is live, and a complex double is 16 bytes: 32 gives 16 KiB
  // per tile, which is already the whole of a 32 KiB L1, so 16 is the largest
  // that leaves room for both. Measured against 8, 16, 32 and 64.
  static Int TileSize() { return 16; }

  // A blocked transpose of `rows x cols` into `cols x rows`.
  //
  // Blocked rather than a plain double loop because one of the two sides is
  // always strided, and at the sizes this runs at -- a slice is tens of
  // thousands of coefficients -- the strided side misses on every element. The
  // tile is sized so that a block of both matrices sits in L1 together.
  template <typename T>
  static void Transpose(const T* in, T* out, Int rows, Int cols,
                        Execution policy) {
    const auto tile = TileSize();
    const auto threads =
        policy.IsParallel() && !omp_in_parallel()
            ? (policy.Threads() > 0 ? policy.Threads() : omp_get_max_threads())
            : 1;

    const auto block = [&](Int i0) {
      const auto iEnd = std::min(i0 + tile, rows);
      for (auto j0 = Int{0}; j0 < cols; j0 += tile) {
        const auto jEnd = std::min(j0 + tile, cols);
        for (auto i = i0; i < iEnd; i++) {
          for (auto j = j0; j < jEnd; j++) {
            out[j * rows + i] = in[i * cols + j];
          }
        }
      }
    };

    if (threads == 1) {
      for (auto i0 = Int{0}; i0 < rows; i0 += tile) block(i0);
    } else {
#pragma omp parallel for schedule(static) num_threads(threads)
      for (Int i0 = 0; i0 < rows; i0 += tile) block(i0);
    }
  }
};

// Apply a radial operator to every line of a radial-major buffer.
//
// The difference from `ApplyRadially` is the whole point of the layout: there
// is no gather and no scatter, because the line is already contiguous. The
// operator writes into a scratch line and the result is copied back, so it
// still need not handle aliasing.
template <typename Stack, typename Op>
void ApplyToLines(const RadialMajor<Stack>& in, RadialMajor<Stack>& out,
                  const Op& op, Execution policy = Execution::Sequential()) {
  using Int = std::ptrdiff_t;
  using Scalar = typename RadialMajor<Stack>::Scalar;

  if (in.NumberOfRadii() != out.NumberOfRadii() ||
      in.NumberOfLines() != out.NumberOfLines()) {
    throw std::invalid_argument(
        "A radial operator acts along the radial axis alone, so its argument "
        "and its result must have the same shape");
  }

  const auto lines = in.NumberOfLines();
  const auto threads =
      policy.IsParallel() && !omp_in_parallel()
          ? (policy.Threads() > 0 ? policy.Threads() : omp_get_max_threads())
          : 1;

  const auto run = [&](Int j) {
    op(std::span<const Scalar>(in.Line(j)), std::span<Scalar>(out.Line(j)));
  };

  if (threads == 1) {
    for (auto j = Int{0}; j < lines; j++) run(j);
  } else {
#pragma omp parallel for schedule(static) num_threads(threads)
    for (Int j = 0; j < lines; j++) run(j);
  }
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_RADIAL_MAJOR_GUARD_H
