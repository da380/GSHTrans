#pragma once

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "../Concepts.hpp"
#include "../OpenMP.hpp"
#include "../Policies.hpp"
#include "../Utility.hpp"
#include "RadialOperator.hpp"

// The pragmas below are for a compiler that was asked for OpenMP. One that was
// not ignores them, which is right, and warns that it has, which in a
// header-only library is a warning in the caller's build. So the warning is
// turned off for this file alone, and only where OpenMP is off: see OpenMP.hpp.
#ifndef _OPENMP
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

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
// This is a buffer and a shape, not a field type. It has no angular grid, no
// upper index and no algebra, because nothing angular is meaningful once the
// angular axis has been shredded into lines. Repack, work on lines, repack
// back. What it does keep is the *radial* grid, since a line is a function of
// radius and an operator applied to it was built on particular radii.

/// @tparam Stack The radius-major stack this is repacked from.
template <typename Stack>
class RadialMajor {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  using Scalar = typename std::remove_cvref_t<
      decltype(std::declval<const Stack&>()
                   .Data())>::value_type;  ///< The stack's value type.
  /// The radial grid type of the stack.
  using RadialGridType =
      std::remove_cvref_t<decltype(std::declval<const Stack&>().Radial())>;

  RadialMajor() = delete;

  /**
   * @brief Repacks a stack into radial lines.
   * @param stack The radius-major stack to transpose.
   * @param policy Whether the transpose may thread.
   */
  explicit RadialMajor(const Stack& stack,
                       Execution policy = Execution::Sequential())
      : radial_{stack.Radial()},
        nR_{stack.NumberOfRadii()},
        lines_{stack.SliceSize()},
        data_(static_cast<std::size_t>(nR_) *
              static_cast<std::size_t>(lines_)) {
    Transpose(stack.Data().data(), data_.data(), nR_, lines_, policy);
  }

  /**
   * @brief An empty buffer of a given shape, to be filled by whoever asked.
   * @details For filling *directly* -- by a transform writing radial lines,
   * which is what ExpandToLines does -- so that the radius-major stack this
   * would otherwise be transposed from need never exist.
   * @param radial The radial grid the lines run along.
   * @param lines How many lines: one per element of a slice of the stack.
   */
  static RadialMajor OfShape(RadialGridType radial, Int lines) {
    if (lines < 1) {
      throw std::invalid_argument(
          "A radial-major buffer needs at least one line");
    }
    const auto nR = radial.NumberOfRadii();
    return RadialMajor(std::move(radial), nR, lines);
  }

  /** @brief The radial grid the lines run along. */
  const RadialGridType& Radial() const { return radial_; }

  /// The buffer as the transform's batch. Element j of radius i is at
  /// j * nR + i, which is `nR` fields interleaved element by element -- so the
  /// transform can read and write radial lines as they stand, and no
  /// transpose is needed to get in or out of this layout.
  auto Batch() const { return GSHTrans::Batch::Interleaved(nR_, nR_); }
  /** @brief How many radii the stack holds. */
  auto NumberOfRadii() const { return nR_; }
  /** @brief How many radial lines there are, one per element of a slice. */
  auto NumberOfLines() const { return lines_; }

  /// A buffer of the same shape, not transposed from anything. Scratch for an
  /// operator to write into: transposing a buffer whose contents are about to
  /// be overwritten is pure waste, and without this there is no way to avoid
  /// it.
  auto SameShape() const { return RadialMajor(radial_, nR_, lines_); }
  /** @brief How many elements are stored. */
  auto Size() const { return static_cast<Int>(data_.size()); }

  /** @brief The underlying buffer. */
  auto Data() { return std::span<Scalar>(data_); }
  /** @brief The underlying buffer. */
  auto Data() const { return std::span<const Scalar>(data_); }

  /// One radial line, contiguous. `j` is the position within a slice: for an
  /// expansion, what `CoefficientIndex(l, m)` returns.
  auto Line(Int j) {
    return Data().subspan(Offset(j), static_cast<std::size_t>(nR_));
  }

  /// The same, read-only.
  auto Line(Int j) const {
    return Data().subspan(Offset(j), static_cast<std::size_t>(nR_));
  }

  /** @brief Every line index. */
  auto LineIndices() const { return std::ranges::views::iota(Int{0}, lines_); }

  /// Refill from a stack of the shape this already has, without allocating.
  ///
  /// The loop this exists for repacks on every operator application, and
  /// allocating a buffer the size of the whole field each time -- with the page
  /// faults that first touching it brings -- would cost more than the transpose
  /// it is there to serve.
  void CopyFrom(const Stack& stack,
                Execution policy = Execution::Sequential()) {
    if (stack.NumberOfRadii() != nR_ || stack.SliceSize() != lines_) {
      throw std::invalid_argument(
          "A radial-major buffer can only be refilled from a stack of the "
          "shape it already has");
    }
    Transpose(stack.Data().data(), data_.data(), nR_, lines_, policy);
  }

  /// Back to radius-major, into a stack of the right shape.
  void CopyInto(Stack& stack,
                Execution policy = Execution::Sequential()) const {
    if (stack.NumberOfRadii() != nR_ || stack.SliceSize() != lines_) {
      throw std::invalid_argument(
          "A radial-major buffer can only be copied back into a stack of the "
          "shape it came from");
    }
    auto target = stack.Data();
    Transpose(data_.data(), target.data(), lines_, nR_, policy);
  }

 private:
  RadialGridType radial_;
  Int nR_;
  Int lines_;
  std::vector<Scalar> data_;

  RadialMajor(RadialGridType radial, Int nR, Int lines)
      : radial_{std::move(radial)},
        nR_{nR},
        lines_{lines},
        data_(static_cast<std::size_t>(nR) * static_cast<std::size_t>(lines)) {}

  std::size_t Offset(Int j) const {
    if (j < 0 || j >= lines_) {
      throw std::invalid_argument("Line index outside the slice");
    }
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(nR_);
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
    const auto threads = policy.TeamSize();

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
// is no gather and no scatter, because the line is already contiguous.
//
// `in` and `out` may be the same buffer, and in the loop this layout exists
// for -- one operator after another on the same lines -- they usually are.
// The operator is then given a scratch line to write into and the result is
// copied over the input afterwards, so it never sees its output alias its
// input and need not be written to cope with that. Two different buffers
// cannot overlap, each owning its storage, so they are handed over directly
// and cost no copy.
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

  RadialDetails::CheckOperatorGrid(op, in.Radial());

  const auto lines = in.NumberOfLines();
  const auto nR = static_cast<std::size_t>(in.NumberOfRadii());
  const auto threads = policy.TeamSize();
  const bool inPlace = static_cast<const void*>(&in) == &out;

  const auto run = [&](Int j) {
    if (!inPlace) {
      op(std::span<const Scalar>(in.Line(j)), std::span<Scalar>(out.Line(j)));
      return;
    }
    thread_local auto scratch = std::vector<Scalar>{};
    if (scratch.size() < nR) scratch.resize(nR);
    op(std::span<const Scalar>(in.Line(j)),
       std::span<Scalar>(scratch.data(), nR));
    std::copy_n(scratch.data(), nR, out.Line(j).data());
  };

  if (threads == 1) {
    for (auto j = Int{0}; j < lines; j++) run(j);
  } else {
    // Nothing may leave the region, and what runs in it is the caller's: see
    // ExceptionCapture. The first line to throw stops the rest being started,
    // and its exception is the one the caller sees, as it would be
    // sequentially.
    auto capture = Details::ExceptionCapture{};
#pragma omp parallel for schedule(static) num_threads(threads)
    for (Int j = 0; j < lines; j++) {
      if (capture.Failed()) continue;
      capture.Run([&] { run(j); });
    }
    capture.Rethrow();
  }
}

//--------------------------------------------------------------------------//
//                  Transforming straight into radial lines                  //
//--------------------------------------------------------------------------//

// Expand a layered field into radial lines, and evaluate radial lines into a
// layered field, without the radius-major expansion ever existing.
//
// `RadialMajor(Expand(field))` makes the expansion and then a transposed copy
// of it. These make only the copy: the transform is told that its coefficient
// side is RadialMajor::Batch() and writes the lines itself. The numbers are
// the same to the last bit, being one transform writing to a different place.
//
// **What this saves is memory, and it should be chosen for that.** One whole
// set of coefficients: `nR * CoefficientSize * 16` bytes a scalar field, which
// is 211 MB at lMax = 256 with 200 radii and 2.1 GB at lMax = 512 with 500.
//
// It does not reliably save time, and in one case it costs it. The transpose
// is five to ten per cent of a transform, so that is the most there is to
// win, and measured on a laptop it is won sequentially and with the matrix
// kernel under threads, by fifteen to twenty-five per cent at lMax = 128 --
// *but not always when nR is a power of two*. The scatter into this layout
// has stride nR, and at lMax = 256 with nR = 64 and 128 under threads it ran
// eight to twenty per cent *slower* than transform-and-transpose: successive
// writes land in the same cache sets, which is the hazard the Fourier stage
// guards against and the reason the transpose above is tiled. With the loop
// kernel the two routes are within a few per cent, and at nR = 128 the direct
// one was the slower by up to a sixth. So there is no rule
// here simple enough to build in, and none is: where time matters more than
// memory, measure both on the machine in question -- the benchmark's `lines`
// section is that measurement.

/**
 * @brief The expansion of a layered field, as radial lines.
 * @param field The field, radius-major as every layered field is.
 * @param lMax The degree to expand to.
 * @param policy Whether the transform may thread.
 * @return What `RadialMajor(Expand(field, lMax))` holds, bit for bit.
 */
template <std::ptrdiff_t N, AngularGrid Grid, RealOrComplexValued Value>
auto ExpandToLines(const LayeredSpinField<N, Grid, Value>& field,
                   std::ptrdiff_t lMax,
                   Execution policy = Execution::Sequential()) {
  using Lines = RadialMajor<LayeredSpinExpansion<N, Grid, Value>>;
  const auto& grid = field.Grid();
  const auto size = std::same_as<Value, RealValued>
                        ? grid.RealCoefficientSize(lMax)
                        : grid.CoefficientSize(lMax, N);
  auto lines = Lines::OfShape(field.Radial(), size);
  auto out = lines.Data();
  grid.ForwardTransformation(lMax, N, field.Data(), field.Batch(), out,
                             lines.Batch(), policy);
  return lines;
}

/**
 * @brief The layered field whose expansion a set of radial lines is.
 * @details A radial-major buffer is a shape and knows neither its angular grid
 * nor its degree, so both are given. They are checked against the shape,
 * which is all that can be checked.
 * @param lines The coefficients, one line per (l, m).
 * @param grid The angular grid to evaluate on.
 * @param lMax The degree the lines were expanded to.
 * @param policy Whether the transform may thread.
 * @throws std::invalid_argument if the lines are not as many as that degree
 * has coefficients.
 */
template <std::ptrdiff_t N, AngularGrid Grid, RealOrComplexValued Value>
auto EvaluateLines(
    const RadialMajor<LayeredSpinExpansion<N, Grid, Value>>& lines,
    const Grid& grid, std::ptrdiff_t lMax,
    Execution policy = Execution::Sequential()) {
  const auto size = std::same_as<Value, RealValued>
                        ? grid.RealCoefficientSize(lMax)
                        : grid.CoefficientSize(lMax, N);
  if (lines.NumberOfLines() != size) {
    throw std::invalid_argument(
        "These radial lines number " + std::to_string(lines.NumberOfLines()) +
        ", and an expansion to degree " + std::to_string(lMax) + " has " +
        std::to_string(size) +
        " coefficients, so they were not expanded to that degree");
  }
  auto field = LayeredSpinField<N, Grid, Value>(lines.Radial(), grid);
  auto out = field.Data();
  grid.InverseTransformation(lMax, N, lines.Data(), lines.Batch(), out,
                             field.Batch(), policy);
  return field;
}

}  // namespace GSHTrans

#ifndef _OPENMP
#pragma GCC diagnostic pop
#endif
