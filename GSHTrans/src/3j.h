#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "NumericConcepts/Numeric.hpp"
#include "NumericConcepts/Ranges.hpp"

namespace GSHTrans {

/**
 * @file Wigner3j.hpp
 * @brief Wigner 3j symbols by stable recursion over the (m1, m3) plane.
 *
 * @details For fixed degrees (l1, l2, l3) the full table of symbols
 *
 *           ( l1  l2  l3 )
 *           ( m1  m2  m3 ),   with m2 = -(m1 + m3),
 *
 * is generated over the axes m1 in [-l1, l1] and m3 in [-l3, l3], one row at
 * a time, by the Schulten-Gordon algorithm.
 *
 * The order recurrence has a growing and a decaying solution, with a
 * classically allowed region between two forbidden ones. Recursing inward
 * from a forbidden end follows the growing solution and is stable; outward
 * follows the decaying one and loses the answer exponentially. So each row is
 * built inward from *both* ends, the two halves are matched where they
 * overlap, and the one remaining constant is fixed from the unitary property
 *
 *     (2 l1 + 1) sum_{m2} g(m2)^2 = 1.
 *
 * No closed-form seed is needed and no factorial is ever formed, so the
 * method holds at large degree: measured against an independent route
 * (column-permutation invariance, which runs the recursion along different
 * lines) it agrees to 1e-16 at (200,200,200) and 1e-15 at (1000,1000,1999).
 *
 * This replaced two earlier schemes, and docs/3j-plan.md records why. A port
 * of Woodhouse's wig2.f ran the three-term recursion in one direction only,
 * and failed exponentially near stretched triangles -- unusable past l = 30
 * there. Racah's closed form was added as a fallback and covers exactly the
 * region that one cannot, but the two together still left a band at
 * intermediate shapes above l = 80 that neither reached. Schulten-Gordon
 * covers all of it, at 1.1x to 1.3x the cost of the scheme it replaces.
 *
 * A second layout is offered beside the plain table -- see CouplingElement()
 * and FillCouplingMatrix -- in which the first order is negated and carries an
 * alternating phase. That is the arrangement in which the symbols appear in
 * normal-mode coupling matrices, and it is what Woodhouse's wig2.f tabulated,
 * so codes written against that routine want it. It is a phase and a
 * relabelling of the table below, not a separate algorithm.
 *
 * All degrees and orders are integers (no half-integer support).
 */

/*------------------------------------------------------------------------*/
/*                          Concepts and helpers                          */
/*------------------------------------------------------------------------*/

/**
 * @brief Concept for a contiguous, sized range of real numbers that can be
 * read and written. Used for the in-place kernels that fill user-provided
 * storage.
 * @tparam R The range type to check.
 */
template <typename R>
concept RealContiguousWritableRange =
    NumericConcepts::RealRange<R> and NumericConcepts::RealWritableRange<R> and
    std::ranges::contiguous_range<R> and std::ranges::sized_range<R>;

/**
 * @brief Checks the triangle inequality |l1 - l3| <= l2 <= l1 + l3.
 * @param l1 The first degree.
 * @param l2 The second degree.
 * @param l3 The third degree.
 * @return true if (l1, l2, l3) can form a triangle.
 */
constexpr bool SatisfiesTriangle(int l1, int l2, int l3) {
  const auto ldel = l1 > l3 ? l1 - l3 : l3 - l1;
  return l2 >= ldel and l2 <= l1 + l3;
}

/**
 * @brief A closed integer interval [Min, Max] used to label an index axis
 * of a table of Wigner symbols (e.g. an order axis m in [-l, l]).
 * @details An Axis is itself a range over its values, so it can be used
 * directly in range-based for loops:
 * @code
 *   for (auto m1 : table.M1Axis()) { ... }
 * @endcode
 */
class Axis {
 public:
  /** @brief Constructs an empty axis. */
  constexpr Axis() = default;

  /**
   * @brief Constructs the axis [min, max].
   * @param min The smallest index on the axis.
   * @param max The largest index on the axis.
   */
  constexpr Axis(int min, int max) : _min{min}, _max{max} {}

  /** @brief Returns the smallest index on the axis. */
  constexpr auto Min() const { return _min; }

  /** @brief Returns the largest index on the axis. */
  constexpr auto Max() const { return _max; }

  /** @brief Returns the number of indices on the axis. */
  constexpr auto Size() const {
    return _max >= _min ? static_cast<std::size_t>(_max - _min + 1)
                        : std::size_t{0};
  }

  /** @brief Returns true if the index i lies on the axis. */
  constexpr auto Contains(int i) const { return i >= _min and i <= _max; }

  /** @brief Returns the zero-based offset of the index i along the axis. */
  constexpr auto Offset(int i) const {
    assert(Contains(i));
    return static_cast<std::size_t>(i - _min);
  }

  /** @brief Returns a view of the indices on the axis. */
  constexpr auto Values() const {
    return std::views::iota(_min, std::max(_min, _max + 1));
  }

  /** @brief Iterator to the first index on the axis. */
  constexpr auto begin() const { return Values().begin(); }

  /** @brief Iterator past the last index on the axis. */
  constexpr auto end() const { return Values().end(); }

  /** @brief Axes compare equal when their bounds agree. */
  constexpr bool operator==(const Axis&) const = default;

 private:
  int _min{0};
  int _max{-1};
};

/*------------------------------------------------------------------------*/
/*                        The computational kernels                        */
/*------------------------------------------------------------------------*/

namespace ThreeJDetails {

/**
 * @brief The two coefficients of the order recurrence.
 * @details With g(m2) = (l1 l2 l3; m1 m2 -m1-m2), the symbols satisfy
 *
 *     A(m2) g(m2) + B(m2) g(m2-1) + A(m2-1) g(m2-2) = 0
 *
 * with A and B as below. A(m2min) is identically zero -- that is the lower
 * endpoint of the run -- which is why the backward pass takes its coefficient
 * one step up from the index it is filling and never forms it.
 */
template <NumericConcepts::Real T>
T RecurrenceA(int l1, int l2, int l3, int m1, int m2) {
  const auto m3 = -m1 - m2;
  return std::sqrt(static_cast<T>(l2 - m2 + 1) * static_cast<T>(l2 + m2) *
                   static_cast<T>(l3 + m3 + 1) * static_cast<T>(l3 - m3));
}

template <NumericConcepts::Real T>
T RecurrenceB(int l1, int l2, int l3, int m1, int m2) {
  const auto m3 = -m1 - m2;
  return static_cast<T>(l1 + l2 + l3 + 1) * static_cast<T>(l2 + l3 - l1) -
         static_cast<T>(l2 - m2 + 1) * static_cast<T>(l3 + m3 + 1) -
         static_cast<T>(l2 + m2 - 1) * static_cast<T>(l3 - m3 - 1);
}

/** @brief The orders m2 for which g(m2) can be non-zero. */
inline std::pair<int, int> OrderRange(int l2, int l3, int m1) {
  return {std::max(-l2, -l3 - m1), std::min(l2, l3 - m1)};
}

/**
 * @brief One row of the (m1, m3) plane by the Schulten-Gordon algorithm.
 * @details Fills g[0 .. n-1] with g(m2) for m2 from m2min to m2max at fixed
 * m1, where n = m2max - m2min + 1.
 *
 * The recurrence has a growing and a decaying solution, and a classically
 * allowed region between two forbidden ones. Recursing *inward* from a
 * forbidden end follows the growing solution, so errors decay; outward
 * follows the decaying one, so they grow exponentially -- which is what
 * destroyed the one-directional scheme this replaces, and which extra
 * precision delays rather than cures.
 *
 * So the run comes inward from both ends and the two halves are matched where
 * they overlap. Each is correct up to its own scale; matching leaves one
 * unknown constant, and the unitary property
 *
 *     (2 l1 + 1) sum_{m2} g(m2)^2 = 1
 *
 * fixes it. No closed-form seed is needed and no factorial is ever formed.
 *
 * The turning point is found by watching the recurrence coefficient rather
 * than by locating it analytically: while |c1| is decreasing the run is
 * heading towards larger values and is stable, and the first increase says to
 * stop and come from the other end.
 *
 * Schulten, K. and Gordon, R. G., J. Math. Phys. 16 (1975) 1961, and the
 * companion at 1971 for the semiclassical picture behind the turning points.
 * The control flow was checked against SLATEC's DRC3JM; docs/3j-plan.md T6
 * records what was changed and why.
 */
template <NumericConcepts::Real T>
int SchultenGordonRow(int l1, int l2, int l3, int m1, std::span<T> g) {
  const auto [m2Min, m2Max] = OrderRange(l2, l3, m1);
  const auto n = m2Max - m2Min + 1;
  if (n < 1) return 0;

  std::fill_n(g.begin(), n, T{0});

  // One order: the recurrence has nowhere to run, and the value is the
  // convention's phase over the square root of the perimeter.
  if (n == 1) {
    g[0] = (std::abs(l2 - l3 - m1) % 2 == 0 ? T{1} : T{-1}) /
           std::sqrt(static_cast<T>(l1 + l2 + l3 + 1));
    return n;
  }

  const auto big = std::sqrt(std::numeric_limits<T>::max() / 20);
  const auto rootBig = std::sqrt(big);
  const auto small = 1 / big;
  const auto rootSmall = 1 / rootBig;

  const auto A = [&](int m2) { return RecurrenceA<T>(l1, l2, l3, m1, m2); };
  const auto B = [&](int m2) { return RecurrenceB<T>(l1, l2, l3, m1, m2); };

  // Forward from m2Min, while |c1| decreases.
  g[0] = rootSmall;
  auto sumAll = small;
  auto sumForward = small;
  auto factor = T{0};
  auto c1 = T{0};
  auto previous = std::numeric_limits<T>::max();
  auto value = T{0};
  auto last = 0;

  for (auto k = 1; k < n; ++k) {
    const auto m2 = m2Min + k;
    const auto older = factor;
    factor = A(m2);
    if (k > 1) previous = std::abs(c1);
    c1 = -B(m2) / factor;

    if (k == 1) {
      value = rootSmall * c1;  // the third term vanishes at the first step
    } else {
      value = c1 * g[k - 1] - (older / factor) * g[k - 2];
    }
    g[k] = value;
    sumForward = sumAll;
    sumAll += value * value;
    last = k;
    if (k == n - 1) break;

    if (std::abs(value) >= rootBig) {
      for (auto i = 0; i <= k; ++i) {
        if (std::abs(g[i]) < rootSmall) g[i] = 0;
        g[i] /= rootBig;
      }
      sumAll /= big;
      sumForward /= big;
      value /= rootBig;
    }
    if (previous <= std::abs(c1)) break;
  }

  auto sumTotal = sumAll;

  // The sign the tail carries, tracked rather than read back from storage.
  //
  // The convention fixes the sign of g(m2Max), and SLATEC recovers it from
  // the computed array. That is not safe here: a near-stretched row at high
  // degree spans a dynamic range of 1e201, so the rescaling above flushes the
  // tail to zero and the sign with it, and whole rows come out negated with
  // every magnitude correct to rounding. Neither the completeness relation
  // nor the recurrence residual can see that ([J8]).
  auto tailSign = (g[n - 1] >= 0) ? T{1} : T{-1};

  if (last < n - 1) {
    // Backward from m2Max, overlapping the forward run at three points.
    const auto x1 = g[last];
    const auto x2 = g[last - 1];
    const auto x3 = g[last - 2];

    g[n - 1] = rootSmall;
    auto sumBack = small;
    auto sumBackward = small;
    factor = 0;
    auto y1 = T{0}, y2 = T{0}, y3 = T{0};

    for (auto j = n - 2; j >= last - 2; --j) {
      const auto older = factor;
      factor = A(m2Min + j + 1);
      c1 = -B(m2Min + j + 2) / factor;

      const auto y = (j == n - 2)
                         ? rootSmall * c1  // the third term vanishes
                         : c1 * g[j + 1] - (older / factor) * g[j + 2];

      if (j == last - 2) {  // the match point: compare, do not store
        y3 = y;
        y2 = g[j + 1];
        y1 = g[j + 2];
        break;
      }

      g[j] = y;
      sumBackward = sumBack;
      sumBack += y * y;

      if (std::abs(y) >= rootBig) {
        for (auto i = j; i < n; ++i) {
          if (std::abs(g[i]) < rootSmall) g[i] = 0;
          g[i] /= rootBig;
        }
        sumBack /= big;
        sumBackward /= big;
      }
    }

    // Least squares over the three overlapping points, which is steadier than
    // matching on one of them.
    auto ratio = (x1 * y1 + x2 * y2 + x3 * y3) /
                 (x1 * x1 + x2 * x2 + x3 * x3);

    if (std::abs(ratio) >= 1) {
      for (auto i = 0; i <= last - 2; ++i) g[i] *= ratio;
      sumTotal = ratio * ratio * sumForward + sumBackward;
      tailSign = 1;  // the backward seed was positive and is untouched
    } else {
      ratio = 1 / ratio;
      for (auto i = last - 1; i < n; ++i) g[i] *= ratio;
      sumTotal = sumForward + ratio * ratio * sumBackward;
      tailSign = (ratio >= 0) ? T{1} : T{-1};
    }
  }

  auto norm = 1 / std::sqrt(static_cast<T>(2 * l1 + 1) * sumTotal);
  const auto wanted = (std::abs(l2 - l3 - m1) % 2 == 0) ? T{1} : T{-1};
  if (tailSign * wanted < 0) norm = -norm;
  for (auto i = 0; i < n; ++i) g[i] *= norm;
  return n;
}

/**
 * @brief The tolerance the recurrence residual is judged against.
 * @details Relative to the row's own largest value, and loose: the residual
 * of a good row is a few epsilon times the number of steps, while a broken
 * match shows up at order one.
 */
template <NumericConcepts::Real T>
T ResidualTolerance(int steps) {
  return static_cast<T>(1000 * (steps + 1)) *
         std::numeric_limits<T>::epsilon();
}

/**
 * @brief Checks a completed row against the recurrence that defines it.
 * @details This is the runtime self-check, and it replaced the completeness
 * relation when the algorithm changed ([J6]). Completeness was the right
 * check for a one-directional recursion seeded from a closed form; it is
 * nearly worthless against Schulten-Gordon, which normalises every row by
 * that very identity, so the sum is one by construction whatever the row
 * holds -- and in particular a mis-scaled join, which is this algorithm's
 * characteristic failure, passes it.
 *
 * The defining recurrence does not pass it. A bad match violates the relation
 * at the join, a bad rescale violates it locally, and a NaN propagates.
 *
 * What it cannot see, since the recurrence is homogeneous: an overall factor
 * or sign on the row. Scale is fixed by the normalisation and sign by the
 * convention, and both are properties of the code rather than of the input,
 * so both are pinned by tests instead.
 */
template <NumericConcepts::Real T>
bool RowSatisfiesRecurrence(int l1, int l2, int l3, int m1,
                            std::span<const T> g, int n) {
  if (n < 3) return true;
  const auto [m2Min, m2Max] = OrderRange(l2, l3, m1);
  (void)m2Max;

  auto scale = T{0};
  for (auto k = 0; k < n; ++k) scale = std::max(scale, std::abs(g[k]));
  if (not(scale > 0)) return false;

  const auto tolerance = ResidualTolerance<T>(n) * scale;
  for (auto k = 2; k < n; ++k) {
    const auto m2 = m2Min + k;
    const auto residual = RecurrenceA<T>(l1, l2, l3, m1, m2) * g[k] +
                          RecurrenceB<T>(l1, l2, l3, m1, m2) * g[k - 1] +
                          RecurrenceA<T>(l1, l2, l3, m1, m2 - 1) * g[k - 2];
    // The coefficients are of order l^2, so the residual is measured against
    // the scale they multiply rather than against the values alone.
    const auto local = std::max(
        {std::abs(RecurrenceA<T>(l1, l2, l3, m1, m2)),
         std::abs(RecurrenceB<T>(l1, l2, l3, m1, m2)),
         std::abs(RecurrenceA<T>(l1, l2, l3, m1, m2 - 1))});
    if (not(std::abs(residual) <= tolerance * local)) return false;
  }
  return true;
}

/**
 * @brief Converts between the plain (m1, m3) table and the coupling layout.
 * @details The coupling layout is
 *
 *     c(m, mp) = (-1)^m ( l1     l2      l3 )
 *                       ( -m   m - mp    mp ),
 *
 * so plain(m1, m3) = (-1)^m1 c(-m1, m3): the rows are mirrored in m1 and
 * carry an alternating phase. It is the array Woodhouse's wig2.f returned.
 *
 * **It is its own inverse**, since mirroring twice is the identity and the
 * phase squares to one, so one routine serves both directions and there is no
 * second convention to keep in step.
 */
template <NumericConcepts::Real T>
void SwapCouplingConvention(const int l1, const int l3, std::span<T> a) {
  const auto columns = static_cast<std::size_t>(2 * l3 + 1);
  for (auto r = 0; r < l1; ++r) {
    const auto phase = ((r + l1) % 2 == 0) ? T{1} : T{-1};
    auto* lower = a.data() + static_cast<std::size_t>(r) * columns;
    auto* upper = a.data() + static_cast<std::size_t>(2 * l1 - r) * columns;
    for (std::size_t j = 0; j < columns; ++j) {
      const auto tmp = lower[j];
      lower[j] = phase * upper[j];
      upper[j] = phase * tmp;
    }
  }
}

/**
 * @brief Fills the whole (m1, m3) plane, row by row.
 * @details One Schulten-Gordon row recursion per m1, scattered into the
 * row-major (m1, m3) table. Each row is normalised independently, which is
 * stronger than one global normalisation would be: an error in one row cannot
 * leak into another.
 */
template <NumericConcepts::Real T>
void Wigner3jPlane(int l1, int l2, int l3, std::span<T> table) {
  std::fill(table.begin(), table.end(), T{0});
  if (not SatisfiesTriangle(l1, l2, l3)) return;

  const auto columns = 2 * l3 + 1;
  auto row = std::vector<T>(static_cast<std::size_t>(2 * l2 + 2));

  for (auto m1 = -l1; m1 <= l1; ++m1) {
    const auto [m2Min, m2Max] = OrderRange(l2, l3, m1);
    if (m2Max < m2Min) continue;
    const auto n = SchultenGordonRow<T>(l1, l2, l3, m1, std::span<T>(row));

    if (not RowSatisfiesRecurrence<T>(l1, l2, l3, m1,
                                      std::span<const T>(row.data(), n), n)) {
      throw std::runtime_error(
          "Wigner3jMatrix: the recurrence is not satisfied for degrees (" +
          std::to_string(l1) + ", " + std::to_string(l2) + ", " +
          std::to_string(l3) + ") at m1 = " + std::to_string(m1) +
          ". The values are not to be trusted. See docs/3j-plan.md.");
    }

    for (auto k = 0; k < n; ++k) {
      const auto m3 = -m1 - (m2Min + k);
      if (std::abs(m3) > l3) continue;
      table[static_cast<std::size_t>(m1 + l1) * columns + (m3 + l3)] = row[k];
    }
  }
}

}  // namespace ThreeJDetails

/*------------------------------------------------------------------------*/
/*                         Free kernel functions                           */
/*------------------------------------------------------------------------*/

/**
 * @brief Fills user-provided storage with the table of Wigner 3j symbols
 * (l1 l2 l3; m1, -(m1+m3), m3) in row-major order, m1 being the row axis
 * and m3 the column axis. Allocation-free.
 * @tparam Range A contiguous, sized, readable and writable range of reals.
 * @param l1 The first degree.
 * @param l2 The second degree.
 * @param l3 The third degree.
 * @param table Storage of size (2 l1 + 1) x (2 l3 + 1).
 */
template <RealContiguousWritableRange Range>
void FillWigner3jMatrix(int l1, int l2, int l3, Range&& table) {
  using T = std::ranges::range_value_t<Range>;
  auto span = std::span<T>(std::ranges::data(table), std::ranges::size(table));
  ThreeJDetails::Wigner3jPlane<T>(l1, l2, l3, span);
}

/**
 * @brief Fills user-provided storage with the coupling layout
 * (-1)^m (l1 l2 l3; -m, m-mp, mp) in row-major order, m being the row axis
 * and mp the column axis: the first order negated, with an alternating phase.
 * This is the arrangement the symbols appear in within normal-mode coupling
 * matrices, and it is exactly the array a(m+l1+1, mp+l3+1) returned by the
 * Fortran routine wig2(l1, l2, l3, a, id1). Allocation-free.
 * @tparam Range A contiguous, sized, readable and writable range of reals.
 * @param l1 The first degree.
 * @param l2 The second degree.
 * @param l3 The third degree.
 * @param table Storage of size (2 l1 + 1) x (2 l3 + 1).
 */
template <RealContiguousWritableRange Range>
void FillCouplingMatrix(int l1, int l2, int l3, Range&& table) {
  using T = std::ranges::range_value_t<Range>;
  auto span = std::span<T>(std::ranges::data(table), std::ranges::size(table));
  ThreeJDetails::Wigner3jPlane<T>(l1, l2, l3, span);
  ThreeJDetails::SwapCouplingConvention<T>(l1, l3, span);
}

/*------------------------------------------------------------------------*/
/*                            Wigner3jMatrix                               */
/*------------------------------------------------------------------------*/

/**
 * @brief Table of Wigner 3j symbols for fixed degrees (l1, l2, l3) over
 * the order axes m1 in [-l1, l1] and m3 in [-l3, l3], with the middle
 * order fixed by the selection rule m2 = -(m1 + m3).
 * @details Storage is a single contiguous std::vector in row-major order
 * (m1 slowest, m3 fastest). Element access is by quantum numbers rather
 * than raw offsets; queries outside the axes return zero, which is the
 * mathematically consistent value of the symbol. The class models
 * std::ranges::random_access_range over its values, and therefore also
 * satisfies NumericConcepts::RealRange.
 * @tparam T The floating-point type used for computation and storage.
 */
template <NumericConcepts::Real T>
class Wigner3jMatrix {
 public:
  using value_type = T;

  /** @brief Default constructor: the trivial table (0 0 0; 0 0 0) = 1. */
  Wigner3jMatrix() : Wigner3jMatrix(0, 0, 0) {}

  /**
   * @brief Computes the table for degrees (l1, l2, l3).
   * @param l1 The first degree.
   * @param l2 The second degree.
   * @param l3 The third degree.
   */
  Wigner3jMatrix(int l1, int l2, int l3)
      : _l1{l1},
        _l2{l2},
        _l3{l3},
        _data(static_cast<std::size_t>(2 * l1 + 1) *
              static_cast<std::size_t>(2 * l3 + 1)) {
    assert(l1 >= 0 and l2 >= 0 and l3 >= 0);
    ThreeJDetails::Wigner3jPlane<T>(_l1, _l2, _l3, std::span<T>(_data));
  }

  /** @brief Returns the first degree. */
  auto L1() const { return _l1; }

  /** @brief Returns the second degree. */
  auto L2() const { return _l2; }

  /** @brief Returns the third degree. */
  auto L3() const { return _l3; }

  /** @brief Returns the degrees as an array {l1, l2, l3}. */
  auto Degrees() const { return std::array{_l1, _l2, _l3}; }

  /** @brief Returns the axis for the first order, m1 in [-l1, l1]. */
  auto M1Axis() const { return Axis(-_l1, _l1); }

  /** @brief Returns the axis for the third order, m3 in [-l3, l3]. */
  auto M3Axis() const { return Axis(-_l3, _l3); }

  /**
   * @brief Returns the symbol (l1 l2 l3; m1, -(m1+m3), m3).
   * @param m1 The first order.
   * @param m3 The third order.
   * @return The value of the symbol; zero if an order is out of range.
   */
  auto operator()(int m1, int m3) const {
    if (not(M1Axis().Contains(m1) and M3Axis().Contains(m3))) return T{0};
    return _data[M1Axis().Offset(m1) * M3Axis().Size() + M3Axis().Offset(m3)];
  }

  /**
   * @brief Returns the symbol (l1 l2 l3; m1 m2 m3) with all three orders
   * given explicitly.
   * @param m1 The first order.
   * @param m2 The second order.
   * @param m3 The third order.
   * @return The value of the symbol; zero if the selection rules
   * m1 + m2 + m3 = 0 and |mi| <= li are not met.
   */
  auto operator()(int m1, int m2, int m3) const {
    if (m1 + m2 + m3 != 0 or std::abs(m2) > _l2) return T{0};
    return (*this)(m1, m3);
  }

  /**
   * @brief Returns the symbol in the coupling layout,
   * (-1)^m (l1 l2 l3; -m, m-mp, mp): the first order negated, with an
   * alternating phase.
   * @details This is how the symbols enter a normal-mode coupling matrix, and
   * it is the quantity the Fortran routine wig2.f tabulated, so it is offered
   * for codes written against that convention. It is a relabelling of the
   * plain table rather than a different calculation.
   * @param m The row order, m in [-l1, l1].
   * @param mp The column order, mp in [-l3, l3].
   */
  auto CouplingElement(int m, int mp) const {
    const auto phase = (m % 2 == 0) ? T{1} : T{-1};
    return phase * (*this)(-m, mp);
  }

  /** @brief Returns a view of the row of symbols with fixed m1. */
  auto Row(int m1) const {
    return std::span<const T>(
        _data.data() + M1Axis().Offset(m1) * M3Axis().Size(), M3Axis().Size());
  }

  /** @brief Returns a flat view of the table, row-major in (m1, m3). */
  auto Data() const { return std::span<const T>(_data); }

  /** @brief Iterator to the start of the flattened table. */
  auto begin() const { return _data.cbegin(); }

  /** @brief Iterator past the end of the flattened table. */
  auto end() const { return _data.cend(); }

  /** @brief Returns the total number of stored values. */
  auto size() const { return _data.size(); }

 private:
  int _l1{0};
  int _l2{0};
  int _l3{0};
  std::vector<T> _data;
};

/*------------------------------------------------------------------------*/
/*                             Wigner3jStack                               */
/*------------------------------------------------------------------------*/

/**
 * @brief Table of Wigner 3j symbols for fixed outer degrees (l1, l3) over
 * three index axes: the middle degree l2, and the orders m1 and m3, with
 * m2 = -(m1 + m3) fixed by the selection rule.
 * @details Stored as one Wigner3jMatrix per value of l2. By default l2
 * runs over the full triangle range [|l1 - l3|, l1 + l3], but any
 * non-negative range may be requested; matrices outside the triangle
 * range are identically zero. The middle degree was chosen as the degree
 * axis because that is the common pattern in applications (e.g. coupling
 * through a structural degree l2 in normal-mode seismology); other
 * arrangements follow from the column-permutation symmetries of the
 * symbols. The stack is itself a range over its matrices.
 * @tparam T The floating-point type used for computation and storage.
 */
template <NumericConcepts::Real T>
class Wigner3jStack {
 public:
  using value_type = T;

  /**
   * @brief Computes tables for all l2 in the triangle range
   * [|l1 - l3|, l1 + l3].
   * @param l1 The first degree.
   * @param l3 The third degree.
   */
  Wigner3jStack(int l1, int l3)
      : Wigner3jStack(l1, l3, std::abs(l1 - l3), l1 + l3) {}

  /**
   * @brief Computes tables for all l2 in [l2Min, l2Max].
   * @param l1 The first degree.
   * @param l3 The third degree.
   * @param l2Min The smallest middle degree.
   * @param l2Max The largest middle degree.
   */
  Wigner3jStack(int l1, int l3, int l2Min, int l2Max)
      : _l1{l1}, _l3{l3}, _l2Axis{l2Min, l2Max} {
    assert(l1 >= 0 and l3 >= 0 and l2Min >= 0 and l2Max >= l2Min);
    _matrices.reserve(_l2Axis.Size());
    for (auto l2 : _l2Axis) _matrices.emplace_back(_l1, l2, _l3);
  }

  /** @brief Returns the first degree. */
  auto L1() const { return _l1; }

  /** @brief Returns the third degree. */
  auto L3() const { return _l3; }

  /** @brief Returns the axis for the middle degree l2. */
  auto L2Axis() const { return _l2Axis; }

  /** @brief Returns the axis for the first order, m1 in [-l1, l1]. */
  auto M1Axis() const { return Axis(-_l1, _l1); }

  /** @brief Returns the axis for the third order, m3 in [-l3, l3]. */
  auto M3Axis() const { return Axis(-_l3, _l3); }

  /**
   * @brief Returns the table of symbols for a given middle degree.
   * @param l2 The middle degree; must lie on L2Axis().
   */
  const Wigner3jMatrix<T>& Matrix(int l2) const {
    assert(_l2Axis.Contains(l2));
    return _matrices[_l2Axis.Offset(l2)];
  }

  /** @brief Equivalent to Matrix(l2), allowing stack[l2](m1, m3). */
  const Wigner3jMatrix<T>& operator[](int l2) const { return Matrix(l2); }

  /**
   * @brief Returns the symbol (l1 l2 l3; m1, -(m1+m3), m3).
   * @param l2 The middle degree.
   * @param m1 The first order.
   * @param m3 The third order.
   * @return The value of the symbol; zero if l2 is off the stored axis or
   * an order is out of range.
   */
  auto operator()(int l2, int m1, int m3) const {
    if (not _l2Axis.Contains(l2)) return T{0};
    return Matrix(l2)(m1, m3);
  }

  /**
   * @brief Returns the symbol (l1 l2 l3; m1 m2 m3) with all three orders
   * given explicitly.
   * @param l2 The middle degree.
   * @param m1 The first order.
   * @param m2 The second order.
   * @param m3 The third order.
   */
  auto operator()(int l2, int m1, int m2, int m3) const {
    if (not _l2Axis.Contains(l2)) return T{0};
    return Matrix(l2)(m1, m2, m3);
  }

  /**
   * @brief Returns the symbol in the coupling layout
   * (-1)^m (l1 l2 l3; -m, m-mp, mp) for a given middle degree.
   * @param l2 The middle degree.
   * @param m The row order.
   * @param mp The column order.
   */
  auto CouplingElement(int l2, int m, int mp) const {
    if (not _l2Axis.Contains(l2)) return T{0};
    return Matrix(l2).CouplingElement(m, mp);
  }

  /** @brief Iterator to the first stored matrix (smallest l2). */
  auto begin() const { return _matrices.cbegin(); }

  /** @brief Iterator past the last stored matrix. */
  auto end() const { return _matrices.cend(); }

  /** @brief Returns the number of stored matrices. */
  auto size() const { return _matrices.size(); }

 private:
  int _l1{0};
  int _l3{0};
  Axis _l2Axis;
  std::vector<Wigner3jMatrix<T>> _matrices;
};

/*------------------------------------------------------------------------*/
/*                          Single-symbol helper                           */
/*------------------------------------------------------------------------*/

/**
 * @brief Convenience function returning a single symbol
 * (l1 l2 l3; m1 m2 m3).
 * @details Builds the full (m1, m3) table internally, at a cost of
 * O((2 l1 + 1)(2 l3 + 1)); when many symbols with shared degrees are
 * needed, construct a Wigner3jMatrix or Wigner3jStack once and query it
 * instead.
 * @tparam T The floating-point type; defaults to double.
 */
template <NumericConcepts::Real T = double>
T Wigner3jSymbol(int l1, int l2, int l3, int m1, int m2, int m3) {
  if (m1 + m2 + m3 != 0) return T{0};
  if (std::abs(m1) > l1 or std::abs(m2) > l2 or std::abs(m3) > l3) return T{0};
  if (not SatisfiesTriangle(l1, l2, l3)) return T{0};
  return Wigner3jMatrix<T>(l1, l2, l3)(m1, m3);
}

/*------------------------------------------------------------------------*/

static_assert(NumericConcepts::RealRange<Wigner3jMatrix<double>>);
static_assert(NumericConcepts::RealRange<Wigner3jMatrix<float>>);

}  // namespace GSHTrans
