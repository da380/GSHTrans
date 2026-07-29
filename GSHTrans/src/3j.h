#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include "NumericConcepts/Numeric.hpp"
#include "NumericConcepts/Ranges.hpp"

namespace GSHTrans {

/**
 * @file Wigner3j.hpp
 * @brief Wigner 3j symbols by stable recursion over the (m1, m3) plane.
 *
 * @details A C++20 port of the Fortran routine wig2.f by J. H. Woodhouse.
 * For fixed degrees (l1, l2, l3) the full table of symbols
 *
 *           ( l1  l2  l3 )
 *           ( m1  m2  m3 ),   with m2 = -(m1 + m3),
 *
 * is generated over the axes m1 in [-l1, l1] and m3 in [-l3, l3] using
 * Woodhouse's recursion scheme:
 *
 *   1. The corner value at (m1, m3) = (l1, -l3) is evaluated in closed form.
 *   2. Two-term recursions propagate the value along the two boundary
 *      edges of the (m1, m3) plane.
 *   3. A three-term recursion sweeps diagonals of constant m2 from the
 *      corner towards the centre of the plane, but only over half of it.
 *   4. The remaining half is filled from the reflection symmetry of the
 *      symbols under (m1, m2, m3) -> (-m1, -m2, -m3).
 *
 * The scheme is numerically stable because the three-term recursion is
 * only ever used in the direction in which the symbols grow in magnitude
 * (from the classically forbidden corner region towards the centre); the
 * decaying tail on the far side is obtained by symmetry rather than by
 * continuing the recursion into an unstable regime. No factorials are
 * formed, so the method remains accurate for large degrees.
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

namespace Internal {

/**
 * @internal
 * @brief Direct port of wig2.f. Fills a row-major (2 l1 + 1) x (2 l3 + 1)
 * buffer with Woodhouse's matrix
 *
 *   a[(m + l1) * (2 l3 + 1) + (mp + l3)]
 *       = (-1)^m ( l1   l2    l3 )
 *                ( -m  m-mp   mp ),
 *
 * for all m in [-l1, l1] and mp in [-l3, l3]. Everything outside the
 * selection rules is set to zero, as is the whole buffer if the triangle
 * inequality fails.
 *
 * To ease comparison with the original, the local variable names follow
 * the Fortran source, in which the degrees (l1, l2, l3) here were called
 * (l1, is, l2), and the orders (m, mp) were called (m1, m2). Unlike the
 * original, integer products are converted to floating point factor by
 * factor, so no intermediate integer overflow can occur at large degrees.
 */
template <NumericConcepts::Real T>
void WoodhouseMatrix(const int l1, const int is, const int l2, std::span<T> a) {
  assert(l1 >= 0 and is >= 0 and l2 >= 0);
  const int nRows = 2 * l1 + 1;
  const int nCols = 2 * l2 + 1;
  assert(a.size() ==
         static_cast<std::size_t>(nRows) * static_cast<std::size_t>(nCols));
  const auto idx = [nCols](int i, int j) {
    return static_cast<std::size_t>(i) * static_cast<std::size_t>(nCols) +
           static_cast<std::size_t>(j);
  };

  std::ranges::fill(a, T{0});
  if (is < std::abs(l1 - l2) or is > l1 + l2) return;

  // Corner value at (m1, m2) = (-l1, -l2), in closed form.
  {
    auto r1 = T{1} / std::sqrt(T(2 * std::max(l1, l2) + 1));
    const auto ldel = std::abs(l1 - l2);
    for (int isc = ldel + 1; isc <= is; ++isc) {
      r1 *= std::sqrt(T(l1 + l2 - isc + 1) / T(l1 + l2 + isc + 1));
    }
    a[idx(0, 0)] = r1;
  }

  // First row: two-term recursion in m2 along the edge m1 = -l1.
  {
    const auto num1 = std::min(is - l1 + l2, l2);
    for (int n = 1; n <= num1; ++n) {
      const auto m2 = -l2 + n;
      a[idx(0, n)] =
          -a[idx(0, n - 1)] * std::sqrt(T(l1 + is + m2) * T(is - l1 - m2 + 1) /
                                        (T(l2 - m2 + 1) * T(l2 + m2)));
    }
  }

  // First column: two-term recursion in m1 along the edge m2 = -l2.
  {
    const auto num2 = is - l2 + l1;
    for (int n = 1; n <= num2; ++n) {
      const auto m1 = -l1 + n;
      a[idx(n, 0)] =
          -a[idx(n - 1, 0)] * std::sqrt(T(l2 + is + m1) * T(is - l2 - m1 + 1) /
                                        (T(l1 + m1) * T(l1 - m1 + 1)));
    }
  }

  // Three-term recursion along south-east diagonals, i.e. in the direction
  // (m1, m2) -> (m1 + 1, m2 + 1) of constant m1 - m2. Only the half-plane
  // m2 <= 0 is computed; the recursion therefore always runs towards
  // growing values and remains stable.
  const auto iss = static_cast<long long>(is) * (is + 1) -
                   static_cast<long long>(l1) * (l1 + 1) -
                   static_cast<long long>(l2) * (l2 + 1);
  const auto numd = std::min(2 * is + 1, is + l1);
  for (int nd = 1; nd <= numd; ++nd) {
    const auto m1b = std::max(-l1, is - l2 - nd + 1);
    const auto m2b = std::max(-l2, -is - l1 + nd - 1);
    const auto numit = std::min(-m2b, l1 - m1b);
    for (int nit = 1; nit <= numit; ++nit) {
      const auto m1 = m1b + nit;
      const auto m2 = m2b + nit;
      const auto i = m1 + l1;
      const auto j = m2 + l2;
      const auto previous = (i < 2 or j < 2) ? T{0} : a[idx(i - 2, j - 2)];
      auto value = -previous * std::sqrt(T(l1 + m1 - 1) * T(l1 - m1 + 2) *
                                         T(l2 - m2 + 2) * T(l2 + m2 - 1));
      value += T(iss + 2LL * (m1 - 1) * (m2 - 1)) * a[idx(i - 1, j - 1)];
      value /=
          std::sqrt(T(l1 + m1) * T(l1 - m1 + 1) * T(l2 - m2 + 1) * T(l2 + m2));
      a[idx(i, j)] = value;
    }
  }

  // Fill the half-plane m2 > 0 using the reflection symmetry of the
  // symbols under negation of all orders.
  {
    const auto sgn = ((l1 + l2 + is) % 2 == 0) ? T{1} : T{-1};
    for (int i = 0; i < nRows; ++i) {
      const auto m1 = i - l1;
      for (int m2 = 1; m2 <= l2; ++m2) {
        a[idx(i, l2 + m2)] = sgn * a[idx(l1 - m1, l2 - m2)];
      }
    }
  }

  // Restore the (-1)^m1 style phase by negating alternate rows.
  for (int i = 0; i < nRows; ++i) {
    if ((i + l2 + is) % 2 != 0) {
      for (int j = 0; j < nCols; ++j) {
        a[idx(i, j)] = -a[idx(i, j)];
      }
    }
  }
}

/**
 * @internal
 * @brief Transforms, in place, Woodhouse's matrix into the table of plain
 * Wigner 3j symbols
 *
 *   plain[(m1 + l1) * (2 l3 + 1) + (m3 + l3)]
 *       = ( l1     l2      l3 )
 *         ( m1  -(m1+m3)   m3 ),
 *
 * using plain(m1, m3) = (-1)^m1 wood(-m1, m3): rows are mirrored in m1 and
 * multiplied by an alternating phase.
 */
template <NumericConcepts::Real T>
void WoodhouseToPlain(const int l1, const int l3, std::span<T> a) {
  const auto nCols = static_cast<std::size_t>(2 * l3 + 1);
  for (int r = 0; r < l1; ++r) {
    const auto phase = ((r + l1) % 2 == 0) ? T{1} : T{-1};
    auto* rowLower = a.data() + static_cast<std::size_t>(r) * nCols;
    auto* rowUpper = a.data() + static_cast<std::size_t>(2 * l1 - r) * nCols;
    for (std::size_t j = 0; j < nCols; ++j) {
      const auto tmp = rowLower[j];
      rowLower[j] = phase * rowUpper[j];
      rowUpper[j] = phase * tmp;
    }
  }
}

}  // namespace Internal

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
  Internal::WoodhouseMatrix<T>(l1, l2, l3, span);
  Internal::WoodhouseToPlain<T>(l1, l3, span);
}

/**
 * @brief Fills user-provided storage with Woodhouse's matrix
 * (-1)^m (l1 l2 l3; -m, m-mp, mp) in row-major order, m being the row axis
 * and mp the column axis. This is exactly the array a(m+l1+1, mp+l3+1)
 * returned by the original Fortran routine wig2(l1, l2, l3, a, id1).
 * Allocation-free.
 * @tparam Range A contiguous, sized, readable and writable range of reals.
 * @param l1 The first degree.
 * @param l2 The second degree.
 * @param l3 The third degree.
 * @param table Storage of size (2 l1 + 1) x (2 l3 + 1).
 */
template <RealContiguousWritableRange Range>
void FillWoodhouseMatrix(int l1, int l2, int l3, Range&& table) {
  using T = std::ranges::range_value_t<Range>;
  auto span = std::span<T>(std::ranges::data(table), std::ranges::size(table));
  Internal::WoodhouseMatrix<T>(l1, l2, l3, span);
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
    Internal::WoodhouseMatrix<T>(_l1, _l2, _l3, _data);
    Internal::WoodhouseToPlain<T>(_l1, _l3, _data);
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
   * @brief Returns Woodhouse's matrix element
   * (-1)^m (l1 l2 l3; -m, m-mp, mp), the quantity tabulated by the
   * original routine wig2.f and common in normal-mode coupling theory.
   * @param m The row order, m in [-l1, l1].
   * @param mp The column order, mp in [-l3, l3].
   */
  auto Woodhouse(int m, int mp) const {
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
   * @brief Returns Woodhouse's matrix element
   * (-1)^m (l1 l2 l3; -m, m-mp, mp) for a given middle degree.
   * @param l2 The middle degree.
   * @param m The row order.
   * @param mp The column order.
   */
  auto Woodhouse(int l2, int m, int mp) const {
    if (not _l2Axis.Contains(l2)) return T{0};
    return Matrix(l2).Woodhouse(m, mp);
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
