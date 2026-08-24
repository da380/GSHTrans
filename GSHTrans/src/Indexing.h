#ifndef GSH_TRANS_INDEXING_GUARD_H
#define GSH_TRANS_INDEXING_GUARD_H

/**
 * @file Indexing.h
 * @brief Where a coefficient lives in a flat buffer.
 *
 * @details A block of coefficients at one upper index @f$n@f$ is stored degree
 * by degree, and within a degree order by order.
 *
 * Two things make the offset of a degree less obvious than it looks. The
 * degrees start at @f$|n|@f$ rather than zero, because @f$d^l_{nm}@f$ vanishes
 * identically below the upper index and there is no coefficient there to hold.
 * And the rows stop widening once the degree passes @p mMax: a block may be
 * truncated in order, which is what a real-valued field's reduced storage and
 * a Wigner table of limited order both are.
 *
 * So a row holds @f$2\min(l, m_{\max}) + 1@f$ coefficients when all orders
 * are stored and @f$\min(l, m_{\max}) + 1@f$ when only the non-negative ones
 * are, and the offset of a degree is the sum of the rows below it. GSHIndices
 * closes that sum rather than accumulating it, which is what makes a row
 * addressable without walking the block, and what lets the transform hand its
 * inner loop a row pointer per degree instead of an iterator across the whole
 * block.
 *
 * Bounds are checked by `assert` and so vanish under `NDEBUG`. That is the
 * library's rule for index arithmetic: the layers above validate degrees,
 * orders and sizes, and throw in every build mode when they are wrong.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <utility>

#include "Concepts.h"

namespace GSHTrans {

/**
 * @brief The orders stored at one degree, and where each sits in its row.
 * @tparam MRange Whether all orders are stored, or only the non-negative ones.
 */
template <OrderIndexRange MRange>
class GSHSubIndices {
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.

 public:
  /**
   * @brief The row for degree @p l, truncated at order @p mMax.
   * @param l The degree.
   * @param mMax The largest order stored anywhere in the block; this row's own
   * largest order is @f$\min(l, m_{\max})@f$.
   */
  constexpr GSHSubIndices(Int l, Int mMax) : _l{l}, _mMax{std::min(l, mMax)} {
    assert(_l >= 0);
    assert(_mMax >= 0);
  }

  /** @brief The degree this row belongs to. */
  constexpr auto Degree() const { return _l; }

  /** @brief The smallest order stored: @f$-m_{\max}@f$, or zero. */
  constexpr auto MinOrder() const {
    if constexpr (std::same_as<MRange, All>) {
      return -_mMax;
    } else {
      return Int{0};
    }
  }

  /** @brief The largest order stored. */
  constexpr auto MaxOrder() const { return _mMax; }

  /** @brief Every order stored, in storage order. */
  constexpr auto Orders() const {
    return std::ranges::views::iota(MinOrder(), MaxOrder() + 1);
  }

  /** @brief The negative orders, empty in reduced storage. */
  constexpr auto NegativeOrders() const {
    return std::ranges::views::iota(MinOrder(), 0);
  }

  /** @brief The orders from zero upwards. */
  constexpr auto NonNegativeOrders() const {
    return std::ranges::views::iota(0, MaxOrder() + 1);
  }

  /** @brief How many coefficients the row holds. */
  constexpr auto Size() const { return MaxOrder() - MinOrder() + 1; }

  /**
   * @brief Where order @p m sits within the row.
   * @param m The order, which must be one the row stores.
   */
  constexpr auto Index(Int m) const {
    if constexpr (std::same_as<MRange, All>) {
      assert(m >= -_l && m <= _l);
      return m + _mMax;
    } else {
      assert(m >= 0 && m <= _l);
      return m;
    }
  }

 private:
  Int _l;
  Int _mMax;
};

/**
 * @brief The layout of a whole coefficient block at one upper index.
 * @tparam MRange Whether all orders are stored, or only the non-negative ones.
 */
template <OrderIndexRange MRange>
class GSHIndices {
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.

 public:
  /** @brief An empty block: degree, order and upper index all zero. */
  constexpr GSHIndices() = default;

  /**
   * @brief A block over degrees @f$|n| \le l \le l_{\max}@f$.
   * @param lMax The largest degree stored.
   * @param mMax The largest order stored; rows stop widening beyond it.
   * @param n The upper index, which fixes the smallest degree at @f$|n|@f$.
   */
  constexpr GSHIndices(Int lMax, Int mMax, Int n)
      : _lMax{lMax}, _mMax{std::min(lMax, mMax)}, _n{n} {
    assert(_lMax >= 0);
    assert(_mMax >= 0);
    assert(std::abs(n) <= _lMax);
  }

  /** @brief The upper index this block belongs to. */
  constexpr auto UpperIndex() const { return _n; }

  /** @brief The largest order stored. */
  constexpr auto MaxOrder() const { return _mMax; }

  /** @brief The smallest degree stored, @f$|n|@f$. */
  constexpr auto MinDegree() const { return std::abs(_n); }
  /** @brief The largest degree stored. */
  constexpr auto MaxDegree() const { return _lMax; }
  /** @brief Every degree stored, in storage order. */
  constexpr auto Degrees() const {
    return std::ranges::views::iota(MinDegree(), MaxDegree() + 1);
  }

  /** @brief Every @f$(l, m)@f$ pair stored, in storage order. */
  constexpr auto Indices() const {
    return Degrees() | std::ranges::views::transform([this](auto l) {
             return std::ranges::views::cartesian_product(
                 std::ranges::views::single(l),
                 GSHSubIndices<MRange>(l, _mMax).Orders());
           }) |
           std::ranges::views::join;
  }

  /**
   * @brief Where the row for degree @p l begins.
   *
   * @details The sum of the rows below it, in closed form. A row holds
   * @f$2\min(l, m_{\max}) + 1@f$ coefficients, so summing from @f$|n|@f$
   * gives @f$l^2 - n^2@f$ while the rows are still widening, and a constant
   * @f$2m_{\max} + 1@f$ per degree once they have stopped. The branch is
   * which of those the degree falls in; the second case is a block whose
   * order limit is below its upper index, where no row ever widens.
   *
   * @param l The degree, which must be one the block stores.
   */
  constexpr auto OffsetForDegree(Int l) const
  requires std::same_as<MRange, All>
  {
    assert(l >= MinDegree() && l <= MaxDegree());
    auto nAbs = std::abs(_n);
    if (_mMax >= nAbs) {
      return l <= _mMax ? l * l - nAbs * nAbs
                        : (_mMax + 1) * (_mMax + 1) - nAbs * nAbs +
                              (l - 1 - _mMax) * (2 * _mMax + 1);
    } else {
      return (l - nAbs) * (2 * _mMax + 1);
    }
  }

  /**
   * @brief Where the row for degree @p l begins, in reduced storage.
   *
   * @details The same sum with rows of @f$\min(l, m_{\max}) + 1@f$ rather
   * than @f$2\min(l, m_{\max}) + 1@f$, so the widening part is a difference
   * of triangular numbers where the full-order form is a difference of
   * squares.
   *
   * @param l The degree, which must be one the block stores.
   */
  constexpr auto OffsetForDegree(Int l) const
  requires std::same_as<MRange, NonNegative>
  {
    assert(l >= MinDegree() && l <= MaxDegree());
    auto nAbs = std::abs(_n);
    if (_mMax >= nAbs) {
      return l <= _mMax
                 ? (l * (l + 1)) / 2 - (nAbs * (nAbs + 1)) / 2
                 : ((_mMax + 1) * (_mMax + 2)) / 2 - (nAbs * (nAbs + 1)) / 2 +
                       (l - 1 - _mMax) * (_mMax + 1);
    } else {
      return (l - nAbs) * (_mMax + 1);
    }
  }

  /** @brief How many coefficients the row for degree @p l holds. */
  constexpr auto SizeForDegree(Int l) const {
    assert(l >= MinDegree() && l <= MaxDegree());
    return GSHSubIndices<MRange>(l, _mMax).Size();
  }

  /**
   * @brief The offset of degree @p l, paired with its row's own indexing.
   * @return The offset, and a GSHSubIndices for the row.
   */
  constexpr auto Index(Int l) const {
    assert(l >= MinDegree() && l <= MaxDegree());
    return std::pair(OffsetForDegree(l), GSHSubIndices<MRange>(l, _mMax));
  }

  /** @brief Where the coefficient at degree @p l and order @p m sits. */
  constexpr auto Index(Int l, Int m) const {
    auto [offset, indices] = Index(l);
    return offset + indices.Index(m);
  }

  /** @brief How many coefficients the whole block holds. */
  constexpr auto Size() const {
    return OffsetForDegree(_lMax) + SizeForDegree(_lMax);
  }

 private:
  Int _lMax{};
  Int _mMax{};
  Int _n{};
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_INDEXING_GUARD_H