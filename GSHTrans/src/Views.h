#ifndef GSH_TRANS_VIEWS_GUARD_H
#define GSH_TRANS_VIEWS_GUARD_H

/**
 * @file Views.h
 * @brief Non-owning views over a coefficient buffer, indexed by degree and
 * order.
 *
 * @details A coefficient buffer is a flat array; these give it the
 * @f$(l, m)@f$ structure without copying or owning anything. Each pairs a
 * pointer with a GSHIndices, which supplies the arithmetic, and adds nothing
 * else — so a view is two words and is meant to be made, used within one
 * expression, and let go.
 *
 * There are four rather than two because constness is carried in the type
 * rather than in a template parameter: GSHView writes, ConstGSHView reads, and
 * the `Sub` forms are one degree's row of orders, which is what an inner loop
 * walks.
 *
 * Indexing is unchecked in release builds, deliberately. These are the
 * innermost accessors in the library — the Legendre stage indexes one per
 * value — and the layer above them has already validated degrees and orders.
 * The `assert`s in GSHIndices are what catch a mistake during development.
 */

#include <cstddef>
#include <iterator>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                              Writable views                               //
//--------------------------------------------------------------------------//

/**
 * @brief A writable view of the orders stored at one degree.
 * @tparam Scalar The buffer's value type.
 * @tparam MRange Whether all orders are stored, or only the non-negative ones.
 */
template <RealOrComplexFloatingPoint Scalar, OrderIndexRange MRange>
class GSHSubView : public GSHSubIndices<MRange> {
 public:
  using Int = std::ptrdiff_t;

  /**
   * @brief Views the row for degree @p l.
   * @param l The degree.
   * @param mMax The largest order stored, which may truncate the row.
   * @param start The first element of the row.
   */
  constexpr GSHSubView(Int l, Int mMax, Scalar* start)
      : GSHSubIndices<MRange>(l, mMax), _start{start} {}

  /** @brief The first stored order. */
  constexpr auto begin() { return _start; }
  /** @brief One past the last stored order. */
  constexpr auto end() { return std::next(begin(), this->Size()); }

  /** @brief The coefficient at order @p m. */
  constexpr auto& operator[](Int m) { return _start[this->Index(m)]; }

 private:
  Scalar* _start;
};

/**
 * @brief A writable view of a whole coefficient block at one upper index.
 * @tparam Scalar The buffer's value type.
 * @tparam MRange Whether all orders are stored, or only the non-negative ones.
 */
template <RealOrComplexFloatingPoint Scalar, OrderIndexRange MRange>
class GSHView : public GSHIndices<MRange> {
 public:
  using Int = std::ptrdiff_t;

  /**
   * @brief Views a block laid out as GSHIndices describes.
   * @param lMax The largest degree stored.
   * @param mMax The largest order stored.
   * @param n The upper index, which fixes the smallest degree at @f$|n|@f$.
   * @param start The first element of the block.
   */
  constexpr GSHView(Int lMax, Int mMax, Int n, Scalar* start)
      : GSHIndices<MRange>(lMax, mMax, n), _start{start} {}

  /** @brief The first coefficient. */
  constexpr auto begin() { return _start; }
  /** @brief One past the last coefficient. */
  constexpr auto end() { return std::next(begin(), this->Size()); }

  /**
   * @brief The row of orders at degree @p l.
   * @details This is the row supplier an inner loop takes once per degree,
   * rather than indexing @f$(l, m)@f$ per element.
   */
  constexpr auto operator[](Int l) {
    return GSHSubView<Scalar, MRange>(
        l, this->MaxOrder(), std::next(begin(), this->OffsetForDegree(l)));
  }

  /** @brief The coefficient at degree @p l and order @p m. */
  constexpr auto& operator[](Int l, Int m) { return _start[this->Index(l, m)]; }

 private:
  Scalar* _start;
};

//--------------------------------------------------------------------------//
//                             Read-only views                               //
//--------------------------------------------------------------------------//

/**
 * @brief A read-only view of the orders stored at one degree.
 * @copydetails GSHSubView
 */
template <RealOrComplexFloatingPoint Scalar, OrderIndexRange MRange>
class ConstGSHSubView : public GSHSubIndices<MRange> {
 public:
  using Int = std::ptrdiff_t;

  /**
   * @brief Views the row for degree @p l.
   * @param l The degree.
   * @param mMax The largest order stored, which may truncate the row.
   * @param start The first element of the row.
   */
  constexpr ConstGSHSubView(Int l, Int mMax, const Scalar* start)
      : GSHSubIndices<MRange>(l, mMax), _start{start} {}

  /** @brief The first stored order. */
  constexpr auto begin() const { return _start; }
  /** @brief One past the last stored order. */
  constexpr auto end() const { return std::next(begin(), this->Size()); }

  /** @brief The coefficient at order @p m. */
  constexpr auto operator[](Int m) const { return _start[this->Index(m)]; }

 private:
  const Scalar* _start;
};

/**
 * @brief A read-only view of a whole coefficient block at one upper index.
 * @copydetails GSHView
 */
template <RealOrComplexFloatingPoint Scalar, OrderIndexRange MRange>
class ConstGSHView : public GSHIndices<MRange> {
 public:
  using Int = std::ptrdiff_t;

  /**
   * @brief Views a block laid out as GSHIndices describes.
   * @param lMax The largest degree stored.
   * @param mMax The largest order stored.
   * @param n The upper index, which fixes the smallest degree at @f$|n|@f$.
   * @param start The first element of the block.
   */
  constexpr ConstGSHView(Int lMax, Int mMax, Int n, const Scalar* start)
      : GSHIndices<MRange>(lMax, mMax, n), _start{start} {}

  /** @brief The first coefficient. */
  constexpr auto begin() const { return _start; }
  /** @brief One past the last coefficient. */
  constexpr auto end() const { return std::next(begin(), this->Size()); }

  /**
   * @brief The row of orders at degree @p l.
   * @details The row supplier the Legendre stage takes once per degree. It is
   * also the substitution point at which a grid generating its Wigner values
   * hands back a row of scratch rather than a row of a stored table: the two
   * have the same type, because a view carries no storage.
   */
  constexpr auto operator[](Int l) const {
    return ConstGSHSubView<Scalar, MRange>(
        l, this->MaxOrder(), std::next(begin(), this->OffsetForDegree(l)));
  }

  /** @brief The coefficient at degree @p l and order @p m. */
  constexpr auto operator[](Int l, Int m) const {
    return _start[this->Index(l, m)];
  }

 private:
  const Scalar* _start;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_VIEWS_GUARD_H
