#ifndef GSH_TRANS_WIGNER_MATRICES_GUARD_H
#define GSH_TRANS_WIGNER_MATRICES_GUARD_H

#include <omp.h>

#include <cassert>
#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

#include "Concepts.h"
#include "Indexing.h"
#include "Views.h"
#include "Wigner.h"

namespace GSHTrans {

// The same d-function values as Wigner, laid out transform-major:
// [n][m][l][theta] rather than [n][theta][(l, m)].
//
// This is core-plan.md step M1, the first half of the transform-major
// restructure of section 11. The matrix kernel writes the Legendre stage as
// one matrix product per order,
//
//     f^n_{lm} = sum_i D^(n,m)_{li} b^(m)_i,     D^(n,m)_{li} = X^n_{lm}(theta_i)
//
// and that needs D^(n,m) contiguous in (l, theta) at fixed (n, m). The table
// Wigner builds is contiguous in (l, m) at fixed (n, theta), which is the
// transpose of what is wanted. This class holds the other one.
//
// The total size is unchanged, which is worth stating because it is not
// obvious: summing (lMax - max(|n|,|m|) + 1) over m gives exactly the same
// count as summing (2 min(l, mMax) + 1) over l. The same triangle, read the
// other way. At lMax = 256 that is 66,045 values per upper index either way.
//
// -- Why a separate class, and not a layout flag on Wigner.
//
// The two have unrelated interfaces. Wigner hands out a (l, m) block for one
// (n, theta), through ConstGSHView; this hands out an (l, theta) matrix for
// one (n, m), as a plain span a BLAS call can take. A single class serving
// both would carry two interfaces and would have to branch in the accessor
// the inner loop calls, which is the one place that cannot afford it. The
// grid already holds its table in an optional -- empty on a generating grid --
// so a second optional beside it is the shape that was already there.
//
// The Storage tag on Wigner (ColumnMajor / RowMajor) orders the (n, theta)
// axes only. This is a finer and different axis order and does not replace it.
//
// -- Why the values are generated rather than transposed.
//
// Building a Wigner table and transposing it would need both live at once:
// 1.3 GB at lMax = 256, nMax = 2, to end with 648 MB. So the recursion is run
// directly into per-thread scratch, one (n, theta) block at a time, and
// scattered into place. The recursion, the seeds, the evaluation order and
// the orthonormalisation are Wigner's own -- WignerDetails::ComputeBlock is
// called with the same arguments -- so the values are bit-identical to the
// stored path by construction rather than by tolerance. M1's test checks
// exactly that.
//
// -- The write pattern, named because it is the one cost here.
//
// A block computed for one (n, theta) scatters across every matrix, at stride
// NumberOfAngles() in the degree. That is one cache line per value written in
// the worst case. It is paid once, at construction, in parallel, and section
// 11's M1 records what it measures; blocking over colatitudes would fix it
// and is not done until something says it needs fixing.
template <RealFloatingPoint _Real, OrderIndexRange _MRange = All,
          IndexRange _NRange = All>
class WignerMatrices {
 public:
  using Int = std::ptrdiff_t;
  using Real = _Real;
  using MRange = _MRange;
  using NRange = _NRange;

  WignerMatrices() = default;

  template <std::ranges::range Range>
  requires RealFloatingPoint<std::ranges::range_value_t<Range>>
  WignerMatrices(Int lMax, Int mMax, Int nMax, Range &&theta)
      : _lMax{lMax},
        _mMax{mMax},
        _nMax{nMax},
        _nTheta(static_cast<Int>(std::ranges::size(theta))) {
    if (lMax < 0) throw std::invalid_argument("Maximum degree must be positive");
    if (mMax < 0 || mMax > lMax) {
      throw std::invalid_argument(
          "Maximum order must lie between zero and the maximum degree");
    }
    if (std::abs(nMax) > lMax) {
      throw std::invalid_argument(
          "Maximum upper index must not exceed the maximum degree");
    }

    // One offset per (n, m), and the matrices laid end to end in that order.
    _offset.reserve(
        static_cast<std::size_t>(NumberOfUpperIndices() * NumberOfOrders()));
    auto size = std::size_t{0};
    for (auto n : UpperIndices()) {
      for (auto m : Orders()) {
        _offset.push_back(size);
        size += static_cast<std::size_t>(NumberOfDegrees(n, m) * _nTheta);
      }
    }
    _data = std::vector<Real>(size);

    ComputeAll(theta);
  }

  // Degrees, orders, upper indices and angles. The upper-index accessors
  // match Wigner's exactly, since a grid hands both the same NRange.
  auto MaxDegree() const { return _lMax; }
  auto MaxOrder() const { return _mMax; }

  auto MinOrder() const {
    if constexpr (std::same_as<MRange, All>) {
      return -_mMax;
    } else {
      return Int{0};
    }
  }

  auto Orders() const {
    return std::ranges::views::iota(MinOrder(), MaxOrder() + 1);
  }

  auto NumberOfOrders() const { return MaxOrder() - MinOrder() + 1; }

  auto MinUpperIndex() const {
    if constexpr (std::same_as<NRange, All>) {
      return -_nMax;
    } else if constexpr (std::same_as<NRange, NonNegative>) {
      return Int{0};
    } else {
      return _nMax;
    }
  }

  auto MaxUpperIndex() const { return _nMax; }

  auto UpperIndices() const {
    return std::ranges::views::iota(MinUpperIndex(), MaxUpperIndex() + 1);
  }

  auto NumberOfUpperIndices() const {
    return MaxUpperIndex() - MinUpperIndex() + 1;
  }

  auto NumberOfAngles() const { return _nTheta; }

  // The lowest degree present at (n, m), and how many there are.
  //
  // A d-function vanishes identically unless l >= |n| and l >= |m|, so the
  // matrix at (n, m) starts at max(|n|, |m|) and its height falls linearly in
  // |m|. That is the source of the load imbalance step M4 has to divide work
  // for rather than count orders.
  auto MinDegree(Int n, Int m) const {
    assert(std::abs(n) <= _nMax);
    assert(m >= MinOrder() && m <= MaxOrder());
    return std::max(std::abs(n), std::abs(m));
  }

  auto NumberOfDegrees(Int n, Int m) const {
    return _lMax - MinDegree(n, m) + 1;
  }

  auto Degrees(Int n, Int m) const {
    return std::ranges::views::iota(MinDegree(n, m), _lMax + 1);
  }

  // The matrix D^(n,m): NumberOfDegrees(n, m) rows by NumberOfAngles()
  // columns, row-major, so the value at (l, iTheta) is at
  //
  //     (l - MinDegree(n, m)) * NumberOfAngles() + iTheta
  //
  // and the leading dimension a BLAS call wants is NumberOfAngles().
  //
  // The index is left to the caller rather than wrapped in an accessor,
  // deliberately: the layout is the whole content of this class, and a test
  // that computes the index itself is testing the layout rather than trusting
  // a member that could be wrong in the same way twice.
  auto operator[](Int n, Int m) const {
    return std::span<const Real>(
        _data.data() + _offset[OffsetIndex(n, m)],
        static_cast<std::size_t>(NumberOfDegrees(n, m) * _nTheta));
  }

 private:
  Int _lMax{};
  Int _mMax{};
  Int _nMax{};
  Int _nTheta{};

  std::vector<Real> _data;
  std::vector<std::size_t> _offset;

  auto OffsetIndex(Int n, Int m) const {
    assert(n >= MinUpperIndex() && n <= MaxUpperIndex());
    assert(m >= MinOrder() && m <= MaxOrder());
    return static_cast<std::size_t>((n - MinUpperIndex()) * NumberOfOrders() +
                                    (m - MinOrder()));
  }

  template <std::ranges::range Range>
  requires RealFloatingPoint<std::ranges::range_value_t<Range>>
  void ComputeAll(Range &&thetaRange) {
    const auto [sqrtInt, sqrtIntInv] =
        WignerDetails::PreComputeTables<Real>(_lMax, _mMax, _nMax);
    const auto sqrtIntView = std::span<const Real>(sqrtInt);
    const auto sqrtIntInvView = std::span<const Real>(sqrtIntInv);

    // The largest (n, theta) block is the one at the smallest |n|, since a
    // block starts at degree |n|. One scratch buffer of that size per thread.
    const auto scratchSize = static_cast<std::size_t>(
        GSHIndices<MRange>(_lMax, _mMax, SmallestUpperIndexModulus()).Size());

    // Flattened to an integer loop and decoded inside, and the region
    // suppressed when one is already open, for the reasons Wigner::ComputeAll
    // gives: OpenMP's canonical loop form, and exactly one level threads.
    const auto count = NumberOfUpperIndices() * _nTheta;
    const auto minUpperIndex = MinUpperIndex();
    const auto nTheta = _nTheta;
    const bool nested = omp_in_parallel();

#pragma omp parallel if (!nested)
    {
      auto scratch = std::vector<Real>(scratchSize);

#pragma omp for schedule(static)
      for (Int index = 0; index < count; index++) {
        const auto n = minUpperIndex + index / nTheta;
        const auto iTheta = index % nTheta;

        // Exactly the call Wigner makes, into scratch instead of into place.
        // Orthonormalisation happens inside it, so what is scattered below is
        // already the stored value sqrt((2l+1)/(4 pi)) d^l_{nm}.
        auto indices = GSHIndices<MRange>(_lMax, _mMax, n);
        WignerDetails::ComputeBlock(
            GSHView<Real, MRange>(_lMax, _mMax, n, scratch.data()), n,
            thetaRange[iTheta], sqrtIntView, sqrtIntInvView);

        // Scatter: (l, m) in the block goes to column iTheta of matrix (n, m).
        for (auto l : indices.Degrees()) {
          auto [blockOffset, sub] = indices.Index(l);
          for (auto m : sub.Orders()) {
            _data[_offset[OffsetIndex(n, m)] +
                  static_cast<std::size_t>((l - MinDegree(n, m)) * nTheta +
                                           iTheta)] =
                scratch[blockOffset + sub.Index(m)];
          }
        }
      }
    }
  }

  // The upper index of smallest modulus this table carries, which is the one
  // whose (n, theta) block is largest.
  auto SmallestUpperIndexModulus() const {
    if (MinUpperIndex() <= 0 && MaxUpperIndex() >= 0) return Int{0};
    return MinUpperIndex() > 0 ? MinUpperIndex() : MaxUpperIndex();
  }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_WIGNER_MATRICES_GUARD_H
