#ifndef GSH_TRANS_SPHERICAL_GRID_GUARD_H
#define GSH_TRANS_SPHERICAL_GRID_GUARD_H

#include <omp.h>

#include <FFTWpp/Core>
#include <FFTWpp/Ranges>
#include <NumericConcepts/Ranges.hpp>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Blas.h"
#include "Concepts.h"
#include "Indexing.h"
#include "Policies.h"
#include "Utility.h"
#include "Wigner.h"
#include "WignerMatrices.h"

namespace GSHTrans {

/**
 * @brief A grid on the sphere, and the transform over it.
 *
 * @details Holds everything that does not depend on *which* quadrature
 * produced the colatitudes: the Wigner values in whichever layout, both
 * Legendre kernels, the Fourier stages, the batch descriptors, chunking,
 * threading and the plan cache. A derived grid supplies the nodes and the
 * weights through the protected constructor and adds no data of its own, so
 * copying one into this base loses only what the derived class added.
 *
 * All a transform needs of a quadrature is a list of colatitudes, a list of
 * weights, a longitude count, and that the longitudes are uniform from zero.
 *
 * Not a CRTP base: one cannot call `Derived().CoLatitudes()` from a base's own
 * constructor, and building the Wigner values at construction is exactly that
 * call. The nodes are passed in and held here instead, which also leaves the
 * point-set helpers nothing to defer.
 *
 * @tparam _Real The precision.
 * @tparam _MRange Whether coefficient blocks hold all orders or only the
 * non-negative ones. A grid over the latter serves real scalars alone.
 * @tparam _NRange Which upper indices the grid covers.
 */
template <RealFloatingPoint _Real, OrderIndexRange _MRange, IndexRange _NRange>
class SphericalGrid {
 public:
  // Public type aliases.
  using Int = std::ptrdiff_t;          ///< Signed index type used throughout.
  using Real = _Real;                  ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  /// Whether all orders are stored, or only the non-negative ones.
  using MRange = _MRange;
  using NRange = _NRange;  ///< Which upper indices are covered.

  // A grid is a value-semantic handle over an immutable, shared
  // implementation: constructing one builds the Wigner values and copying one
  // copies a pointer.
  //
  // That is not a convenience. The table at lMax = 256, nMax = 2 with all
  // upper indices is about 648 MB, against 2.1 MB for a complex field on the
  // same grid, and the members used to be held by value with defaulted copy:
  // copying a grid by accident was not a performance wart but an
  // out-of-memory event.
  //
  // There is no default constructor: a grid without nodes is not a grid, and
  // reading MaxDegree() on a default-constructed one was undefined.
  SphericalGrid() = delete;

 protected:
  /// The contract a derived grid meets, and it is short.
  ///
  /// `coLatitudes` must be strictly increasing and lie strictly inside
  /// (0, pi); `weights` must be the quadrature weights for those nodes and the
  /// same length. Both are checked, because a grid that gets them wrong fails
  /// in the Wigner recursion rather than here.
  ///
  /// The interior condition is load-bearing beyond this class: Interpolate's
  /// polar padding rests on it, and it is
  /// where a scheme whose grid contains the pole would need thinking.
  ///
  /// `nPhi` defaults rather than being demanded, since it is about resolving
  /// orders |m| <= lMax and not about the quadrature; a scheme that needs a
  /// particular longitude count passes one.
  SphericalGrid(Int lMax, Int nMax, std::vector<Real> coLatitudes,
                std::vector<Real> coLatitudeWeights,
                FFTWpp::Flag flag = FFTWpp::Measure,
                Chunking chunking = Chunking::Automatic(),
                WignerValues values = WignerValues::Stored(),
                TransformKernel kernel = TransformKernel::Loop(), Int nPhi = 0)
      : _impl{std::make_shared<const Impl>(
            lMax, nMax, std::move(coLatitudes), std::move(coLatitudeWeights),
            nPhi > 0 ? nPhi : FastFFTSize(2 * lMax + 1), values, kernel)},
        _chunking{chunking},
        _flag{flag} {
    assert(flag != FFTWpp::WisdomOnly);
  }

 public:
  /** @brief Copy constructor. */
  SphericalGrid(const SphericalGrid&) = default;

  /** @brief Move constructor. */
  SphericalGrid(SphericalGrid&&) = default;

  /** @brief Copy assignment. */
  SphericalGrid& operator=(const SphericalGrid&) = default;

  /** @brief Move assignment. */
  SphericalGrid& operator=(SphericalGrid&&) = default;

  // Two grids are the same grid when they share an implementation. This is
  // the test the field layer's binary nodes make at construction. Structural
  // comparison is deliberately not offered: two separately built grids with
  // equal parameters hold different quadrature objects and different wisdom,
  // and treating them as interchangeable would make a node's operands
  // silently disagree about the buffers they index.
  /** @brief Identity, for deciding whether two share an implementation. */
  auto Identity() const { return _impl.get(); }

  /// The same grid with a different chunking policy, or a different planner
  /// flag: a pointer copy and a scalar, sharing one table.
  ///
  /// Offered for these two and for nothing else. WignerValues and
  /// TransformKernel each decide what the table *is*, so changing one means a
  /// different table -- which is a different grid, and the constructor is
  /// where you say so.
  ///
  /// **The result shares Identity() with its parent, and that is correct
  /// rather than a leak.** Identity is the field layer's test that two
  /// operands index the same buffers, and they do: same points, same degrees,
  /// same table, fields interchangeable. A chunk is how the inner loop
  /// schedules itself and is not observable in any result -- the batched tests
  /// demand exact equality against unbatched calls, which is the standing
  /// check that it is not. Sharing identity is also what makes this useful,
  /// since a tuned grid has to stay compatible with fields already built on
  /// the untuned one.
  auto With(Chunking chunking) const {
    auto grid = *this;
    grid._chunking = chunking;
    return grid;
  }

  /// The same grid under a different planner flag, sharing one table.
  auto With(FFTWpp::Flag flag) const {
    assert(flag != FFTWpp::WisdomOnly);
    auto grid = *this;
    grid._flag = flag;
    return grid;
  }

  /** @brief How many fields of a batch the inner loop takes at once. */
  auto ChunkingPolicy() const { return _chunking; }
  /** @brief How hard FFTW is asked to work at planning. */
  auto PlannerFlag() const { return _flag; }

  //------------------------------------------------//
  //                 The grid itself                 //
  //------------------------------------------------//
  /** @brief The largest degree stored. */
  auto MaxDegree() const { return _impl->lMax; }
  /** @brief The largest upper index covered. */
  auto MaxUpperIndex() const { return _impl->nMax; }

  /** @brief The colatitudes, strictly increasing in @f$(0,\pi)@f$. */
  auto CoLatitudes() const {
    return std::ranges::views::all(_impl->coLatitudes);
  }
  /** @brief The quadrature weights over colatitude. */
  auto CoLatitudeWeights() const {
    return std::ranges::views::all(_impl->coLatitudeWeights);
  }

  /** @brief The longitudes, uniform from zero. */
  auto Longitudes() const {
    const auto nPhi = NPhi();
    const auto dPhi = 2 * std::numbers::pi_v<Real> / static_cast<Real>(nPhi);
    return std::ranges::views::iota(Int{0}, nPhi) |
           std::ranges::views::transform([dPhi](auto i) { return i * dPhi; });
  }
  /** @brief The quadrature weights over longitude, which are uniform. */
  auto LongitudeWeights() const {
    const auto nPhi = NPhi();
    const auto dPhi = 2 * std::numbers::pi_v<Real> / static_cast<Real>(nPhi);
    return std::ranges::views::repeat(dPhi, nPhi);
  }

  //------------------------------------------------//
  //          The point set, and its sizes           //
  //------------------------------------------------//

  /** @brief The smallest upper index covered. */
  auto MinUpperIndex() const {
    if constexpr (std::same_as<NRange, All>) {
      return -MaxUpperIndex();
    } else if constexpr (std::same_as<NRange, NonNegative>) {
      return Int{0};
    } else {
      return MaxUpperIndex();
    }
  }

  /** @brief Every upper index covered. */
  auto UpperIndices() const {
    return std::ranges::views::iota(MinUpperIndex(), MaxUpperIndex() + 1);
  }

  /** @brief How many colatitudes the grid has. */
  auto NumberOfCoLatitudes() const { return CoLatitudes().size(); }
  /** @brief Indices of the colatitudes. */
  auto CoLatitudeIndices() const {
    return std::ranges::views::iota(Int{0},
                                    static_cast<Int>(NumberOfCoLatitudes()));
  }

  /** @brief How many longitudes the grid has. */
  auto NumberOfLongitudes() const { return Longitudes().size(); }
  /** @brief Indices of the longitudes. */
  auto LongitudeIndices() const {
    return std::ranges::views::iota(Int{0},
                                    static_cast<Int>(NumberOfLongitudes()));
  }

  /** @brief Every @f$(\theta,\phi)@f$ point, in storage order. */
  auto Points() const {
    return std::ranges::views::cartesian_product(CoLatitudes(), Longitudes());
  }

  /** @brief Every point index pair, in storage order. */
  auto PointIndices() const {
    return std::ranges::views::cartesian_product(CoLatitudeIndices(),
                                                 LongitudeIndices());
  }

  /** @brief The quadrature weight at every point. */
  auto Weights() const {
    return std::ranges::views::cartesian_product(CoLatitudeWeights(),
                                                 LongitudeWeights()) |
           std::ranges::views::transform(
               [](auto pair) { return std::get<0>(pair) * std::get<1>(pair); });
  }

  /**
   * @brief A callable sampled at every point of the grid, in storage order.
   *
   * @details Const because a value-semantic grid should be as usable through
   * a const handle as through a mutable one, and this reads nothing but the
   * point set.
   *
   * @param f The callable, invoked as `f(theta, phi)`.
   */
  template <typename Function>
  auto ProjectFunction(Function f) const {
    return Points() | std::ranges::views::transform([f](auto pair) {
             auto [theta, phi] = pair;
             return f(theta, phi);
           });
  }

  /** @brief How many samples one angular field holds. */
  auto FieldSize() const {
    return NumberOfCoLatitudes() * NumberOfLongitudes();
  }

  /// Number of coefficients of a complex-valued field of degree lMax at upper
  /// index n. This is the full (all orders) storage.
  auto CoefficientSize(Int lMax, Int n) const {
    return GSHIndices<All>(lMax, lMax, n).Size();
  }

  /// The same at the grid's own maximum degree.
  auto CoefficientSize(Int n) const { return CoefficientSize(MaxDegree(), n); }

  /// Number of stored coefficients of a real-valued field of degree lMax,
  /// which uses the reduced m >= 0 storage. There is no upper-index argument
  /// because real-valued fields exist only at upper index zero; the reduced
  /// storage is a statement about n = 0 alone.
  auto RealCoefficientSize(Int lMax) const {
    return GSHIndices<NonNegative>(lMax, lMax, 0).Size();
  }

  /// The same at the grid's own maximum degree.
  auto RealCoefficientSize() const { return RealCoefficientSize(MaxDegree()); }

  //-----------------------------------------------------//
  //          Forward transformation for ranges          //
  //-----------------------------------------------------//
  /**
   * @brief Analysis: a batch of fields to their coefficients.
   *
   * @details The two descriptors are separate because even when the caller's
   * arrangement is the same kind on both sides its `dist` differs — the fields
   * are FieldSize apart and the coefficient blocks CoefficientSize apart.
   * Equal counts is the only thing tying them together, and it is checked.
   *
   * The sizes are checked against Batch::Span rather than for equality, since
   * an interleaved batch is a window onto a larger range whose other elements
   * are no business of this call.
   *
   * The output is assigned rather than accumulated into, so transforming twice
   * into the same buffer gives the same answer twice.
   *
   * @param lMax The largest degree to resolve, at most the grid's own.
   * @param n The upper index, which every field of the batch shares.
   * @param in The samples, real or complex.
   * @param inBatch How the fields are arranged in @p in.
   * @param out Where the coefficients go.
   * @param outBatch How the blocks are arranged in @p out.
   * @param policy Whether the transform may thread.
   */
  template <NumericConcepts::RealOrComplexRange InRange,
            NumericConcepts::ComplexWritableRange OutRange>
  requires requires() {
    // A complex-valued field needs all orders in the coefficient storage; a
    // real-valued one uses the reduced m >= 0 storage and does not.
    requires std::same_as<_MRange, All> or NumericConcepts::RealRange<InRange>;
    // Field and coefficients both carry this grid's precision.
    requires std::same_as<NumericConcepts::RangePrecision<InRange>, Real>;
    requires std::same_as<std::ranges::range_value_t<OutRange>, Complex>;
  }
  void ForwardTransformation(Int lMax, Int n, InRange&& in, Batch inBatch,
                             OutRange& out, Batch outBatch,
                             Execution policy = Execution::Sequential()) const {
    // Get scalar type for field.
    using Scalar = std::ranges::range_value_t<InRange>;

    ValidateTransformRequest<Scalar>(lMax, n);

    const auto fieldSize = static_cast<Int>(this->FieldSize());
    const auto coefficientSize =
        static_cast<Int>(CoefficientSizeFor<Scalar>(lMax, n));

    // Two descriptors, not one. Even when the caller's arrangement is the
    // same kind on both sides, its dist differs: the fields are FieldSize
    // apart and the coefficient blocks are CoefficientSize apart. Equal
    // counts is the only thing tying them together, and it is checked.
    CheckBatch(inBatch, outBatch, fieldSize, coefficientSize);

    // A span check, not the equality the unbatched entry makes. An
    // interleaved batch is a window onto a larger range whose other elements
    // are no business of this call, so requiring an exact size would reject
    // precisely the layout the descriptor exists to accept.
    CheckSpan(std::ranges::size(in), inBatch.Span(fieldSize), "field");
    CheckSpan(std::ranges::size(out), outBatch.Span(coefficientSize),
              "coefficient");

    const auto count = inBatch.Count();
    auto outFirst = std::ranges::begin(out);
    auto inFirst = std::ranges::begin(in);

    // A one-point grid needs no FFT.
    if (_impl->lMax == 0) {
      for (auto k = Int{0}; k < count; k++) {
        outFirst[outBatch.Offset(0, k)] = inFirst[inBatch.Offset(0, k)] *
                                          static_cast<Real>(2) /
                                          std::numbers::inv_sqrtpi_v<Real>;
      }
      return;
    }

    const auto nPhi = static_cast<Int>(this->NumberOfLongitudes());
    const auto nTheta = static_cast<Int>(this->NumberOfCoLatitudes());
    const auto scaleFactor = static_cast<Real>(2) * std::numbers::pi_v<Real> /
                             static_cast<Real>(nPhi);

#ifdef GSHTRANS_HAVE_BLAS
    // The matrix kernel is a different arrangement of the same sum, not a
    // different sum: all the FFTs first, then one matrix product per order.
    // It is reached only when the grid was built for it, and the loop
    // below is untouched by its existence.
    //
    // The `if constexpr` is what keeps a long double grid compiling: BLAS has
    // no such precision, so the body must not be instantiated there. The
    // constructor has already refused the combination, so the discarded
    // branch is unreachable as well as uninstantiated.
    if constexpr (BlasDetails::BlasReal<Real>)
      if (_impl->kernel.IsMatrix()) {
        ForwardMatrixKernel<Scalar>(lMax, n, in, inBatch, outFirst, outBatch,
                                    count, nPhi, nTheta, scaleFactor, policy);
        return;
      }
#endif

    ForwardLoopKernel<Scalar>(lMax, n, inFirst, inBatch, outFirst, outBatch,
                              count, nPhi, nTheta, scaleFactor, coefficientSize,
                              policy);
  }

  // The single field, which is the batched primitive at count = 1.
  //
  // The size checks here are equalities rather than spans, because this
  // entry's contract is that the range *is* the field: a caller who hands
  // over a range of the wrong length has made a mistake, whereas a batched
  // caller may legitimately hand over a window onto a larger one.
  /**
   * @brief Analysis: one field to its coefficients.
   *
   * @details The single field, which is the batched primitive at
   * `count = 1`. The size checks here are equalities rather than spans,
   * because this entry's contract is that the range *is* the field.
   *
   * @param lMax The largest degree to resolve, at most the grid's own.
   * @param n The upper index of the field.
   * @param in The samples, real or complex.
   * @param out Where the coefficients go.
   * @param policy Whether the transform may thread.
   */
  template <NumericConcepts::RealOrComplexRange InRange,
            NumericConcepts::ComplexWritableRange OutRange>
  requires requires() {
    requires std::same_as<_MRange, All> or NumericConcepts::RealRange<InRange>;
    requires std::same_as<NumericConcepts::RangePrecision<InRange>, Real>;
    requires std::same_as<std::ranges::range_value_t<OutRange>, Complex>;
  }
  void ForwardTransformation(Int lMax, Int n, InRange&& in, OutRange& out,
                             Execution policy = Execution::Sequential()) const {
    using Scalar = std::ranges::range_value_t<InRange>;
    ValidateTransformRequest<Scalar>(lMax, n);
    const auto fieldSize = static_cast<Int>(this->FieldSize());
    const auto coefficientSize =
        static_cast<Int>(CoefficientSizeFor<Scalar>(lMax, n));
    CheckSize(std::ranges::size(in), fieldSize, "field");
    CheckSize(std::ranges::size(out), coefficientSize, "coefficient");
    ForwardTransformation(lMax, n, in, Batch::One(fieldSize), out,
                          Batch::One(coefficientSize), policy);
  }

  //------------------------------------------------//
  //       Inverse  transformation for ranges       //
  //------------------------------------------------//
  /**
   * @brief Synthesis: a batch of coefficient blocks to their fields.
   *
   * @details The mirror of ForwardTransformation, and the same remarks apply
   * to the two descriptors and to the size checks.
   *
   * @param lMax The largest degree present, at most the grid's own.
   * @param n The upper index, which every block of the batch shares.
   * @param in The coefficients.
   * @param inBatch How the blocks are arranged in @p in.
   * @param out Where the samples go.
   * @param outBatch How the fields are arranged in @p out.
   * @param policy Whether the transform may thread.
   */
  template <NumericConcepts::ComplexRange InRange,
            NumericConcepts::RealOrComplexWritableRange OutRange>
  requires requires() {
    // As above, read the other way round: the field is the output here.
    requires std::same_as<_MRange, All> or
                 NumericConcepts::RealWritableRange<OutRange>;
    requires std::same_as<NumericConcepts::RangePrecision<OutRange>, Real>;
    requires std::same_as<std::ranges::range_value_t<InRange>, Complex>;
  }
  void InverseTransformation(Int lMax, Int n, InRange&& in, Batch inBatch,
                             OutRange& out, Batch outBatch,
                             Execution policy = Execution::Sequential()) const {
    // Get scalar type for field.
    using Scalar = std::ranges::range_value_t<OutRange>;

    ValidateTransformRequest<Scalar>(lMax, n);

    const auto fieldSize = static_cast<Int>(this->FieldSize());
    const auto coefficientSize =
        static_cast<Int>(CoefficientSizeFor<Scalar>(lMax, n));

    CheckBatch(inBatch, outBatch, coefficientSize, fieldSize);
    CheckSpan(std::ranges::size(in), inBatch.Span(coefficientSize),
              "coefficient");
    CheckSpan(std::ranges::size(out), outBatch.Span(fieldSize), "field");

    const auto count = inBatch.Count();
    auto inFirst = std::ranges::begin(in);
    auto outFirst = std::ranges::begin(out);

    // A one-point grid needs no FFT.
    if (_impl->lMax == 0) {
      for (auto k = Int{0}; k < count; k++) {
        const auto value = inFirst[inBatch.Offset(0, k)] *
                           std::numbers::inv_sqrtpi_v<Real> /
                           static_cast<Real>(2);
        if constexpr (RealFloatingPoint<Scalar>) {
          outFirst[outBatch.Offset(0, k)] = std::real(value);
        } else {
          outFirst[outBatch.Offset(0, k)] = value;
        }
      }
      return;
    }

    const auto nPhi = static_cast<Int>(this->NumberOfLongitudes());
    const auto nTheta = static_cast<Int>(this->NumberOfCoLatitudes());

#ifdef GSHTRANS_HAVE_BLAS
    if constexpr (BlasDetails::BlasReal<Real>)
      if (_impl->kernel.IsMatrix()) {
        InverseMatrixKernel<Scalar>(lMax, n, inFirst, inBatch, out, outBatch,
                                    count, nPhi, nTheta, policy);
        return;
      }
#endif

    InverseLoopKernel<Scalar>(lMax, n, inFirst, inBatch, outFirst, outBatch,
                              count, nPhi, nTheta, coefficientSize, policy);
  }

  // The single field, which is the batched primitive at count = 1.
  /**
   * @brief Synthesis: one coefficient block to its field.
   *
   * @details The single field, which is the batched primitive at
   * `count = 1`; the size checks are equalities rather than spans.
   *
   * @param lMax The largest degree present, at most the grid's own.
   * @param n The upper index of the field.
   * @param in The coefficients.
   * @param out Where the samples go.
   * @param policy Whether the transform may thread.
   */
  template <NumericConcepts::ComplexRange InRange,
            NumericConcepts::RealOrComplexWritableRange OutRange>
  requires requires() {
    requires std::same_as<_MRange, All> or
                 NumericConcepts::RealWritableRange<OutRange>;
    requires std::same_as<NumericConcepts::RangePrecision<OutRange>, Real>;
    requires std::same_as<std::ranges::range_value_t<InRange>, Complex>;
  }
  void InverseTransformation(Int lMax, Int n, InRange&& in, OutRange& out,
                             Execution policy = Execution::Sequential()) const {
    using Scalar = std::ranges::range_value_t<OutRange>;
    ValidateTransformRequest<Scalar>(lMax, n);
    const auto fieldSize = static_cast<Int>(this->FieldSize());
    const auto coefficientSize =
        static_cast<Int>(CoefficientSizeFor<Scalar>(lMax, n));
    CheckSize(std::ranges::size(in), coefficientSize, "coefficient");
    CheckSize(std::ranges::size(out), fieldSize, "field");
    InverseTransformation(lMax, n, in, Batch::One(coefficientSize), out,
                          Batch::One(fieldSize), policy);
  }

  // -- The matrix kernel's Fourier stage.
  //
  /// How many longitudinal Fourier coefficients a field of this scalar type
  /// has: nPhi for a complex field, nPhi / 2 + 1 for a real one, whose
  /// negative orders are the conjugates of its positive ones.
  template <RealOrComplexFloatingPoint Scalar>
  auto FourierSize() const {
    const auto nPhi = static_cast<Int>(this->NumberOfLongitudes());
    return static_cast<Int>(FFTWpp::DataSize<Scalar, Complex>(nPhi).second);
  }

  /// The size of buffer ForwardFourierStage fills for `count` fields.
  template <RealOrComplexFloatingPoint Scalar>
  auto ForwardFourierStageSize(Int count) const {
    return FourierSize<Scalar>() *
           static_cast<Int>(this->NumberOfCoLatitudes()) * count;
  }

  /// The longitudinal DFT of `count` fields of a batch, starting at field
  /// `first`, written m-major:
  ///
  ///     out[m * (nTheta * count) + iTheta * count + k]
  ///
  /// so that for one order m the (theta, k) block is contiguous with k
  /// fastest. That block is what the per-order matrix product
  /// multiplies: a complex (nTheta x count) matrix, which -- because
  /// std::complex stores its parts adjacently -- is also a real
  /// (nTheta x 2 count) one, and so can go to dgemm with N = 2 count.
  ///
  /// The loop kernel interleaves this stage with the Legendre stage, one
  /// colatitude at a time, because that keeps its working set to a single row.
  /// A per-order product cannot: it needs every colatitude of one order at
  /// once, so all the FFTs have to run first. That is the trade being made,
  /// and it is why this is a separate entry point rather than a change to the
  /// loop kernel, which is kept exactly as it is.
  ///
  /// -- No scaling is applied here. The quadrature weight and the 2 pi / nPhi
  /// both belong to the Legendre stage, which is where the loop kernel applies
  /// them too, so the two kernels can be compared value by value at this seam.
  ///
  /// -- `thetaBlock` is how many colatitudes are transformed per FFTW call, and
  /// zero asks the library to choose. It is a hint in the sense Chunking is:
  /// the library rounds it to something it will not regret.
  ///
  /// Section 12 of the reference note says the transpose is free, on the
  /// grounds that one plan_many with output stride howMany and output distance
  /// 1 lands the data in the order above at no cost, exactly as tier-1
  /// batching landed it in [m][k] order. **Measured, that is true only
  /// away from one specific hazard, and false at it.** The output write has
  /// stride `howMany` complex doubles, so when `howMany * 16` is a power of two
  /// the writes for successive orders collide in the same cache sets, and the
  /// FFT costs four to eight times what it costs at howMany plus or minus two.
  /// Measured directly at nPhi = 520, per transform: 0.80 at howMany = 62,
  /// 2.96 at 64, 0.83 at 66; 0.90 at 254, 7.32 at 256, 0.88 at 258; 1.76 at
  /// 2046, 8.53 at 2048, 1.65 at 2050.
  ///
  /// This is the hypothesis raised for RadialMajor and **rejected** there,
  /// because its tiling already handled it. Here
  /// nothing tiles: FFTW writes straight through, so the same hazard bites.
  /// The guard is to shrink the block until the product is not a power of two,
  /// which costs nothing and also reduces the workspace.
  ///
  /// Blocking also bounds the FFT buffers, which `out` does not: at
  /// lMax = 256 with eight fields a full-height call wants two of about 16 MB
  /// each, per thread, while a block of three wants 0.4 MB. Section 11.3
  /// priced the intermediate at 17 MB and missed that there is a second buffer
  /// of the same size on the input side.
  ///
  /// The cost of blocking is one contiguous copy per order per block -- a
  /// memcpy, not a gather, since both sides are contiguous in (theta, k) -- and
  /// it is measured at under a fifth of the stage.
  template <std::ranges::input_range InRange>
  void ForwardFourierStage(InRange&& in, Batch inBatch, Int first, Int count,
                           std::span<Complex> out, Int thetaBlock = 0,
                           Execution policy = Execution::Sequential()) const {
    using Scalar = std::ranges::range_value_t<InRange>;
    static_assert(RealOrComplexFloatingPoint<Scalar>);

    const auto nPhi = static_cast<Int>(this->NumberOfLongitudes());
    const auto nTheta = static_cast<Int>(this->NumberOfCoLatitudes());
    const auto nFourier = FourierSize<Scalar>();

    if (count < 1) {
      throw std::invalid_argument("Field count must be positive");
    }
    if (first < 0 || first + count > inBatch.Count()) {
      throw std::invalid_argument("Requested fields lie outside the batch");
    }
    CheckSpan(std::ranges::size(in),
              inBatch.Span(static_cast<Int>(this->FieldSize())), "field");
    if (static_cast<Int>(out.size()) !=
        ForwardFourierStageSize<Scalar>(count)) {
      throw std::invalid_argument("Fourier stage output has the wrong size");
    }

    const auto block = ChooseThetaBlock(thetaBlock, count, nTheta);
    auto inFirst = std::ranges::begin(in);

    // Threaded over blocks. Each block reads its own colatitudes of
    // the caller's fields and writes its own run within every order's output,
    // so the blocks are disjoint on both sides; the workspace and its plan are
    // thread_local already, so a thread finds or makes its own.
    //
    // Left sequential at first, which made it a 29 per cent Amdahl term as
    // soon
    // as the Legendre stage threaded -- and the matrix kernel stopped scaling
    // past two threads because of it.
    const auto blocks = (nTheta + block - 1) / block;
    const bool parallel = RunInParallel(policy);

#pragma omp parallel for schedule(static) \
    num_threads(ThreadCount(policy)) if (parallel)
    for (Int b = 0; b < blocks; b++) {
      const auto theta0 = b * block;
      const auto rows = std::min(block, nTheta - theta0);
      auto& work = GetWorkspace<Scalar, true>(nPhi, rows * count, _flag);

      // Pack (theta, k) rows in that order, so the FFT's own transform index
      // runs theta-major with k fastest -- which is precisely the order the
      // output block wants.
      for (auto iRow = Int{0}; iRow < rows; iRow++) {
        for (auto k = Int{0}; k < count; k++) {
          PackRow(std::next(inFirst,
                            inBatch.Offset((theta0 + iRow) * nPhi, first + k)),
                  nPhi, inBatch.Stride(),
                  std::next(work.in.begin(), (iRow * count + k) * nPhi));
        }
      }
      work.plan.Execute();

      // One contiguous run per order. The workspace holds
      // [m][thetaLocal][k] and the target holds [m][theta][k], so the run for
      // order m starts at m * rows * count in one and at
      // m * nTheta * count + theta0 * count in the other, and has the same
      // length in both.
      const auto run = rows * count;
      for (auto m = Int{0}; m < nFourier; m++) {
        const auto* source = work.out.data() + m * run;
        std::copy_n(source, run,
                    out.begin() + m * nTheta * count + theta0 * count);
      }
    }
  }

  /// The inverse of ForwardFourierStage: an m-major intermediate in, `count`
  /// fields of a batch out, starting at field `first`.
  ///
  ///     in[m * (nTheta * count) + iTheta * count + k]
  ///
  /// The caller owns the intermediate and must have set **every** order it
  /// holds, including the ones no coefficient reaches: the transform is over
  /// nPhi orders whatever the degree, and the band between lMax and
  /// nPhi - lMax carries no coefficient but is still read. The matrix kernel
  /// zeroes exactly that band; the loop kernel zeroes its whole row buffer for
  /// the same reason.
  ///
  /// Blocked and guarded exactly as the forward stage is, and for the same
  /// reasons -- see there, including the power-of-two hazard, which is a
  /// property of the strided *read* here rather than the strided write.
  template <typename OutRange>
  void InverseFourierStage(std::span<const Complex> in, OutRange& out,
                           Batch outBatch, Int first, Int count,
                           Int thetaBlock = 0,
                           Execution policy = Execution::Sequential()) const {
    using Scalar = std::ranges::range_value_t<OutRange>;
    static_assert(RealOrComplexFloatingPoint<Scalar>);

    const auto nPhi = static_cast<Int>(this->NumberOfLongitudes());
    const auto nTheta = static_cast<Int>(this->NumberOfCoLatitudes());
    const auto nFourier = FourierSize<Scalar>();

    if (count < 1) {
      throw std::invalid_argument("Field count must be positive");
    }
    if (first < 0 || first + count > outBatch.Count()) {
      throw std::invalid_argument("Requested fields lie outside the batch");
    }
    CheckSpan(std::ranges::size(out),
              outBatch.Span(static_cast<Int>(this->FieldSize())), "field");
    if (static_cast<Int>(in.size()) != ForwardFourierStageSize<Scalar>(count)) {
      throw std::invalid_argument("Fourier stage input has the wrong size");
    }

    const auto block = ChooseThetaBlock(thetaBlock, count, nTheta);
    auto outFirst = std::ranges::begin(out);

    // Threaded over blocks, as the forward stage is and for the same reason.
    const auto blocks = (nTheta + block - 1) / block;
    const bool parallel = RunInParallel(policy);

#pragma omp parallel for schedule(static) \
    num_threads(ThreadCount(policy)) if (parallel)
    for (Int b = 0; b < blocks; b++) {
      const auto theta0 = b * block;
      const auto rows = std::min(block, nTheta - theta0);
      auto& work = GetWorkspace<Scalar, false>(nPhi, rows * count, _flag);

      const auto run = rows * count;
      for (auto m = Int{0}; m < nFourier; m++) {
        std::copy_n(in.begin() + m * nTheta * count + theta0 * count, run,
                    work.in.begin() + m * run);
      }
      work.plan.Execute();

      for (auto iRow = Int{0}; iRow < rows; iRow++) {
        for (auto k = Int{0}; k < count; k++) {
          UnpackRow(
              std::next(work.out.begin(), (iRow * count + k) * nPhi), nPhi,
              std::next(outFirst,
                        outBatch.Offset((theta0 + iRow) * nPhi, first + k)),
              outBatch.Stride());
        }
      }
    }
  }

 private:
  // Exactly one level threads. A transform asked to run in parallel from
  // inside an existing parallel region runs sequentially instead, so that a
  // caller parallelising over slices, components or realisations cannot nest
  // with this, and neither can Wigner::ComputeAll.
  static bool RunInParallel(Execution policy) {
    return policy.IsParallel() && !omp_in_parallel();
  }

  static int ThreadCount(Execution policy) {
    return policy.Threads() > 0 ? policy.Threads() : omp_get_max_threads();
  }

  // How many colatitudes ForwardFourierStage transforms per FFTW call.
  //
  // The default is small, which is the opposite of what the reference note
  // assumes and is what measurement says: at lMax = 256 with
  // eight fields the stage costs 3.1 ms at a block of three against 6.0 ms
  // with every colatitude in one call, and wants 0.4 MB of workspace rather
  // than 33 MB. A strided write over a few hundred bytes stays inside a
  // couple of cache lines; one over tens of kilobytes does not.
  //
  // Then the guard: shrink while the output stride in bytes is a power of two,
  // which is the cache-set collision measured at the call site above. Only
  // ever downwards, so it terminates and never asks for more memory than was
  // requested. A caller who forces a block gets it rounded the same way, since
  // the rounding costs nothing and the alternative is a silent factor of two.
  //
  // The final partial block takes whatever colatitudes are left and is not
  // adjusted -- it cannot be, since the remainder is what it is. It is at most
  // one block out of many.
  static Int ChooseThetaBlock(Int requested, Int count, Int nTheta) {
    constexpr auto defaultBlock = Int{3};
    auto block = requested > 0 ? std::min(requested, nTheta)
                               : std::min(defaultBlock, nTheta);
    const auto aliases = [count](Int rows) {
      const auto bytes = rows * count * static_cast<Int>(sizeof(Complex));
      return bytes >= 512 && (bytes & (bytes - 1)) == 0;
    };
    while (block > 1 && aliases(block)) block--;
    return block;
  }

  // A per-thread accumulator for the forward transform's partial sums, kept
  // between calls for the same reason the work buffers are: this is the size
  // of the coefficient array, and allocating it per call would put back the
  // an allocation back into every call. It only ever grows.
  static std::vector<Complex>& Accumulator(std::size_t size) {
    thread_local auto buffer = std::vector<Complex>{};
    if (buffer.size() < size) buffer.resize(size);
    return buffer;
  }

  // The aligned buffers a transform works in, and the FFTW plan bound to them.
  //
  // Held per thread rather than shared. The buffers are scratch and must be
  // per-thread whatever else happens; binding the plan to them keeps the two
  // together, lets execution use the plan's own buffers rather than
  // FFTW's new-array form, and leaves Impl immutable, which is what makes a
  // shared grid safe to use concurrently without a lock.
  // The cost is planning once per thread per shape instead of once per shape,
  // which after the first is a wisdom lookup.
  template <RealOrComplexFloatingPoint Scalar, bool IsForward>
  struct Workspace {
    using In = std::conditional_t<IsForward, Scalar, Complex>;
    using Out = std::conditional_t<IsForward, Complex, Scalar>;

    // `count` colatitude rows at a time, one per field of a batch chunk. The
    // single-field case is count = 1 and plans the same transform it always
    // did: with howMany = 1 the stride and dist below are unreachable.
    Workspace(Int nPhi, Int count, FFTWpp::Flag flag)
        : count{count},
          in(count * FFTWpp::DataSize<In, Out>(nPhi).first),
          out(count * FFTWpp::DataSize<In, Out>(nPhi).second),
          plan(MakePlan(in, out, nPhi, count, flag)) {}

    // The layouts are the whole of "batching the FFT is free": FFTWpp
    // passes rank, howMany, embed, stride and dist straight through to
    // plan_many, so the batch costs a descriptor rather than a repack.
    //
    // The asymmetry between the two sides is deliberate. On the *field* side
    // the rows are packed one after another, because that is what a copy from
    // the caller's storage produces however the caller strides. On the
    // *coefficient* side the batch index runs fastest, so the data lands in
    // [m][k] order and the Legendre stage's inner loop over the batch is
    // unit-stride. That layout costs nothing to ask for here and is the one
    // the batching measurement was made on.
    static auto MakePlan(FFTWpp::vector<In>& in, FFTWpp::vector<Out>& out,
                         Int nPhi, Int count, FFTWpp::Flag flag) {
      const auto sizes = FFTWpp::DataSize<In, Out>(nPhi);
      const auto inSize = static_cast<int>(sizes.first);
      const auto outSize = static_cast<int>(sizes.second);
      const auto howMany = static_cast<int>(count);

      // FFTWpp reads N as the extent of each transform's own array, which for
      // the half-spectrum side of a real transform is nPhi / 2 + 1 rather
      // than nPhi; DataSize already returns that, so both sides use it.
      const auto inN = std::array<int, 1>{inSize};
      const auto outN = std::array<int, 1>{outSize};

      auto inLayout =
          IsForward ? FFTWpp::Ranges::Layout(1, inN, howMany, inN, 1, inSize)
                    : FFTWpp::Ranges::Layout(1, inN, howMany, inN, howMany, 1);
      auto outLayout =
          IsForward
              ? FFTWpp::Ranges::Layout(1, outN, howMany, outN, howMany, 1)
              : FFTWpp::Ranges::Layout(1, outN, howMany, outN, 1, outSize);

      auto inView = FFTWpp::Ranges::View(in, inLayout);
      auto outView = FFTWpp::Ranges::View(out, outLayout);
      // No lock here. FFTW's planner is not re-entrant, but FFTWpp now takes
      // its own PlannerMutex inside every planner entry point, so serialising
      // again on ours would only add a second, coarser lock over the same
      // critical section -- and one that unrelated FFTWpp users could not see.
      // Execution stays unlocked in both: FFTW allows a plan to be executed
      // concurrently, and each thread executes on its own workspace anyway.
      if constexpr (std::same_as<In, Out>) {
        return FFTWpp::Ranges::Plan(
            inView, outView, flag,
            IsForward ? FFTWpp::Forward : FFTWpp::Backward);
      } else {
        return FFTWpp::Ranges::Plan(inView, outView, flag);
      }
    }

    // How far apart the same order of successive fields sits in the
    // coefficient-side buffer, and how far apart successive orders sit. The
    // Legendre stage reads both rather than assuming either.
    Int Count() const { return count; }

    Int count;
    FFTWpp::vector<In> in;
    FFTWpp::vector<Out> out;
    decltype(MakePlan(std::declval<FFTWpp::vector<In>&>(),
                      std::declval<FFTWpp::vector<Out>&>(), Int{}, Int{},
                      std::declval<FFTWpp::Flag>())) plan;
  };

  // Plans and buffers used to be created on every call -- two allocations and
  // a plan per transform, with FFTW's non-re-entrant planner run each time
  //. They are now made once per thread per shape and kept. Held by
  // pointer so that the plan's reference to its buffers survives any
  // rehashing of the cache.
  template <RealOrComplexFloatingPoint Scalar, bool IsForward>
  static Workspace<Scalar, IsForward>& GetWorkspace(Int nPhi, Int count,
                                                    FFTWpp::Flag flag) {
    using Entry = Workspace<Scalar, IsForward>;
    thread_local auto cache =
        std::map<std::tuple<Int, Int, unsigned>, std::unique_ptr<Entry>>{};
    const auto key = std::tuple{nPhi, count, static_cast<unsigned>(flag)};
    auto found = cache.find(key);
    if (found == cache.end()) {
      found =
          cache.emplace(key, std::make_unique<Entry>(nPhi, count, flag)).first;
    }
    return *found->second;
  }

  // Move one colatitude row between the caller's field and the plan's own
  // buffer.
  //
  // These are the only points at which caller storage is touched during a
  // transform. The FFT plans are created on the grid's own fftw_malloc'd
  // buffers and are executed on those buffers alone; FFTW's new-array execute
  // is valid only for storage with the same alignment characteristics as the
  // planning buffers, and neither FFTW nor FFTWpp checks.
  // The forward transform used to take that path on any writable input and
  // the inverse took it unconditionally on the caller's output, which made
  // alignment an obligation propagating outward into every stride the field
  // and tensor layers might choose. Copying instead costs one pass over a row
  // against an FFT of the same row, and it is what makes the field plan's
  // promise -- that a slice target needs no alignment beyond Scalar's --
  // true rather than aspirational.
  //
  // The caller's stride enters at these two points and at no other. A batch
  // may be interleaved with components this call knows nothing about, and
  // everything downstream of here works on our own contiguous buffers.
  template <typename Iterator, typename Target>
  static void PackRow(Iterator first, Int nPhi, Int stride, Target target) {
    if (stride == 1) {
      std::copy_n(first, nPhi, target);
      return;
    }
    for (auto i = Int{0}; i < nPhi; i++) {
      *target++ = *first;
      std::advance(first, stride);
    }
  }

  template <typename Source, typename Iterator>
  static void UnpackRow(Source source, Int nPhi, Iterator first, Int stride) {
    if (stride == 1) {
      std::copy_n(source, nPhi, first);
      return;
    }
    for (auto i = Int{0}; i < nPhi; i++) {
      *first = *source++;
      std::advance(first, stride);
    }
  }

  // How many fields the inner loop takes at once: the grid's chunking policy,
  // asked about this particular call.
  //
  // The two directions ask different questions, which is why there are two
  // entry points rather than one. What the policy needs to know is how many
  // copies of the coefficient block will be live at once, and that is a
  // property of the decomposition, not of the call:
  //
  //   forward -- every thread accumulates into a private buffer of
  //              chunk * coefficientSize, so the copies are the threads;
  //   inverse -- one block is gathered, shared and read-only, so there is one
  //              copy however many threads read it.
  //
  // Serving both with the thread count is what the policy used to do, and it
  // starved the inverse: at lMax = 256 and k = 8 on eight threads it returned
  // a chunk of one where the whole batch fits, and taking the whole batch
  // measured 2.2x faster. The forward's rule is
  // unchanged, and the same measurement is the evidence for that too -- the
  // whole batch there is 2x *slower*, because eight private accumulators of
  // 8.4 MB ask for 68 MB of a 16 MB cache.
  //
  // The thread count is one whenever the policy is sequential or a parallel
  // region is already open, since that is what the call will really run on.
  Int ForwardChunkSize(Int coefficientSize, Execution policy) const {
    const auto threads = RunInParallel(policy) ? ThreadCount(policy) : 1;
    return _chunking.Count(coefficientSize * static_cast<Int>(sizeof(Complex)),
                           threads);
  }

  Int InverseChunkSize(Int coefficientSize) const {
    return _chunking.Count(coefficientSize * static_cast<Int>(sizeof(Complex)),
                           1);
  }

  // The loop kernel of the forward transform: a colatitude at a time, over a
  // Wigner table contiguous in (l, m) at fixed (n, theta).
  //
  // A named member rather than the body of the public entry point, so that
  // the two kernels read as the peers they are: the entry point validates,
  // chooses and delegates, and both kernels are the same size of thing.
  template <RealOrComplexFloatingPoint Scalar, typename InIterator,
            typename OutIterator>
  void ForwardLoopKernel(Int lMax, Int n, InIterator inFirst, Batch inBatch,
                         OutIterator outFirst, Batch outBatch, Int count,
                         Int nPhi, Int nTheta, Real scaleFactor,
                         Int coefficientSize, Execution policy) const {
    // One colatitude's contribution from a chunk of `c` fields, accumulated
    // into a [coefficient][field] scratch buffer.
    //
    // The batch index runs fastest in both the FFT output and the scratch, so
    // the innermost loop is an axpy of length c over contiguous memory. That
    // is the whole of what batching buys at tier 1: the Wigner row for this
    // (n, iTheta, l) is read once and used c times, rather than re-streamed
    // per field. The Wigner traffic per field falls by c; nothing else about
    // the arithmetic changes.
    auto AccumulateRow = [&](Int iTheta, Int first, Int c, Complex* scratch,
                             auto& work) {
      // Copy this colatitude's row out of each field. Caller stride enters
      // here and at the unpack, and nowhere else.
      for (auto k = Int{0}; k < c; k++) {
        PackRow(std::next(inFirst, inBatch.Offset(iTheta * nPhi, first + k)),
                nPhi, inBatch.Stride(), std::next(work.in.begin(), k * nPhi));
      }
      work.plan.Execute();

      // Get the Wigner values and quadrature weight.
      auto d = WignerBlock(n, iTheta, lMax);
      const auto w =
          _impl->coLatitudeWeights[static_cast<std::size_t>(iTheta)] *
          scaleFactor;
      const auto orders = static_cast<Int>(work.out.size()) / c;

      // Loop over the spherical harmonic coefficients, taking the Wigner
      // values one degree at a time.
      //
      // The row pointer comes from d[l] rather than from a single iterator
      // walked across the whole block. That is the supplier seam: the only
      // thing this loop needs is a contiguous run of values in (l, m) order,
      // one run per degree, and asking for it per degree is what lets a
      // generating grid substitute a supplier that builds the row into
      // per-thread scratch for one that points into the stored table.
      // OffsetForDegree is closed-form, so the seam costs a few integer
      // operations per degree against a loop of length (2l+1) * c.
      //
      // The two constraints it carries, both already satisfied here: degrees
      // are visited in ascending contiguous order from |n|, and each
      // (n, iTheta) is visited once per pass. A generated row cannot be
      // revisited without re-running the recursion.
      auto* target = scratch;
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });

      for (auto l : degrees) {
        auto dl = d[l];
        auto wigIter = dl.begin();

        if constexpr (ComplexFloatingPoint<Scalar>) {
          // Negative orders live at the top of the FFT output, c apart.
          const auto* source = work.out.data() + (orders - dl.MaxOrder()) * c;
          for ([[maybe_unused]] auto m : dl.NegativeOrders()) {
            const auto a = *wigIter++ * w;
            for (auto k = Int{0}; k < c; k++) target[k] += a * source[k];
            target += c;
            source += c;
          }
          source = work.out.data();
          for ([[maybe_unused]] auto m : dl.NonNegativeOrders()) {
            const auto a = *wigIter++ * w;
            for (auto k = Int{0}; k < c; k++) target[k] += a * source[k];
            target += c;
            source += c;
          }
        } else {
          const auto* source = work.out.data();
          if constexpr (std::same_as<_MRange, All>) {
            std::advance(wigIter, dl.MaxOrder());
          }
          for ([[maybe_unused]] auto m : dl.NonNegativeOrders()) {
            const auto a = *wigIter++ * w;
            for (auto k = Int{0}; k < c; k++) target[k] += a * source[k];
            target += c;
            source += c;
          }
        }
      }
    };

    // Write a completed [coefficient][field] block to wherever the caller
    // keeps it. Assignment rather than accumulation, which is what makes the
    // routine own its output: transforming twice into the same buffer used to
    // double the answer, because the colatitude loop accumulates and nothing
    // initialised the destination. Zeroing the range first
    // would be both redundant and wrong here, since a range holding an
    // interleaved batch also holds components this call must not touch.
    auto Scatter = [&](const Complex* scratch, Int first, Int c, Int fromJ,
                       Int toJ) {
      for (auto j = fromJ; j < toJ; j++) {
        for (auto k = Int{0}; k < c; k++) {
          outFirst[outBatch.Offset(j, first + k)] = scratch[j * c + k];
        }
      }
    };

    const auto chunk = ForwardChunkSize(coefficientSize, policy);
    for (auto first = Int{0}; first < count; first += chunk) {
      const auto c = std::min(chunk, count - first);
      const auto scratchSize = static_cast<std::size_t>(coefficientSize * c);

      if (!RunInParallel(policy)) {
        auto& work = GetWorkspace<Scalar, true>(nPhi, c, _flag);
        auto& scratch = CoefficientScratch(scratchSize);
        std::fill_n(scratch.begin(), scratchSize, Complex{});
        for (auto iTheta = Int{0}; iTheta < nTheta; iTheta++) {
          AccumulateRow(iTheta, first, c, scratch.data(), work);
        }
        Scatter(scratch.data(), first, c, 0, coefficientSize);
        continue;
      }

      // Every colatitude contributes to every coefficient, so the colatitudes
      // cannot simply be divided between threads writing into `out`. Each
      // thread accumulates into a private buffer and the partial sums are
      // added at the end.
      //
      // That addition is *partitioned*, not serialised. Each thread owns one
      // block of degrees and sums every thread's partials for that block
      // alone, then scatters it, so the reduction runs in parallel and no
      // thread waits on another. It used to be a critical section in which
      // each thread added a whole coefficient array in turn, which costs one
      // serialised pass per thread: invisible against the colatitude loop at
      // eight threads, and 128 MB of serialised adds per transform at 128.
      // The decomposition itself is unchanged -- thread-private
      // accumulators over colatitudes are the wrong shape well before 128
      // threads, but choosing what replaces them needs a machine this was
      // not measured on.
      const auto threadCount = ThreadCount(policy);

      // The reduction reads every thread's accumulator, so the thread-local
      // buffers have to be published to the team. Written before the implicit
      // barrier at the end of the colatitude loop and read after it.
      auto partials = std::vector<Complex*>(threadCount, nullptr);

#pragma omp parallel num_threads(threadCount)
      {
        const auto thread = static_cast<Int>(omp_get_thread_num());
        const auto threads = static_cast<Int>(omp_get_num_threads());
        auto& work = GetWorkspace<Scalar, true>(nPhi, c, _flag);
        auto& partial = Accumulator(scratchSize);
        std::fill_n(partial.begin(), scratchSize, Complex{});
        partials[thread] = partial.data();

#pragma omp for schedule(static)
        for (Int iTheta = 0; iTheta < nTheta; iTheta++) {
          AccumulateRow(iTheta, first, c, partial.data(), work);
        }

        const auto fromJ = coefficientSize * thread / threads;
        const auto toJ = coefficientSize * (thread + 1) / threads;
        auto* base = partials[0];
        for (auto t = Int{1}; t < threads; t++) {
          const auto* p = partials[t];
          for (auto i = fromJ * c; i < toJ * c; i++) base[i] += p[i];
        }
        Scatter(base, first, c, fromJ, toJ);
      }
    }
  }

  // The loop kernel of the inverse transform, the mirror of the forward one.
  //
  // Unlike the forward direction the colatitudes write disjoint output and
  // share only read-only input, so threading over them needs no reduction.
  template <RealOrComplexFloatingPoint Scalar, typename InIterator,
            typename OutIterator>
  void InverseLoopKernel(Int lMax, Int n, InIterator inFirst, Batch inBatch,
                         OutIterator outFirst, Batch outBatch, Int count,
                         Int nPhi, Int nTheta, Int coefficientSize,
                         Execution policy) const {
    // One colatitude, synthesised into its own row of each field. Unlike the
    // forward transform, the colatitudes here write disjoint output and share
    // only read-only input, so threading over them needs no reduction.
    auto SynthesiseRow = [&](Int iTheta, Int first, Int c,
                             const Complex* scratch, auto& work) {
      std::ranges::for_each(work.in, [](auto& x) { return x = 0; });

      // Get the Wigner values.
      auto d = WignerBlock(n, iTheta, lMax);
      const auto orders = static_cast<Int>(work.in.size()) / c;

      // Loop over the coefficients, one degree at a time. As in the forward
      // direction, the row pointer comes from d[l]: this is the same supplier
      // seam, and a generating grid substitutes at the same point. The
      // inner loop is again an axpy of length c, over a batch index that
      // runs fastest on both sides.
      const auto* source = scratch;
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });
      for (auto l : degrees) {
        auto dl = d[l];
        auto wigIter = dl.begin();
        if constexpr (ComplexFloatingPoint<Scalar>) {
          auto* target = work.in.data() + (orders - dl.MaxOrder()) * c;
          for ([[maybe_unused]] auto m : dl.NegativeOrders()) {
            const auto a = *wigIter++;
            for (auto k = Int{0}; k < c; k++) target[k] += a * source[k];
            target += c;
            source += c;
          }
          target = work.in.data();
          for ([[maybe_unused]] auto m : dl.NonNegativeOrders()) {
            const auto a = *wigIter++;
            for (auto k = Int{0}; k < c; k++) target[k] += a * source[k];
            target += c;
            source += c;
          }
        } else {
          auto* target = work.in.data();
          if constexpr (std::same_as<_MRange, All>) {
            std::advance(wigIter, dl.MaxOrder());
          }
          for ([[maybe_unused]] auto m : dl.NonNegativeOrders()) {
            const auto a = *wigIter++;
            for (auto k = Int{0}; k < c; k++) target[k] += a * source[k];
            target += c;
            source += c;
          }
        }
      }

      // Perform FFT to recover the field at this colatitude, through the
      // plan's own buffers, then hand each row to the caller.
      work.plan.Execute();
      for (auto k = Int{0}; k < c; k++) {
        UnpackRow(
            std::next(work.out.begin(), k * nPhi), nPhi,
            std::next(outFirst, outBatch.Offset(iTheta * nPhi, first + k)),
            outBatch.Stride());
      }
    };

    const auto chunk = InverseChunkSize(coefficientSize);
    for (auto first = Int{0}; first < count; first += chunk) {
      const auto c = std::min(chunk, count - first);
      const auto scratchSize = static_cast<std::size_t>(coefficientSize * c);

      // Gather this chunk's coefficients into [coefficient][field] order
      // once, rather than reaching through the caller's stride on every
      // colatitude. The Legendre stage works on our own buffers, whose layout
      // we choose, and the choice is the one the batched loop above wants
      // above. It costs one pass over the coefficients against a colatitude
      // loop that reads the whole Wigner block, which measures
      // invisible for the same reason in the other direction.
      //
      // Gathered by the calling thread before the parallel region opens, and
      // read-only inside it.
      auto& scratch = CoefficientScratch(scratchSize);
      for (auto k = Int{0}; k < c; k++) {
        for (auto j = Int{0}; j < coefficientSize; j++) {
          scratch[j * c + k] = inFirst[inBatch.Offset(j, first + k)];
        }
      }
      const auto* gathered = scratch.data();

      if (!RunInParallel(policy)) {
        auto& work = GetWorkspace<Scalar, false>(nPhi, c, _flag);
        for (auto iTheta = Int{0}; iTheta < nTheta; iTheta++) {
          SynthesiseRow(iTheta, first, c, gathered, work);
        }
        continue;
      }

#pragma omp parallel num_threads(ThreadCount(policy))
      {
        auto& work = GetWorkspace<Scalar, false>(nPhi, c, _flag);
#pragma omp for schedule(static)
        for (Int iTheta = 0; iTheta < nTheta; iTheta++) {
          SynthesiseRow(iTheta, first, c, gathered, work);
        }
      }
    }
  }

#ifdef GSHTRANS_HAVE_BLAS
  // Scratch for the m-major Fourier intermediate and for one order's product.
  // Kept per thread and grown, never allocated per call: an allocation
  // inside a loop it was meant to
  // serve costing 3 to 5 times the operation itself.
  static std::vector<Complex>& MatrixScratch(std::size_t size) {
    thread_local auto buffer = std::vector<Complex>{};
    if (buffer.size() < size) buffer.resize(size);
    return buffer;
  }

  static std::vector<Complex>& OrderScratch(std::size_t size) {
    thread_local auto buffer = std::vector<Complex>{};
    if (buffer.size() < size) buffer.resize(size);
    return buffer;
  }

  // Run `body(m, scratch)` for every order from minOrder to lMax, threaded or
  // not.
  //
  // **schedule(dynamic), and the reason is that the obvious static split is
  // wrong twice over.** Section 12 of the reference note warns once: the work
  // per order is not constant, since n_L(m) falls linearly in |m|, so counting
  // orders is about twice as unbalanced as it looks and the schedule has to
  // divide work instead. M3b found the second and sharper reason. A static
  // split weighted by n_L would assume time proportional to n_L -- but the
  // inverse's inner dimension *is* n_L, and a GEMM's efficiency falls with its
  // inner dimension, so time is superlinear in n_L and a linear model
  // mis-splits in the same direction it was correcting. Dynamic needs no model
  // and so cannot hold a wrong one.
  //
  // Determinism is unaffected: the orders write disjoint output, so the answer
  // does not depend on which thread took which or in what sequence.
  //
  // The scratch is per thread and is taken *inside* the region for that
  // reason. Anything shared -- the table, the intermediate -- is captured by
  // reference from outside it, and must be: OrderScratch and MatrixScratch are
  // both thread_local, so a buffer filled before the region is the master
  // thread's and reaching for it again inside would find an empty one.
  template <typename Body>
  void OverOrders(Int minOrder, Int lMax, Int c, Execution policy,
                  Body&& body) const {
    // Room for a paired product and a paired right-hand side at once: step
    // The reflection puts the orders +m and -m through one GEMM, so the
    // widest case is
    // two columns of each, and the right-hand side is nTheta rows deep.
    const auto scratchSize = static_cast<std::size_t>(
        2 * c * (lMax + 1 + static_cast<Int>(this->NumberOfCoLatitudes())));
    const auto orders = lMax - minOrder + 1;

    // **A team is opened even for the sequential case, and that is the point
    // rather than an accident of the code.** M3a measured that a threaded
    // BLAS loses on these products -- 23.7 ms against 17.0 at
    // lMax = 256, k = 8, with GSHTrans sequential and so with no nesting to
    // blame -- because the 513 products are skinny and thread launch
    // dominates. The library cannot set a BLAS's thread count portably: there
    // is no standard call, and OPENBLAS_NUM_THREADS is a property of one
    // implementation and is inert in its OpenMP build.
    //
    // What it can do is make sure every GEMM is issued from *inside* an
    // OpenMP region. A BLAS on the same runtime is then nested, and the
    // default of one active level makes it serial without anyone being asked
    // to set anything. A BLAS on its own pthread pool is unaffected and still
    // needs the caller's environment, which is what the CMake comment says.
    //
    // Found by falling into it: the benchmark's sequential rows reported
    // table traffic at twice single-core bandwidth, because "GSHTrans
    // sequential" had been letting the BLAS take the whole machine.
    const auto threads = RunInParallel(policy) ? ThreadCount(policy) : 1;

#pragma omp parallel num_threads(threads)
    {
      auto& scratch = OrderScratch(scratchSize);
#pragma omp for schedule(dynamic)
      for (Int i = 0; i < orders; i++) body(minOrder + i, scratch.data());
    }
  }

  // The forward transform as one matrix product per order.
  //
  // Reordered so that the sum over colatitudes is innermost, the Legendre
  // stage is, at fixed upper index and order,
  //
  //     f^n_{lm} = sum_i D^(n,m)_{li} b^(m)_i,    b^(m)_i = w_i F_m(theta_i)
  //
  // and over a chunk of c fields the right-hand side is a matrix, so this is
  // a GEMM of shape (nL x nTheta) times (nTheta x c). D is real while the
  // data are complex, and that is a gift rather than an obstacle: a complex
  // (nTheta x c) block with the batch index fastest *is* a real
  // (nTheta x 2c) one, because std::complex stores its parts adjacently. So
  // the whole stage is dgemm with N = 2c, and the doubling helps most exactly
  // where c is small.
  //
  // Degrees below max(|n|, |m|) do not exist, so the matrix at each order
  // starts there and its height falls linearly in |m|. A transform at a
  // degree below the grid's takes a contiguous *prefix* of the rows, at the
  // same leading dimension, so nothing is copied for a truncated call.
  template <RealOrComplexFloatingPoint Scalar, std::ranges::input_range InRange,
            typename OutIterator>
  void ForwardMatrixKernel(Int lMax, Int n, InRange&& in, Batch inBatch,
                           OutIterator outFirst, Batch outBatch, Int count,
                           Int nPhi, Int nTheta, Real scaleFactor,
                           Execution policy) const {
    using MRangeForScalar =
        std::conditional_t<RealFloatingPoint<Scalar>, NonNegative, All>;

    const auto& matrices = *_impl->wignerMatrices;
    const auto indices = GSHIndices<MRangeForScalar>(lMax, lMax, n);
    const auto coefficientSize = static_cast<Int>(indices.Size());
    const auto nFourier = FourierSize<Scalar>();

    // A real field has no negative orders to pair with: its coefficients at
    // -m are the conjugates of those at +m and are not stored at all.
    constexpr bool paired = !RealFloatingPoint<Scalar>;

    const auto chunk = InverseChunkSize(coefficientSize);

    for (auto first = Int{0}; first < count; first += chunk) {
      const auto c = std::min(chunk, count - first);

      // All the FFTs, landing [m][theta][k].
      auto& stage =
          MatrixScratch(static_cast<std::size_t>(nFourier * nTheta * c));
      auto stageSpan = std::span<Complex>(
          stage.data(), static_cast<std::size_t>(nFourier * nTheta * c));
      ForwardFourierStage(in, inBatch, first, c, stageSpan, 0, policy);

      // Everything one order needs, and nothing another order touches. That
      // disjointness is the whole of the threading argument: the
      // products at different orders read a shared, read-only table and a
      // shared, read-only intermediate, and write coefficients no other order
      // writes. **No accumulator and no reduction, in either direction** --
      // which is the forward loop kernel's weak point deleted rather than
      // tuned.
      // One order and, where there is one, its negation -- both against the
      // *same* matrix.
      //
      //     f^n_{l,-m} = (-1)^{l+n} sum_j w_j D^(n,m)_{lj} F_{-m}(theta_j-bar)
      //
      // so the product at -m is the product at +m applied to the -m Fourier
      // data in reversed colatitude order, with a sign on the output rows.
      // Two consequences, and the second is the one that pays.
      //
      // The table halves, because only m >= 0 is stored. And the two
      // right-hand sides go into **one** GEMM rather than two, which doubles
      // N from 2c to 4c -- N being the skinniest dimension in the problem and
      // the one measured as limiting, at 55 to 95 Gflop/s against peak.
      //
      // The arithmetic does *not* halve: the same products are still done,
      // of the same shapes.
      auto DoOrder = [&](Int m, Complex* scratch) {
        const auto lMin = std::max(std::abs(n), std::abs(m));
        const auto rows = lMax - lMin + 1;
        const auto pairs = paired && m > 0 ? Int{2} : Int{1};

        const auto* a = matrices[n, m].data();
        auto* plus = stage.data() + m * nTheta * c;

        // The +m half is scaled in place and multiplied where it lies, as it
        // was before the reflection existed. Only the -m half is copied, and
        // it has to be: the reflection reverses the colatitude, and no BLAS
        // takes a negative stride.
        for (auto iTheta = Int{0}; iTheta < nTheta; iTheta++) {
          const auto w =
              _impl->coLatitudeWeights[static_cast<std::size_t>(iTheta)] *
              scaleFactor;
          for (auto k = Int{0}; k < c; k++) plus[iTheta * c + k] *= w;
        }

        BlasDetails::RowMajorGemm(
            static_cast<int>(rows), static_cast<int>(2 * c),
            static_cast<int>(nTheta), Real{1}, a, static_cast<int>(nTheta),
            reinterpret_cast<const Real*>(plus), static_cast<int>(2 * c),
            Real{0}, reinterpret_cast<Real*>(scratch), static_cast<int>(2 * c));

        // Scatter: the coefficient block is triangular, so the stride between
        // consecutive degrees at one order is not constant and the product
        // cannot be written into it directly.
        for (auto l = lMin; l <= lMax; l++) {
          const auto j = indices.Index(l, m);
          for (auto k = Int{0}; k < c; k++) {
            outFirst[outBatch.Offset(j, first + k)] =
                scratch[(l - lMin) * c + k];
          }
        }

        if (pairs == 1) return;

        auto* rhs = scratch + rows * c;
        const auto* minus = stage.data() + (nPhi - m) * nTheta * c;
        for (auto iTheta = Int{0}; iTheta < nTheta; iTheta++) {
          const auto w =
              _impl->coLatitudeWeights[static_cast<std::size_t>(iTheta)] *
              scaleFactor;
          const auto mirror = nTheta - 1 - iTheta;
          for (auto k = Int{0}; k < c; k++) {
            rhs[iTheta * c + k] = minus[mirror * c + k] * w;
          }
        }

        BlasDetails::RowMajorGemm(
            static_cast<int>(rows), static_cast<int>(2 * c),
            static_cast<int>(nTheta), Real{1}, a, static_cast<int>(nTheta),
            reinterpret_cast<const Real*>(rhs), static_cast<int>(2 * c),
            Real{0}, reinterpret_cast<Real*>(scratch), static_cast<int>(2 * c));

        for (auto l = lMin; l <= lMax; l++) {
          const auto sign = std::remove_cvref_t<decltype(matrices)>::Sign(l, n);
          const auto jMinus = indices.Index(l, -m);
          for (auto k = Int{0}; k < c; k++) {
            outFirst[outBatch.Offset(jMinus, first + k)] =
                sign * scratch[(l - lMin) * c + k];
          }
        }
      };

      OverOrders(Int{0}, lMax, c, policy, DoOrder);
    }
  }

  // The inverse transform, the same matrices read the other way.
  //
  //     F_m(theta_i) = sum_l D^(n,m)_{li} f^n_{lm}
  //
  // which is D transposed applied to the coefficients, so **one stored matrix
  // serves both directions** and the transpose flag is the whole of the
  // difference. Nothing is copied and no second table exists.
  //
  // No quadrature weight here, which is why it could not have been folded
  // into the matrix in the forward direction.
  template <RealOrComplexFloatingPoint Scalar, typename InIterator,
            std::ranges::range OutRange>
  void InverseMatrixKernel(Int lMax, Int n, InIterator inFirst, Batch inBatch,
                           OutRange& out, Batch outBatch, Int count, Int nPhi,
                           Int nTheta, Execution policy) const {
    using MRangeForScalar =
        std::conditional_t<RealFloatingPoint<Scalar>, NonNegative, All>;

    const auto& matrices = *_impl->wignerMatrices;
    const auto indices = GSHIndices<MRangeForScalar>(lMax, lMax, n);
    const auto coefficientSize = static_cast<Int>(indices.Size());
    const auto nFourier = FourierSize<Scalar>();

    // A real field has no negative orders to pair with: its coefficients at
    // -m are the conjugates of those at +m and are not stored at all.
    constexpr bool paired = !RealFloatingPoint<Scalar>;

    const auto chunk = InverseChunkSize(coefficientSize);

    for (auto first = Int{0}; first < count; first += chunk) {
      const auto c = std::min(chunk, count - first);

      auto& stage =
          MatrixScratch(static_cast<std::size_t>(nFourier * nTheta * c));

      // Only the orders a coefficient reaches are written by the products
      // below; the band between them carries nothing and has to be zeroed,
      // because the FFT reads every order whatever the degree, and because
      // this buffer is kept between calls and so holds whatever the last one
      // left there. A transform below the grid's degree leaves the widest
      // band and is where removing this shows up first.
      //
      // The band alone, not the whole intermediate: at lMax = 256 the whole
      // is 16 MiB and this is the part of it no product touches.
      const auto zeroFrom = lMax + 1;
      const auto zeroTo = RealFloatingPoint<Scalar> ? nFourier : nPhi - lMax;
      for (auto m = zeroFrom; m < zeroTo; m++) {
        std::fill_n(stage.data() + m * nTheta * c, nTheta * c, Complex{0, 0});
      }

      // Disjoint in the same way the forward direction is, and for the same
      // reason: each order reads coefficients no other order reads and writes
      // the one block of the intermediate that its own FFT order occupies.
      // The mirror of the forward direction, order by order. With
      // g_l = (-1)^{l+n} f^n_{l,-m}, the reflection gives
      //
      //     F_{-m}(theta_i) = sum_l D^(n,m)_{l,i-bar} g_l
      //
      // so the -m result comes out of the same product, in reversed
      // colatitude order, from coefficients that carry the sign on the way in.
      auto DoOrder = [&](Int m, Complex* gathered) {
        const auto lMin = std::max(std::abs(n), std::abs(m));
        const auto rows = lMax - lMin + 1;
        const auto pairs = paired && m > 0 ? Int{2} : Int{1};

        const auto* a = matrices[n, m].data();

        // The +m product writes straight into its block of the intermediate,
        // as it did before the reflection existed. Only the -m product needs
        // somewhere to land first, because its rows come out mirrored.
        for (auto l = lMin; l <= lMax; l++) {
          const auto j = indices.Index(l, m);
          for (auto k = Int{0}; k < c; k++) {
            gathered[(l - lMin) * c + k] =
                inFirst[inBatch.Offset(j, first + k)];
          }
        }

        BlasDetails::RowMajorGemmTransposed(
            static_cast<int>(nTheta), static_cast<int>(2 * c),
            static_cast<int>(rows), Real{1}, a, static_cast<int>(nTheta),
            reinterpret_cast<const Real*>(gathered), static_cast<int>(2 * c),
            Real{0}, reinterpret_cast<Real*>(stage.data() + m * nTheta * c),
            static_cast<int>(2 * c));

        if (pairs == 1) return;

        // The sign rides in on the coefficients, so the product is the -m
        // field with its colatitudes reversed.
        for (auto l = lMin; l <= lMax; l++) {
          const auto sign = std::remove_cvref_t<decltype(matrices)>::Sign(l, n);
          const auto jMinus = indices.Index(l, -m);
          for (auto k = Int{0}; k < c; k++) {
            gathered[(l - lMin) * c + k] =
                sign * inFirst[inBatch.Offset(jMinus, first + k)];
          }
        }

        auto* result = gathered + rows * c;
        BlasDetails::RowMajorGemmTransposed(
            static_cast<int>(nTheta), static_cast<int>(2 * c),
            static_cast<int>(rows), Real{1}, a, static_cast<int>(nTheta),
            reinterpret_cast<const Real*>(gathered), static_cast<int>(2 * c),
            Real{0}, reinterpret_cast<Real*>(result), static_cast<int>(2 * c));

        auto* minus = stage.data() + (nPhi - m) * nTheta * c;
        for (auto iTheta = Int{0}; iTheta < nTheta; iTheta++) {
          // Row iTheta of the product is F_{-m} at the *mirrored* angle.
          const auto mirror = nTheta - 1 - iTheta;
          const auto* row = result + iTheta * c;
          for (auto k = Int{0}; k < c; k++) minus[mirror * c + k] = row[k];
        }
      };

      OverOrders(Int{0}, lMax, c, policy, DoOrder);

      auto stageSpan = std::span<const Complex>(
          stage.data(), static_cast<std::size_t>(nFourier * nTheta * c));
      InverseFourierStage(stageSpan, out, outBatch, first, c, 0, policy);
    }
  }
#endif

  // The Wigner values for one (n, iTheta), over the degrees |n| .. lMax.
  //
  // This is the supplier seam, and it returns the same type on
  // both paths: ConstGSHView carries no storage, being (lMax, mMax, n,
  // const Real*) over the index arithmetic it inherits from GSHIndices, so a
  // view of generated scratch and a view into the table are indistinguishable
  // to the loops that consume them. Neither transform changes by a line.
  //
  // On a stored grid this hands back a pointer into the table. On a generating
  // one it runs the recursion into this thread's scratch and returns a view of
  // that -- the same recursion, the same order, the same values, so the two
  // paths agree bit for bit rather than to a tolerance.
  //
  // The generated block is built to the *call's* degree, not the grid's, so a
  // truncated transform generates only the degrees it uses. The stored path
  // cannot do that: its rows are laid out for the grid's maximum degree
  // whatever a call asks for. The orders are asked for to the same bound,
  // which makes every degree's row full-width -- min(l, lMax) is l -- and so
  // lays the block out exactly as the table lays out its own prefix.
  //
  // The view points into thread_local scratch and is valid until this thread
  // asks for another block. Both consumers use it within one colatitude and
  // then let it go.
  auto WignerBlock(Int n, Int iTheta, Int lMax) const {
    if (_impl->wigner) return (*_impl->wigner)[n, iTheta];

    const auto size =
        static_cast<std::size_t>(GSHIndices<_MRange>(lMax, lMax, n).Size());
    auto& scratch = WignerScratch(size);
    WignerDetails::ComputeBlock(
        GSHView<Real, _MRange>(lMax, lMax, n, scratch.data()), n,
        _impl->coLatitudes[static_cast<std::size_t>(iTheta)],
        std::span<const Real>(_impl->sqrtInt),
        std::span<const Real>(_impl->sqrtIntInv));
    return ConstGSHView<Real, _MRange>(lMax, lMax, n, scratch.data());
  }

  // Where a generating grid puts the block it has just computed.
  //
  // Per thread and grow-only, like the accumulator and the coefficient
  // scratch, and for the same reason: allocating it per colatitude would put
  // an allocation back into every call. It is 528 KB at lMax = 256
  // against a table of 648 MB, and 8.4 MB at lMax = 1024 against 43 GB.
  static std::vector<Real>& WignerScratch(std::size_t size) {
    thread_local auto buffer = std::vector<Real>{};
    if (buffer.size() < size) buffer.resize(size);
    return buffer;
  }

  // Scratch for one chunk's coefficients in [coefficient][field] order.
  //
  // Kept per thread and only grown, for the same reason the work buffers and
  // the accumulator are: allocating it per call would put an allocation back
  // into every call. The forward transform accumulates into it and
  // the inverse gathers into it.
  static std::vector<Complex>& CoefficientScratch(std::size_t size) {
    thread_local auto buffer = std::vector<Complex>{};
    if (buffer.size() < size) buffer.resize(size);
    return buffer;
  }

  // What must be true of a pair of batch descriptors before anything is read
  // through them.
  static void CheckBatch(Batch in, Batch out, Int inSize, Int outSize) {
    if (in.Count() != out.Count()) {
      throw std::invalid_argument(
          "Transform batch counts differ: " + std::to_string(in.Count()) +
          " fields in and " + std::to_string(out.Count()) + " out");
    }
    // The transform writes every element of every field it is given, so
    // members that overlap would produce a wrong answer rather than an error.
    if (!in.Disjoint(inSize) || !out.Disjoint(outSize)) {
      throw std::invalid_argument(
          "Transform batch members overlap at this size");
    }
  }

  static void CheckSpan(std::size_t given, std::integral auto needed,
                        const char* what) {
    if (given < static_cast<std::size_t>(needed)) {
      throw std::invalid_argument(std::string("Transform ") + what +
                                  " range has size " + std::to_string(given) +
                                  ", but this batch spans " +
                                  std::to_string(needed));
    }
  }

  // The coefficient count for a field of the given scalar type: reduced
  // m >= 0 storage for a real field, all orders for a complex one.
  //
  // Named distinctly from CoefficientSize above rather than overloading it,
  // since the two are selected by different things -- that one by the upper
  // index, this one by the scalar -- and an overload set spanning both would
  // be resolved by argument count alone.
  template <RealOrComplexFloatingPoint Scalar>
  auto CoefficientSizeFor(Int lMax, Int n) const {
    if constexpr (RealFloatingPoint<Scalar>) {
      return GSHIndices<NonNegative>(lMax, lMax, n).Size();
    } else {
      return GSHIndices<All>(lMax, lMax, n).Size();
    }
  }

  // Size mismatches were assert-only, so under NDEBUG a short output range was
  // a silent heap overflow. Checked in all build modes.
  static void CheckSize(std::size_t given, std::integral auto expected,
                        const char* what) {
    if (given != static_cast<std::size_t>(expected)) {
      throw std::invalid_argument(std::string("Transform ") + what +
                                  " range has size " + std::to_string(given) +
                                  ", but this request needs " +
                                  std::to_string(expected));
    }
  }

  // The longitude quadrature is the trapezoid rule on nPhi equally spaced
  // points, which is exact for exp(i (m - m') phi) only when |m - m'| < nPhi.
  // Resolving orders |m| <= lMax therefore needs nPhi >= 2 * lMax + 1, not
  // 2 * lMax: at 2 * lMax the orders m = +lMax and m = -lMax are the same
  // discrete mode and cannot be separated, which is why the transform used to
  // zero the (lMax, lMax) coefficient rather than compute it. The smallest
  // fast FFT length at or above the bound is used, so that
  // the fix does not land on a length with a large prime factor.
  auto NPhi() const { return _impl->nPhi; }

  template <RealOrComplexFloatingPoint Scalar>
  void ValidateTransformRequest(Int lMax, Int n) const {
    if (lMax < 0 || lMax > MaxDegree()) {
      throw std::invalid_argument(
          "Transform degree must be between zero and the grid maximum degree");
    }
    if (std::abs(n) > lMax || !std::ranges::contains(this->UpperIndices(), n)) {
      throw std::invalid_argument(
          "Transform upper index is not supported at the requested degree");
    }

    // A spin-weighted field of nonzero upper index cannot be real-valued:
    // real-valuedness is not preserved by the frame rotation
    // e_{+-} -> e^{-+ i psi} e_{+-}, so it is not a property any component of
    // any tensor can have at N != 0; see section 7, item 5, of the theory
    // note, docs/canonical-components.tex. The
    // reduced m >= 0 coefficient storage that a real transform uses assumes
    // the self-relation f^N_{l,-m} = (-1)^{m-N} conj(f^N_{lm}), which holds
    // only when f is its own conjugate, i.e. only at N = 0. n is a runtime
    // argument, so this is a throw rather than a static_assert.
    if constexpr (RealFloatingPoint<Scalar>) {
      if (n != 0) {
        throw std::invalid_argument(
            "Real-valued fields exist only at upper index zero");
      }
    }
  }

  // Everything a grid owns *that decides its table*, built once and never
  // mutated afterwards. Shared by every copy of the handle, which is what
  // makes copying cheap and what makes concurrent use safe: an immutable
  // object behind a shared_ptr needs no synchronisation.
  //
  // The planner flag and the chunking policy are deliberately **not** here.
  // They are read per call and neither decides the table,
  // so keeping them beside a 648 MB object meant that changing either
  // rebuilt it -- and sweeping four candidate chunks, which is what a tuner
  // and the benchmark harness both do, built four tables to choose an
  // integer. They live on the handle instead, where With() can change one
  // for the price of a pointer copy.
  struct Impl {
    Impl(Int lMaxIn, Int nMaxIn, std::vector<Real> coLatitudesIn,
         std::vector<Real> coLatitudeWeightsIn, Int nPhiIn,
         WignerValues valuesIn, TransformKernel kernelIn)
        : lMax{lMaxIn},
          nMax{nMaxIn},
          values{valuesIn},
          kernel{kernelIn},
          coLatitudes{std::move(coLatitudesIn)},
          coLatitudeWeights{std::move(coLatitudeWeightsIn)},
          nPhi{nPhiIn} {
      assert(lMax >= 0);
      assert(std::abs(nMax) <= lMax);

      // The derived grid's contract, checked rather than trusted: a grid
      // that gets the nodes wrong otherwise fails inside the Wigner
      // recursion, where the message would be about something else
      // entirely.
      if (coLatitudes.size() != coLatitudeWeights.size()) {
        throw std::invalid_argument(
            "A grid needs one quadrature weight per colatitude");
      }
      if (coLatitudes.empty()) {
        throw std::invalid_argument("A grid needs at least one colatitude");
      }
      constexpr auto pi = std::numbers::pi_v<Real>;
      if (!(coLatitudes.front() > 0) || !(coLatitudes.back() < pi)) {
        throw std::invalid_argument(
            "A grid's colatitudes must lie strictly inside (0, pi)");
      }
      for (std::size_t i = 1; i < coLatitudes.size(); ++i) {
        if (!(coLatitudes[i] > coLatitudes[i - 1])) {
          throw std::invalid_argument(
              "A grid's colatitudes must be strictly increasing");
        }
      }
      if (nPhi < 1) {
        throw std::invalid_argument("A grid needs at least one longitude");
      }

      // An MRange = NonNegative grid stores only m >= 0, so it cannot serve a
      // complex-valued transform at all, and its real-valued transforms exist
      // only at upper index zero. Such a grid with nMax != 0 could serve no
      // transform whatever, so it is a configuration error rather than a
      // wasteful but usable choice: this is the "real scalar grid" reading
      // of MRange.
      if constexpr (std::same_as<_MRange, NonNegative>) {
        if (nMax != 0) {
          throw std::invalid_argument(
              "A grid storing only non-negative orders serves real-valued "
              "transforms at upper index zero, so its maximum upper index "
              "must be zero");
        }
      }

      // The matrix kernel needs values in an order the recursion cannot
      // produce one order at a time, so the two policies do not compose. See
      // Policies.h at TransformKernel: this is a fact about the recursion and
      // not an unimplemented case, which is why it is
      // refused here rather than worked around.
      // BLAS offers single and double precision and nothing wider, so a
      // long double grid cannot have the matrix kernel whatever else is true.
      // Refused here, where Real is known, rather than in a call.
      //
      // Guarded, because BlasDetails does not exist without a BLAS -- and
      // without one neither does TransformKernel::Matrix(), so there is
      // nothing left to refuse.
#ifdef GSHTRANS_HAVE_BLAS
      if constexpr (!BlasDetails::BlasReal<Real>) {
        if (kernel.IsMatrix()) {
          throw std::invalid_argument(
              "The matrix transform kernel needs a BLAS, and BLAS offers "
              "single and double precision only, so it is unavailable at this "
              "grid's precision");
        }
      }
#endif

      if (kernel.IsMatrix() && !values.AreStored()) {
        throw std::invalid_argument(
            "The matrix transform kernel cannot generate its Wigner values on "
            "the fly: the recursion produces every order of one colatitude "
            "together, so a single (n, m) block cannot be had without either "
            "keeping the whole table or repeating the recursion for every "
            "order. Ask for WignerValues::Stored(), or for "
            "TransformKernel::Loop()");
      }

      // A generating grid builds no table. What it needs instead is the two
      // square-root tables the recursion indexes, which are 2 lMax + 1
      // entries each against the table's 648 MB at lMax = 256.
      //
      // A matrix grid builds the same values in the transform-major layout
      // instead, and only that one: the two are the same size and holding
      // both would double 648 MB for nothing.
      if (!values.AreStored()) {
        std::tie(sqrtInt, sqrtIntInv) =
            WignerDetails::PreComputeTables<Real>(lMax, lMax, nMax);
      } else if (kernel.IsMatrix()) {
        // Reflected: non-negative orders only, halving 648 MB to 324 at
        // lMax = 256 with nMax = 2. What makes it available is that the nodes
        // are symmetric about the equator, and WignerMatrices checks that
        // rather than taking it on trust -- so a grid whose nodes are not
        // symmetric simply does not get the halved table. That is a
        // property of the node set and not of Gauss-Legendre, which is why
        // the check lives here.
        wignerMatrices = WignerMatrices<Real, _MRange, _NRange>::Reflected(
            lMax, lMax, nMax, coLatitudes);
      } else {
        wigner = Wigner<Real, _MRange, _NRange, Multiple>(lMax, lMax, nMax,
                                                          coLatitudes);
      }

      // The planner flag is kept as the caller gave it. It used to be used to
      // pre-generate wisdom for exactly two shapes and then replaced by
      // WisdomOnly, which meant that any shape the constructor had not
      // anticipated -- every batched shape, in particular -- would fail to
      // plan rather than fall back. Shapes are planned on first use and
      // cached, so there is nothing to anticipate.
    }

    Int lMax;
    Int nMax;
    WignerValues values;
    TransformKernel kernel;
    std::vector<Real> coLatitudes;
    std::vector<Real> coLatitudeWeights;
    Int nPhi;

    // Empty on a generating grid, which is the whole of what that grid saves,
    // and on a matrix grid, which holds the other layout instead.
    std::optional<Wigner<Real, _MRange, _NRange, Multiple>> wigner;

    // Empty unless this is a matrix grid. Exactly one of these two is ever
    // occupied, and on a generating grid neither is.
    std::optional<WignerMatrices<Real, _MRange, _NRange>> wignerMatrices;

    // Empty on a stored grid, whose table already carries what these are for.
    std::vector<Real> sqrtInt;
    std::vector<Real> sqrtIntInv;
  };

  std::shared_ptr<const Impl> _impl;

  // Read per call and shared with nothing. See Impl above for why they sit
  // here rather than in it.
  //
  // The flag is what an uncached plan shape is planned with. It was once
  // followed by the constructor generating wisdom and then setting
  // WisdomOnly, which meant that any shape the constructor had not
  // anticipated -- every batched shape, in particular -- would fail to plan
  // rather than fall back. Shapes are planned on first use and cached, so
  // there is nothing to anticipate.
  Chunking _chunking;
  FFTWpp::Flag _flag;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SPHERICAL_GRID_GUARD_H
