#ifndef GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
#define GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H

#include <omp.h>

#include <FFTWpp/Core>
#include <FFTWpp/Ranges>
#include <GaussQuad/All>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <numeric>
#include <array>
#include <map>
#include <tuple>
#include <mutex>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <NumericConcepts/Ranges.hpp>

#include "Concepts.h"
#include "GridBase.h"
#include "Indexing.h"
#include "Utility.h"
#include "Wigner.h"

namespace GSHTrans {

template <RealFloatingPoint _Real, OrderIndexRange _MRange, IndexRange _NRange>
class GaussLegendreGrid
    : public GridBase<GaussLegendreGrid<_Real, _MRange, _NRange>> {
 public:
  // Public type aliases.
  using Int = std::ptrdiff_t;
  using Real = _Real;
  using Complex = std::complex<Real>;
  using MRange = _MRange;
  using NRange = _NRange;

  // A grid is a value-semantic handle over an immutable, shared
  // implementation: constructing one builds the quadrature and the Wigner
  // table, and copying one copies a pointer.
  //
  // That is not a convenience. The Wigner table at lMax = 256, nMax = 2 with
  // all upper indices is about 679 MB, against 2.1 MB for a complex field on
  // the same grid, and the members used to be held by value with defaulted
  // copy: copying a grid by accident was not a performance wart but an
  // out-of-memory event (core-plan.md F9). Putting the indirection inside the
  // grid rather than leaving each consumer to wrap it in a shared_ptr also
  // makes the question unaskable, and gives step E's plan cache somewhere to
  // live that is shared by construction rather than by convention.
  //
  // There is no default constructor: a grid without a quadrature is not a
  // grid, and reading MaxDegree() on a default-constructed one was undefined
  // (F4).
  GaussLegendreGrid() = delete;

  // The chunking policy is a property of the machine rather than of the call,
  // which is why it is set here alongside the planner flag and not on every
  // transform. `Automatic` assumes a modest cache; a caller who knows their
  // machine passes `Chunking::ForCache(bytes)`, and one who has measured
  // their own optimum passes `Chunking::Fixed(count)`.
  GaussLegendreGrid(Int lMax, Int nMax, FFTWpp::Flag flag = FFTWpp::Measure,
                    Chunking chunking = Chunking::Automatic(),
                    WignerValues values = WignerValues::Stored())
      : _impl{std::make_shared<const Impl>(lMax, nMax, flag, chunking,
                                           values)} {}

  // A grid for working with fields of maximum degree lBand, with quadrature
  // headroom for degree oversampling * lBand.
  //
  // The distinction this expresses is between the *band* of the fields being
  // worked with and the *resolution* of the grid, which the library could not
  // previously say: the transforms already take lMax per call, so an
  // oversampled grid with truncated transforms works, but there was no way to
  // ask for one. Headroom is wanted whenever a quantity of higher degree than
  // the fields themselves is formed on the grid -- a product of two band-L
  // fields has band 2L, and integrating |f|^2 for band-L f integrates a
  // degree-2L quantity. The 3/2 rule is oversampling = 1.5; oversampling = 2
  // is exact for a single product.
  static auto ForBand(Int lBand, Int nMax, Real oversampling = 1,
                      FFTWpp::Flag flag = FFTWpp::Measure,
                      Chunking chunking = Chunking::Automatic(),
                      WignerValues values = WignerValues::Stored()) {
    if (lBand < 0) {
      throw std::invalid_argument("Band must be non-negative");
    }
    if (!(oversampling >= 1)) {
      throw std::invalid_argument("Oversampling factor must be at least one");
    }
    const auto lGrid = static_cast<Int>(
        std::ceil(oversampling * static_cast<Real>(lBand)));
    return GaussLegendreGrid(lGrid, nMax, flag, chunking, values);
  }

  GaussLegendreGrid(const GaussLegendreGrid&) = default;

  GaussLegendreGrid(GaussLegendreGrid&&) = default;

  GaussLegendreGrid& operator=(const GaussLegendreGrid&) = default;

  GaussLegendreGrid& operator=(GaussLegendreGrid&&) = default;

  // Two grids are the same grid when they share an implementation. This is
  // the test the field layer's binary nodes make at construction. Structural
  // comparison is deliberately not offered: two separately built grids with
  // equal parameters hold different quadrature objects and different wisdom,
  // and treating them as interchangeable would make a node's operands
  // silently disagree about the buffers they index.
  auto Identity() const { return _impl.get(); }

  //------------------------------------------------//
  //    Methods needed to inherit from GridBase     //
  //------------------------------------------------//
  auto MaxDegree() const { return _impl->lMax; }
  auto MaxUpperIndex() const { return _impl->nMax; }

  auto CoLatitudes() const { return std::ranges::views::all(_impl->quad.Points()); }
  auto CoLatitudeWeights() const {
    return std::ranges::views::all(_impl->quad.Weights());
  }

  auto Longitudes() const {
    const auto nPhi = NPhi();
    const auto dPhi = 2 * std::numbers::pi_v<Real> / static_cast<Real>(nPhi);
    return std::ranges::views::iota(Int{0}, nPhi) |
           std::ranges::views::transform([dPhi](auto i) { return i * dPhi; });
  }
  auto LongitudeWeights() const {
    const auto nPhi = NPhi();
    const auto dPhi = 2 * std::numbers::pi_v<Real> / static_cast<Real>(nPhi);
    return std::ranges::views::repeat(dPhi, nPhi);
  }

  //-----------------------------------------------------//
  //          Forward transformation for ranges          //
  //-----------------------------------------------------//
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
      const auto w = _impl->quad.W(iTheta) * scaleFactor;
      const auto orders = static_cast<Int>(work.out.size()) / c;

      // Loop over the spherical harmonic coefficients, taking the Wigner
      // values one degree at a time.
      //
      // The row pointer comes from d[l] rather than from a single iterator
      // walked across the whole block. That is the supplier seam of
      // core-plan.md [C10]: the only thing this loop needs is a contiguous
      // run of values in (l, m) order, one run per degree, and asking for it
      // per degree is what lets step F' substitute a supplier that generates
      // the row into per-thread scratch for one that points into the stored
      // table. OffsetForDegree is closed-form, so the seam costs a few
      // integer operations per degree against a loop of length (2l+1) * c.
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
    // initialised the destination (core-plan.md F1). Zeroing the range first
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
        auto& work = GetWorkspace<Scalar, true>(nPhi, c, _impl->flag);
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
      // eight threads, and 128 MB of serialised adds per transform at 128
      // (core-plan.md P8, and step H's own suggested fix). The decomposition
      // itself is unchanged -- thread-private accumulators over colatitudes
      // are the wrong shape well before 128 threads, but choosing what
      // replaces them needs a machine this was not measured on ([C11]).
      const auto threadCount = ThreadCount(policy);

      // The reduction reads every thread's accumulator, so the thread-local
      // buffers have to be published to the team. Written before the implicit
      // barrier at the end of the colatitude loop and read after it.
      auto partials = std::vector<Complex*>(threadCount, nullptr);

#pragma omp parallel num_threads(threadCount)
      {
        const auto thread = static_cast<Int>(omp_get_thread_num());
        const auto threads = static_cast<Int>(omp_get_num_threads());
        auto& work = GetWorkspace<Scalar, true>(nPhi, c, _impl->flag);
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

  // The single field, which is the batched primitive at count = 1 ([C3]).
  //
  // The size checks here are equalities rather than spans, because this
  // entry's contract is that the range *is* the field: a caller who hands
  // over a range of the wrong length has made a mistake, whereas a batched
  // caller may legitimately hand over a window onto a larger one.
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
      // seam, and step F' substitutes at the same point. The inner loop is
      // again an axpy of length c, over a batch index that runs fastest on
      // both sides.
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
        UnpackRow(std::next(work.out.begin(), k * nPhi), nPhi,
                  std::next(outFirst,
                            outBatch.Offset(iTheta * nPhi, first + k)),
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
      // ([C9]). It costs one pass over the coefficients against a colatitude
      // loop that reads the whole Wigner block, which T4 measured to be
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
        auto& work = GetWorkspace<Scalar, false>(nPhi, c, _impl->flag);
        for (auto iTheta = Int{0}; iTheta < nTheta; iTheta++) {
          SynthesiseRow(iTheta, first, c, gathered, work);
        }
        continue;
      }

#pragma omp parallel num_threads(ThreadCount(policy))
      {
        auto& work = GetWorkspace<Scalar, false>(nPhi, c, _impl->flag);
#pragma omp for schedule(static)
        for (Int iTheta = 0; iTheta < nTheta; iTheta++) {
          SynthesiseRow(iTheta, first, c, gathered, work);
        }
      }
    }
  }

  // The single field, which is the batched primitive at count = 1 ([C3]).
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

  // A per-thread accumulator for the forward transform's partial sums, kept
  // between calls for the same reason the work buffers are: this is the size
  // of the coefficient array, and allocating it per call would put back the
  // per-call allocation step E removed. It only ever grows.
  static std::vector<Complex>& Accumulator(std::size_t size) {
    thread_local auto buffer = std::vector<Complex>{};
    if (buffer.size() < size) buffer.resize(size);
    return buffer;
  }

  // FFTW's planner is not re-entrant, so plan creation everywhere in the
  // process is serialised on this. Execution is not: FFTW does allow a plan to
  // be executed concurrently, and in any case each thread executes on its own
  // workspace.
  static std::mutex& PlannerMutex() {
    static std::mutex mutex;
    return mutex;
  }

  // The aligned buffers a transform works in, and the FFTW plan bound to them.
  //
  // Held per thread rather than shared. The buffers are scratch and must be
  // per-thread whatever else happens; binding the plan to them keeps the two
  // together, lets execution use the plan's own buffers rather than
  // FFTW's new-array form, and leaves Impl immutable, which is what makes a
  // shared grid safe to use concurrently without a lock (core-plan.md step B).
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

    // The layouts are the whole of P6's "batching the FFT is free": FFTWpp
    // passes rank, howMany, embed, stride and dist straight through to
    // plan_many, so the batch costs a descriptor rather than a repack.
    //
    // The asymmetry between the two sides is deliberate. On the *field* side
    // the rows are packed one after another, because that is what a copy from
    // the caller's storage produces however the caller strides. On the
    // *coefficient* side the batch index runs fastest, so the data lands in
    // [m][k] order and the Legendre stage's inner loop over the batch is
    // unit-stride. That layout costs nothing to ask for here and is the one
    // P2's batching measurement was made on.
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
          IsForward
              ? FFTWpp::Ranges::Layout(1, inN, howMany, inN, 1, inSize)
              : FFTWpp::Ranges::Layout(1, inN, howMany, inN, howMany, 1);
      auto outLayout =
          IsForward
              ? FFTWpp::Ranges::Layout(1, outN, howMany, outN, howMany, 1)
              : FFTWpp::Ranges::Layout(1, outN, howMany, outN, 1, outSize);

      auto inView = FFTWpp::Ranges::View(in, inLayout);
      auto outView = FFTWpp::Ranges::View(out, outLayout);
      auto lock = std::scoped_lock(PlannerMutex());
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
  // (P7). They are now made once per thread per shape and kept. Held by
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
  // planning buffers, and neither FFTW nor FFTWpp checks (core-plan.md F3).
  // The forward transform used to take that path on any writable input and
  // the inverse took it unconditionally on the caller's output, which made
  // alignment an obligation propagating outward into every stride the field
  // and tensor layers might choose. Copying instead costs one pass over a row
  // against an FFT of the same row, and it is what makes the field plan's
  // promise -- that a slice target needs no alignment beyond Scalar's --
  // true rather than aspirational.
  //
  // They are also the seam step F widens: the (stride, dist) of a batch
  // ([C9]) enters here and nowhere else.
  // The caller's stride enters at these two points and at no other. A batch
  // may be interleaved with components this call knows nothing about
  // ([C9]), and everything downstream of here works on our own contiguous
  // buffers.
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
  // measured 2.2x faster (core-plan.md section 10). The forward's rule is
  // unchanged, and the same measurement is the evidence for that too -- the
  // whole batch there is 2x *slower*, because eight private accumulators of
  // 8.4 MB ask for 68 MB of a 16 MB cache.
  //
  // The thread count is one whenever the policy is sequential or a parallel
  // region is already open, since that is what the call will really run on.
  Int ForwardChunkSize(Int coefficientSize, Execution policy) const {
    const auto threads = RunInParallel(policy) ? ThreadCount(policy) : 1;
    return _impl->chunking.Count(
        coefficientSize * static_cast<Int>(sizeof(Complex)), threads);
  }

  Int InverseChunkSize(Int coefficientSize) const {
    return _impl->chunking.Count(
        coefficientSize * static_cast<Int>(sizeof(Complex)), 1);
  }

  // The Wigner values for one (n, iTheta), over the degrees |n| .. lMax.
  //
  // This is the seam of core-plan.md [C10], and it returns the same type on
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
        _impl->quad.X(static_cast<int>(iTheta)),
        std::span<const Real>(_impl->sqrtInt),
        std::span<const Real>(_impl->sqrtIntInv));
    return ConstGSHView<Real, _MRange>(lMax, lMax, n, scratch.data());
  }

  // Where a generating grid puts the block it has just computed.
  //
  // Per thread and grow-only, like the accumulator and the coefficient
  // scratch, and for the same reason: allocating it per colatitude would put
  // back the per-call allocation step E removed. It is 528 KB at lMax = 256
  // against a table of 648 MB, and 8.4 MB at lMax = 1024 against 43 GB.
  static std::vector<Real>& WignerScratch(std::size_t size) {
    thread_local auto buffer = std::vector<Real>{};
    if (buffer.size() < size) buffer.resize(size);
    return buffer;
  }

  // Scratch for one chunk's coefficients in [coefficient][field] order.
  //
  // Kept per thread and only grown, for the same reason the work buffers and
  // the accumulator are: allocating it per call would put back the per-call
  // allocation step E removed. The forward transform accumulates into it and
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
      throw std::invalid_argument(
          std::string("Transform ") + what + " range has size " +
          std::to_string(given) + ", but this batch spans " +
          std::to_string(needed));
    }
  }

  // The coefficient count for a field of the given scalar type: reduced
  // m >= 0 storage for a real field, all orders for a complex one.
  // Named distinctly from GridBase::CoefficientSize, which it would otherwise
  // hide.
  template <RealOrComplexFloatingPoint Scalar>
  auto CoefficientSizeFor(Int lMax, Int n) const {
    if constexpr (RealFloatingPoint<Scalar>) {
      return GSHIndices<NonNegative>(lMax, lMax, n).Size();
    } else {
      return GSHIndices<All>(lMax, lMax, n).Size();
    }
  }

  // Size mismatches were assert-only, so under NDEBUG a short output range was
  // a silent heap overflow (core-plan.md F5). Checked in all build modes.
  static void CheckSize(std::size_t given, std::integral auto expected,
                        const char* what) {
    if (given != static_cast<std::size_t>(expected)) {
      throw std::invalid_argument(
          std::string("Transform ") + what + " range has size " +
          std::to_string(given) + ", but this request needs " +
          std::to_string(expected));
    }
  }

  // The longitude quadrature is the trapezoid rule on nPhi equally spaced
  // points, which is exact for exp(i (m - m') phi) only when |m - m'| < nPhi.
  // Resolving orders |m| <= lMax therefore needs nPhi >= 2 * lMax + 1, not
  // 2 * lMax: at 2 * lMax the orders m = +lMax and m = -lMax are the same
  // discrete mode and cannot be separated, which is why the transform used to
  // zero the (lMax, lMax) coefficient rather than compute it (core-plan.md
  // F2). The smallest fast FFT length at or above the bound is used, so that
  // the fix does not land on a length with a large prime factor.
  auto NPhi() const { return FastFFTSize(2 * _impl->lMax + 1); }

  template <RealOrComplexFloatingPoint Scalar>
  void ValidateTransformRequest(Int lMax, Int n) const {
    if (lMax < 0 || lMax > MaxDegree()) {
      throw std::invalid_argument(
          "Transform degree must be between zero and the grid maximum degree");
    }
    if (std::abs(n) > lMax ||
        !std::ranges::contains(this->UpperIndices(), n)) {
      throw std::invalid_argument(
          "Transform upper index is not supported at the requested degree");
    }

    // A spin-weighted field of nonzero upper index cannot be real-valued:
    // real-valuedness is not preserved by the frame rotation
    // e_{+-} -> e^{-+ i psi} e_{+-}, so it is not a property any component of
    // any tensor can have at N != 0 (theory note section 7, item 5). The
    // reduced m >= 0 coefficient storage that a real transform uses assumes
    // the self-relation f^N_{l,-m} = (-1)^{m-N} conj(f^N_{lm}), which holds
    // only when f is its own conjugate, i.e. only at N = 0. n is a runtime
    // argument, so this is a throw rather than a static_assert
    // (core-plan.md step A, [C2]).
    if constexpr (RealFloatingPoint<Scalar>) {
      if (n != 0) {
        throw std::invalid_argument(
            "Real-valued fields exist only at upper index zero");
      }
    }
  }

  // Everything a grid owns, built once and never mutated afterwards. Shared
  // by every copy of the handle, which is what makes copying cheap and what
  // makes concurrent use safe: an immutable object behind a shared_ptr needs
  // no synchronisation. Step E's plan cache belongs here, and will be the one
  // mutable member, with its own lock.
  struct Impl {
    Impl(Int lMaxIn, Int nMaxIn, FFTWpp::Flag flagIn, Chunking chunkingIn,
         WignerValues valuesIn)
        : lMax{lMaxIn},
          nMax{nMaxIn},
          flag{flagIn},
          chunking{chunkingIn},
          values{valuesIn} {
      assert(lMax >= 0);
      assert(std::abs(nMax) <= lMax);
      assert(flag != FFTWpp::WisdomOnly);

      // An MRange = NonNegative grid stores only m >= 0, so it cannot serve a
      // complex-valued transform at all, and its real-valued transforms exist
      // only at upper index zero. Such a grid with nMax != 0 could serve no
      // transform whatever, so it is a configuration error rather than a
      // wasteful but usable choice. See core-plan.md step A and [C1]: this is
      // the "real scalar grid" reading of MRange.
      if constexpr (std::same_as<_MRange, NonNegative>) {
        if (nMax != 0) {
          throw std::invalid_argument(
              "A grid storing only non-negative orders serves real-valued "
              "transforms at upper index zero, so its maximum upper index "
              "must be zero");
        }
      }

      quad = GaussQuad::LegendrePolynomial<Real>{}.GaussQuadrature(lMax + 1);
      quad.Transform([](auto x) { return std::acos(-x); },
                     [](auto x) -> Real { return 1; });

      // A generating grid builds no table. What it needs instead is the two
      // square-root tables the recursion indexes, which are 2 lMax + 1
      // entries each against the table's 648 MB at lMax = 256.
      if (values.AreStored()) {
        wigner = Wigner<Real, _MRange, _NRange, Multiple, ColumnMajor>(
            lMax, lMax, nMax, quad.Points());
      } else {
        std::tie(sqrtInt, sqrtIntInv) =
            WignerDetails::PreComputeTables<Real>(lMax, lMax, nMax);
      }

      // The planner flag is kept as the caller gave it. It used to be used to
      // pre-generate wisdom for exactly two shapes and then replaced by
      // WisdomOnly, which meant that any shape the constructor had not
      // anticipated -- every batched shape, in particular -- would fail to
      // plan rather than fall back (core-plan.md P7). Shapes are now planned
      // on first use and cached, so there is nothing to anticipate.
    }

    Int lMax;
    Int nMax;
    FFTWpp::Flag flag;
    Chunking chunking;
    WignerValues values;
    GaussQuad::Quadrature1D<Real> quad;

    // Empty on a generating grid, which is the whole of what that grid saves.
    std::optional<Wigner<Real, _MRange, _NRange, Multiple, ColumnMajor>> wigner;

    // Empty on a stored grid, whose table already carries what these are for.
    std::vector<Real> sqrtInt;
    std::vector<Real> sqrtIntInv;
  };

  std::shared_ptr<const Impl> _impl;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
