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
#include <numeric>
#include <map>
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

  GaussLegendreGrid(Int lMax, Int nMax, FFTWpp::Flag flag = FFTWpp::Measure)
      : _impl{std::make_shared<const Impl>(lMax, nMax, flag)} {}

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
                      FFTWpp::Flag flag = FFTWpp::Measure) {
    if (lBand < 0) {
      throw std::invalid_argument("Band must be non-negative");
    }
    if (!(oversampling >= 1)) {
      throw std::invalid_argument("Oversampling factor must be at least one");
    }
    const auto lGrid = static_cast<Int>(
        std::ceil(oversampling * static_cast<Real>(lBand)));
    return GaussLegendreGrid(lGrid, nMax, flag);
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
  void ForwardTransformation(Int lMax, Int n, InRange&& in, OutRange& out,
                             Execution policy = Execution::Sequential()) const {
    // Get scalar type for field.
    using Scalar = std::ranges::range_value_t<InRange>;

    ValidateTransformRequest<Scalar>(lMax, n);
    CheckSize(std::ranges::size(in), this->FieldSize(), "field");
    CheckSize(std::ranges::size(out), CoefficientSizeFor<Scalar>(lMax, n),
              "coefficient");

    // The colatitude loop below is the quadrature sum, so it accumulates into
    // out. Nothing else initialises it, and every caller happened to arrive
    // with a zeroed buffer, which made "out arrives zeroed" an unstated and
    // unchecked precondition: transforming twice into the same buffer doubled
    // the answer (core-plan.md F1). Zero it here and the routine owns its
    // output. Accumulation, if ever wanted, is a separate named entry point.
    std::ranges::fill(out, Complex{});

    // A one-point grid needs no FFT.
    if (_impl->lMax == 0) {
      out[0] =
          in[0] * static_cast<Real>(2) / std::numbers::inv_sqrtpi_v<Real>;
      return;
    }

    // Pre compute some constants.
    const auto nPhi = this->NumberOfLongitudes();
    const auto scaleFactor = static_cast<Real>(2) * std::numbers::pi_v<Real> /
                             static_cast<Real>(nPhi);

    // One colatitude's contribution, accumulated wherever the caller says.
    // Taking the destination as a parameter is what lets the parallel path
    // give each thread a private accumulator without a second copy of the
    // loop.
    auto AccumulateRow = [&](Int iTheta, auto outBegin, auto& work) {
      // FFT the current data slice, through the plan's own buffer.
      PackRow(std::next(in.begin(), iTheta * nPhi), nPhi, work.in);
      work.plan.Execute();

      // Get the Wigner values and quadrature weight.
      auto d = _impl->wigner[n, iTheta];
      auto w = _impl->quad.W(iTheta) * scaleFactor;

      // Loop over the spherical harmonic coefficients, taking the Wigner
      // values one degree at a time.
      //
      // The row pointer comes from d[l] rather than from a single iterator
      // walked across the whole block. That is the supplier seam of
      // core-plan.md [C10]: the only thing this loop needs is a contiguous run
      // of values in (l, m) order, one run per degree, and asking for it per
      // degree is what lets step F' substitute a supplier that generates the
      // row into per-thread scratch for one that points into the stored table.
      // OffsetForDegree is closed-form, so the seam costs a few integer
      // operations per degree against a loop of length 2l+1.
      //
      // The two constraints it carries, both already satisfied here: degrees
      // are visited in ascending contiguous order from |n|, and each
      // (n, iTheta) is visited once per pass. A generated row cannot be
      // revisited without re-running the recursion.
      auto outIter = outBegin;
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });

      for (auto l : degrees) {
        auto dl = d[l];
        auto wigIter = dl.begin();

        if constexpr (ComplexFloatingPoint<Scalar>) {
          auto workIter = std::prev(work.out.end(), dl.MaxOrder());
          for (auto m : dl.NegativeOrders()) {
            *outIter++ += *wigIter++ * *workIter++ * w;
          }
          workIter = work.out.begin();
          for (auto m : dl.NonNegativeOrders()) {
            *outIter++ += *wigIter++ * *workIter++ * w;
          }
        } else {
          auto workIter = work.out.begin();
          if constexpr (std::same_as<_MRange, All>) {
            std::advance(wigIter, dl.MaxOrder());
          }
          for (auto m : dl.NonNegativeOrders()) {
            *outIter++ += *wigIter++ * *workIter++ * w;
          }
        }
      }
    };

    const auto nTheta = static_cast<Int>(this->NumberOfCoLatitudes());

    if (!RunInParallel(policy)) {
      auto& work = GetWorkspace<Scalar, true>(nPhi, _impl->flag);
      for (auto iTheta = Int{0}; iTheta < nTheta; iTheta++) {
        AccumulateRow(iTheta, out.begin(), work);
      }
      return;
    }

    // Every colatitude contributes to every coefficient, so the colatitudes
    // cannot simply be divided between threads writing into `out`. Each thread
    // accumulates into a private buffer and the partial sums are added at the
    // end.
    //
    // That addition is *partitioned*, not serialised. Each thread owns one
    // block of the coefficient array and sums every thread's partial sums for
    // that block alone, so the reduction runs in parallel and no thread waits
    // on another. It used to be a critical section in which each thread added
    // a whole coefficient array in turn, which costs one serialised pass per
    // thread: invisible against the colatitude loop at eight threads, and
    // 128 MB of serialised adds per transform at 128 (core-plan.md P8, and
    // step H's own suggested fix). The decomposition itself is unchanged --
    // thread-private accumulators over colatitudes are the wrong shape well
    // before 128 threads, but choosing what replaces them needs a machine this
    // was not measured on ([C11]).
    const auto size = static_cast<std::size_t>(std::ranges::size(out));
    const auto threadCount = ThreadCount(policy);

    // The reduction reads every thread's accumulator, so the thread-local
    // buffers have to be published to the team. Written before the implicit
    // barrier at the end of the colatitude loop and read after it.
    auto partials = std::vector<Complex*>(threadCount, nullptr);

#pragma omp parallel num_threads(threadCount)
    {
      const auto thread = static_cast<std::size_t>(omp_get_thread_num());
      const auto threads = static_cast<std::size_t>(omp_get_num_threads());
      auto& work = GetWorkspace<Scalar, true>(nPhi, _impl->flag);
      auto& partial = Accumulator(size);
      std::fill_n(partial.begin(), size, Complex{});
      partials[thread] = partial.data();

#pragma omp for schedule(static)
      for (Int iTheta = 0; iTheta < nTheta; iTheta++) {
        AccumulateRow(iTheta, partial.begin(), work);
      }

      const auto first = size * thread / threads;
      const auto last = size * (thread + 1) / threads;
      const auto outFirst = std::next(std::ranges::begin(out), first);
      for (auto t = std::size_t{0}; t < threads; t++) {
        const auto* p = partials[t];
        auto outIter = outFirst;
        for (auto i = first; i < last; i++) *outIter++ += p[i];
      }
    }
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
  void InverseTransformation(Int lMax, Int n, InRange&& in, OutRange& out,
                             Execution policy = Execution::Sequential()) const {
    // Get scalar type for field.
    using Scalar = std::ranges::range_value_t<OutRange>;

    ValidateTransformRequest<Scalar>(lMax, n);
    CheckSize(std::ranges::size(in), CoefficientSizeFor<Scalar>(lMax, n),
              "coefficient");
    CheckSize(std::ranges::size(out), this->FieldSize(), "field");

    // A one-point grid needs no FFT.
    if (_impl->lMax == 0) {
      if constexpr (RealFloatingPoint<Scalar>) {
        out[0] = std::real(in[0]) * std::numbers::inv_sqrtpi_v<Real> /
                 static_cast<Real>(2);
      } else {
        out[0] = in[0] * std::numbers::inv_sqrtpi_v<Real> /
                 static_cast<Real>(2);
      }
      return;
    }

    // Precompute constants
    const auto nPhi = this->NumberOfLongitudes();

    // One colatitude, synthesised into its own row of the field. Unlike the
    // forward transform, the colatitudes here write disjoint output and share
    // only read-only input, so threading over them needs no reduction.
    auto SynthesiseRow = [&](Int iTheta, auto& work) {
      std::ranges::for_each(work.in, [](auto& x) { return x = 0; });

      // Get the Wigner values.
      auto d = _impl->wigner[n, iTheta];

      // Loop over the coefficients, one degree at a time. As in the forward
      // direction, the row pointer comes from d[l]: this is the same supplier
      // seam, and step F' substitutes at the same point.
      auto inIter = in.begin();
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });
      for (auto l : degrees) {
        auto dl = d[l];
        auto wigIter = dl.begin();
        if constexpr (ComplexFloatingPoint<Scalar>) {
          auto workIter = std::prev(work.in.end(), dl.MaxOrder());
          for (auto m : dl.NegativeOrders()) {
            *workIter++ += *inIter++ * *wigIter++;
          }
          workIter = work.in.begin();
          for (auto m : dl.NonNegativeOrders()) {
            *workIter++ += *inIter++ * *wigIter++;
          }
        } else {
          auto workIter = work.in.begin();
          if constexpr (std::same_as<_MRange, All>) {
            std::advance(wigIter, dl.MaxOrder());
          }
          for (auto m : dl.NonNegativeOrders()) {
            *workIter++ += *inIter++ * *wigIter++;
          }
        }
      }

      // Perform FFT to recover field at the colatitude, through the plan's
      // own buffer, then hand the row to the caller.
      work.plan.Execute();
      UnpackRow(work.out, nPhi, std::next(out.begin(), iTheta * nPhi));
    };

    const auto nTheta = static_cast<Int>(this->NumberOfCoLatitudes());

    if (!RunInParallel(policy)) {
      auto& work = GetWorkspace<Scalar, false>(nPhi, _impl->flag);
      for (auto iTheta = Int{0}; iTheta < nTheta; iTheta++) {
        SynthesiseRow(iTheta, work);
      }
      return;
    }

#pragma omp parallel num_threads(ThreadCount(policy))
    {
      auto& work = GetWorkspace<Scalar, false>(nPhi, _impl->flag);
#pragma omp for schedule(static)
      for (Int iTheta = 0; iTheta < nTheta; iTheta++) {
        SynthesiseRow(iTheta, work);
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

    Workspace(Int nPhi, FFTWpp::Flag flag)
        : in(FFTWpp::DataSize<In, Out>(nPhi).first),
          out(FFTWpp::DataSize<In, Out>(nPhi).second),
          plan(MakePlan(in, out, flag)) {}

    static auto MakePlan(FFTWpp::vector<In>& in, FFTWpp::vector<Out>& out,
                         FFTWpp::Flag flag) {
      auto inView = FFTWpp::Ranges::View(in);
      auto outView = FFTWpp::Ranges::View(out);
      auto lock = std::scoped_lock(PlannerMutex());
      if constexpr (std::same_as<In, Out>) {
        return FFTWpp::Ranges::Plan(
            inView, outView, flag,
            IsForward ? FFTWpp::Forward : FFTWpp::Backward);
      } else {
        return FFTWpp::Ranges::Plan(inView, outView, flag);
      }
    }

    FFTWpp::vector<In> in;
    FFTWpp::vector<Out> out;
    decltype(MakePlan(std::declval<FFTWpp::vector<In>&>(),
                      std::declval<FFTWpp::vector<Out>&>(),
                      std::declval<FFTWpp::Flag>())) plan;
  };

  // Plans and buffers used to be created on every call -- two allocations and
  // a plan per transform, with FFTW's non-re-entrant planner run each time
  // (P7). They are now made once per thread per shape and kept. Held by
  // pointer so that the plan's reference to its buffers survives any
  // rehashing of the cache.
  template <RealOrComplexFloatingPoint Scalar, bool IsForward>
  static Workspace<Scalar, IsForward>& GetWorkspace(Int nPhi,
                                                    FFTWpp::Flag flag) {
    using Entry = Workspace<Scalar, IsForward>;
    thread_local auto cache =
        std::map<std::pair<Int, unsigned>, std::unique_ptr<Entry>>{};
    const auto key = std::pair{nPhi, static_cast<unsigned>(flag)};
    auto found = cache.find(key);
    if (found == cache.end()) {
      found = cache.emplace(key, std::make_unique<Entry>(nPhi, flag)).first;
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
  template <typename Iterator, typename Buffer>
  static void PackRow(Iterator first, Int nPhi, Buffer& work) {
    std::copy_n(first, nPhi, work.begin());
  }

  template <typename Buffer, typename Iterator>
  static void UnpackRow(const Buffer& work, Int nPhi, Iterator first) {
    std::copy_n(work.begin(), nPhi, first);
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
    Impl(Int lMaxIn, Int nMaxIn, FFTWpp::Flag flagIn)
        : lMax{lMaxIn}, nMax{nMaxIn}, flag{flagIn} {
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

      wigner = Wigner<Real, _MRange, _NRange, Multiple, ColumnMajor>(
          lMax, lMax, nMax, quad.Points());

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
    GaussQuad::Quadrature1D<Real> quad;
    Wigner<Real, _MRange, _NRange, Multiple, ColumnMajor> wigner;
  };

  std::shared_ptr<const Impl> _impl;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
