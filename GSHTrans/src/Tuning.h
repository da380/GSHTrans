#ifndef GSH_TRANS_TUNING_GUARD_H
#define GSH_TRANS_TUNING_GUARD_H

// Choosing a policy by measuring it, on the caller's own problem.
//
// The plan is core-plan.md section 12, and the decisions it records are
// referred to below as [C18] to [C23]. The premise, which that section argues
// from the record rather than from principle, is that some of this library's
// choices cannot be settled by reasoning and vary by machine: the chunk
// optimum "moves between runs", the batched forward gives 2.3x sequentially
// and 1.12x threaded at lMax = 256, direction-aware chunking was worth 2.0x
// and was found only by measuring. A caller on a machine neither of us has
// will otherwise run with the wrong setting and have no way to know.
//
// Three things this is not.
//
//   -- Not an autotuner. FFTW searches a space of plans it generates; this
//      times a handful of named alternatives.
//   -- Not a configured grid. Tune hands back *values* ([C19]). A tuner is
//      the piece of machinery most likely to substitute something silently --
//      wisdom naming a kernel the build cannot offer, a candidate that failed
//      to construct, a tie resolved in favour of the incumbent -- and
//      returning values makes every one of those visible in a variable the
//      caller can print. That is [C17]'s argument, and it is load-bearing
//      under [C12]: the whole justification for carrying two kernels is being
//      able to compare them, and a mechanism that reports one under the
//      other's name destroys it.
//   -- Not a promise of the optimum. There may not be one resolvable: the
//      lMax = 128 chunk peak swapped between two runs of the same binary. The
//      property a caller needs, and the one that is testable, is that tuning
//      never returns something worse than the default by more than the margin.

#include <FFTWpp/Core>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Concepts.h"
#include "GaussLegendreGrid.h"
#include "Policies.h"

namespace GSHTrans {

namespace TuningDetails {

using Clock = std::chrono::steady_clock;

// Best of several windows, not the mean of one.
//
// The development machine's clock scales under load: repeating one
// measurement seconds later moves it by tens of per cent, which is enough to
// invent a difference that is not there. The benchmark harness has said so
// since step E and this is the same discipline made a constant of the
// library rather than a paragraph in a plan.
//
// The action is run once before timing starts, because the first call plans,
// faults its pages and fills its caches, and none of that is what is being
// compared.
template <typename Action>
double BestSeconds(Action&& action, int windows = 3) {
  action();
  auto best = std::numeric_limits<double>::max();
  for (auto window = 0; window < windows; ++window) {
    const auto start = Clock::now();
    action();
    const auto elapsed =
        std::chrono::duration<double>(Clock::now() - start).count();
    best = std::min(best, elapsed);
  }
  return best;
}

// One forward and one inverse over the caller's shape, on scratch of the
// right size. What is being measured is a schedule; a schedule does not
// depend on the values.
template <typename GridType, typename Complex>
double TimeRound(const GridType& grid, std::ptrdiff_t lMax, std::ptrdiff_t n,
                 std::ptrdiff_t count, Execution policy, int windows = 3) {
  using Int = std::ptrdiff_t;
  const auto fieldSize = static_cast<Int>(grid.FieldSize());
  const auto coefficientSize = static_cast<Int>(grid.CoefficientSize(lMax, n));

  auto fields = FFTWpp::vector<Complex>(
      static_cast<std::size_t>(count * fieldSize), Complex{1, 0});
  auto coefficients = FFTWpp::vector<Complex>(
      static_cast<std::size_t>(count * coefficientSize));
  const auto fieldBatch = Batch::Contiguous(count, fieldSize);
  const auto coefficientBatch = Batch::Contiguous(count, coefficientSize);

  return BestSeconds(
      [&] {
        grid.ForwardTransformation(lMax, n, fields, fieldBatch, coefficients,
                                   coefficientBatch, policy);
        grid.InverseTransformation(lMax, n, coefficients, coefficientBatch,
                                   fields, fieldBatch, policy);
      },
      windows);
}

}  // namespace TuningDetails

//--------------------------------------------------------------------------//
//                              What Tune returns                            //
//--------------------------------------------------------------------------//

// The chunking policy a measurement chose, and enough of the measurement to
// see why.
//
// `conclusive` is [C22]. The noise floor on the development machine is about
// ten per cent and several of the differences at stake are smaller, so a
// candidate must beat the incumbent by more than a margin to displace it, and
// the default is the incumbent. An entry recording an inconclusive comparison
// is one a later, quieter run may usefully revisit; one recording a three per
// cent win looks like knowledge and is not.
struct TunedChunking {
  Chunking chunking = Chunking::Automatic();
  bool conclusive = false;
  double seconds = 0;         // the chosen candidate
  double defaultSeconds = 0;  // the incumbent, for comparison
  int candidates = 0;         // distinct schedules there were to choose from

  // How much better the choice is than doing nothing. One, or a little under,
  // when the comparison was inconclusive.
  double Speedup() const {
    return seconds > 0 ? defaultSeconds / seconds : 1;
  }
};

//--------------------------------------------------------------------------//
//                             Tuning the chunk                              //
//--------------------------------------------------------------------------//

// The margin a candidate must beat the incumbent by, per [C22]. Ten per cent
// is the measured noise floor of the development machine, and figures below
// it have repeatedly meant nothing here.
inline constexpr double TuningMargin = 0.10;

// The cache figures swept, in bytes. Powers of two from 1 MiB to 64 MiB,
// which spans a laptop's shared L3 and a server CCD's, and is small enough
// that the whole sweep costs a fraction of building the table it runs beside.
//
// **The cache figure is what is tuned, not the chunk** ([C18], and section
// 12.4). Chunking::Count already takes the per-field size and the number of
// live copies and computes a chunk from them, and it already tells the two
// directions apart through `copies`. So fitting the one number the formula
// does not know means one answer serves every batch size and both
// directions, where a chunk fitted for one batch would not.
inline std::vector<std::ptrdiff_t> TuningCacheCandidates() {
  auto candidates = std::vector<std::ptrdiff_t>{};
  for (auto megabytes = std::ptrdiff_t{1}; megabytes <= 64; megabytes *= 2) {
    candidates.push_back(megabytes << 20);
  }
  return candidates;
}

// Time the caller's own problem under each candidate and return the best.
//
// Both directions are timed and their times added, because a grid carries one
// chunking policy and both directions read it. The transforms run on scratch
// of the right shape rather than on the caller's data: what is being measured
// is a schedule, and a schedule does not depend on the values.
//
// Cheap because of [C18]: every candidate is `grid.With(...)`, a pointer copy
// sharing one table, where before this would have built a table per
// candidate -- 648 MB and 0.16 s apiece at lMax = 256, to choose an integer.
//
// Complex fields, deliberately, since that is the general case; a caller
// whose work is entirely real scalars has half the coefficient block and may
// find a different optimum, which is what the `count` and `n` arguments are
// for -- tune the shape you will run.
//
// **Two things about the method, both of which the first version got wrong
// and the numbers caught.**
//
// *The candidates are interleaved, not run one after another.* core-plan.md
// section 8 is explicit that comparing two versions requires running them
// alternately in one session, because this machine's clock scales under load;
// timing the incumbent to completion and then each candidate to completion is
// exactly what that rule forbids. Measured, the first pass over a fixed grid
// runs about nine per cent slow against the sixth -- just under the margin,
// and systematically against whichever candidate goes first, which is the
// incumbent. So every candidate is warmed up before any is timed, and the
// windows go round-robin.
//
// *Candidates that would run the same schedule are collapsed.* Two cache
// figures that yield the same chunk are the same experiment, and timing both
// only gives noise two chances to beat the margin. The clearest case is
// `count == 1`, where every chunk of one or more processes the whole batch in
// one go, so **there is nothing to tune at all** -- and the first version
// duly reported a 2x "win" there, which is definitionally impossible.
//
// The chunk a candidate would use is worked out here rather than asked of the
// grid, mirroring ForwardChunkSize and InverseChunkSize. If that mirror ever
// drifts, the cost is a few redundant measurements and not a wrong answer:
// deduplication is a guard against spurious wins, not a correctness
// requirement.
template <typename GridType>
TunedChunking TuneChunking(const GridType& grid, std::ptrdiff_t lMax,
                           std::ptrdiff_t n, std::ptrdiff_t count,
                           Execution policy = Execution::Sequential(),
                           int windows = 3) {
  using Int = std::ptrdiff_t;
  using Complex = typename GridType::Complex;

  if (count < 1) {
    throw std::invalid_argument("Tuning: the batch count must be positive");
  }

  const auto fieldSize = static_cast<Int>(grid.FieldSize());
  const auto coefficientSize = static_cast<Int>(grid.CoefficientSize(lMax, n));

  auto fields = FFTWpp::vector<Complex>(
      static_cast<std::size_t>(count * fieldSize), Complex{1, 0});
  auto coefficients = FFTWpp::vector<Complex>(
      static_cast<std::size_t>(count * coefficientSize));

  const auto fieldBatch = Batch::Contiguous(count, fieldSize);
  const auto coefficientBatch = Batch::Contiguous(count, coefficientSize);

  const auto Round = [&](const GridType& g) {
    g.ForwardTransformation(lMax, n, fields, fieldBatch, coefficients,
                            coefficientBatch, policy);
    g.InverseTransformation(lMax, n, coefficients, coefficientBatch, fields,
                            fieldBatch, policy);
  };

  // The schedule a policy would actually run, as a pair: the forward's
  // copies are its threads and the inverse's are one, and a chunk beyond the
  // batch is the batch.
  const auto blockBytes = coefficientSize * static_cast<Int>(sizeof(Complex));
  const auto threads =
      policy.IsParallel()
          ? (policy.Threads() > 0 ? policy.Threads() : omp_get_max_threads())
          : 1;
  const auto Schedule = [&](const Chunking& chunking) {
    return std::pair(std::min(chunking.Count(blockBytes, threads), count),
                     std::min(chunking.Count(blockBytes, 1), count));
  };

  struct Candidate {
    Chunking chunking;
    double best = std::numeric_limits<double>::max();
  };

  auto seen = std::vector<std::pair<Int, Int>>{};
  auto candidates = std::vector<Candidate>{};
  const auto Consider = [&](Chunking chunking) {
    const auto schedule = Schedule(chunking);
    if (std::find(seen.begin(), seen.end(), schedule) != seen.end()) return;
    seen.push_back(schedule);
    candidates.push_back(Candidate{chunking});
  };

  // The incumbent goes first, so that it is candidates[0] and holds any tie.
  Consider(grid.ChunkingPolicy());
  for (auto bytes : TuningCacheCandidates()) {
    Consider(Chunking::ForCache(bytes));
  }

  auto tuned = std::vector<GridType>{};
  tuned.reserve(candidates.size());
  for (const auto& candidate : candidates) {
    tuned.push_back(grid.With(candidate.chunking));
  }

  // Warm every candidate before timing any: the first pass over a shape
  // plans, faults its pages and fills its caches, and paying that once per
  // candidate would charge it to whichever ran first.
  for (const auto& g : tuned) Round(g);

  for (auto window = 0; window < windows; ++window) {
    for (auto i = std::size_t{0}; i < candidates.size(); ++i) {
      const auto start = TuningDetails::Clock::now();
      Round(tuned[i]);
      const auto elapsed =
          std::chrono::duration<double>(TuningDetails::Clock::now() - start)
              .count();
      candidates[i].best = std::min(candidates[i].best, elapsed);
    }
  }

  auto result = TunedChunking{};
  result.chunking = candidates.front().chunking;
  result.defaultSeconds = candidates.front().best;
  result.seconds = result.defaultSeconds;

  // [C22]: the incumbent holds unless a candidate beats it by the margin, and
  // the result records whether one did. A tuner that picks the nominal winner
  // of a seven per cent difference is picking noise, and will pick
  // differently next run.
  for (auto i = std::size_t{1}; i < candidates.size(); ++i) {
    if (candidates[i].best < result.seconds * (1 - TuningMargin)) {
      result.chunking = candidates[i].chunking;
      result.seconds = candidates[i].best;
      result.conclusive = true;
    }
  }

  // How many distinct schedules there were to choose between. One means the
  // question did not arise -- at count == 1 it never does.
  result.candidates = static_cast<int>(candidates.size());

  return result;
}

//--------------------------------------------------------------------------//
//                            Tuning the kernel                              //
//--------------------------------------------------------------------------//

// Which of the two Legendre kernels a machine should use for a given problem.
//
// [C12] keeps both permanently, and this is the customer that decision was
// taken for: unlike every other knob here, the kernel has **two complete
// implementations that compute the same answer**, so timing both on the
// caller's actual problem is a well-posed measurement rather than a
// heuristic. thoughts.md section 10 calls it the strongest case for the whole
// mechanism, and the measured differences -- 2x to 6x, against the chunk's
// 1.1x to 1.2x -- say so.
//
// Grid *parameters* rather than a grid, because a grid already has a kernel
// baked into its table and the point is to build one of each.
struct TunedKernel {
  TransformKernel kernel = TransformKernel::Loop();
  bool conclusive = false;
  bool matrixTried = false;
  double loopSeconds = 0;
  double matrixSeconds = 0;

  // Why the matrix kernel was not measured, empty when it was. [C19]: a tuner
  // is the machinery most likely to substitute silently, so every reason it
  // did not do what was asked is a string the caller can print.
  std::string skipped;

  double Speedup() const {
    if (!matrixTried || matrixSeconds <= 0) return 1;
    return loopSeconds / matrixSeconds;
  }
};

// The matrix kernel was not available, so only the loop is timed. Its number
// is still reported, because a caller comparing machines wants it.
template <typename GridType>
TunedKernel TuneKernelLoopOnly(TunedKernel result, std::ptrdiff_t lMax,
                               std::ptrdiff_t nMax, std::ptrdiff_t n,
                               std::ptrdiff_t count, Execution policy,
                               FFTWpp::Flag flag, Chunking chunking,
                               WignerValues values, int rounds) {
  using Complex = typename GridType::Complex;
  const auto grid = GridType(lMax, nMax, flag, chunking, values,
                             TransformKernel::Loop());
  auto best = std::numeric_limits<double>::max();
  for (auto round = 0; round < rounds; ++round) {
    best = std::min(best, TuningDetails::TimeRound<GridType, Complex>(
                              grid, lMax, n, count, policy));
  }
  result.loopSeconds = best;
  return result;
}

// Build one grid of each kernel in turn, time the caller's problem on each,
// and return the winner.
//
// **Sequential construction, per [C20], and the reason is memory.** section
// 11.4's M5 records that comparing kernels means two grids and that at
// lMax = 256 both live costs 1.3 GB. Each is therefore built, timed and
// destroyed before the next, so the peak is one table.
//
// **But the order alternates between rounds, which [C20] did not ask for and
// W3's numbers say it needs.** Sequential comparison has an ordering bias:
// the first pass over a fixed grid measures about nine per cent slow against
// the sixth on this machine, whichever code it is running, so whatever is
// measured second wins a little for free. Running the pair twice with the
// order reversed makes the bias symmetric, at the price of four table builds
// rather than [C20]'s three -- about 0.9 s at lMax = 256, against a decision
// worth several times the transform.
//
// The bias is well under the differences at stake here, so this is insurance
// rather than a correction. It matters for a marginal result, and a marginal
// result is exactly the one [C22]'s margin is there to refuse.
template <typename GridType>
TunedKernel TuneKernel(std::ptrdiff_t lMax, std::ptrdiff_t nMax,
                       std::ptrdiff_t n, std::ptrdiff_t count,
                       Execution policy = Execution::Sequential(),
                       FFTWpp::Flag flag = FFTWpp::Measure,
                       Chunking chunking = Chunking::Automatic(),
                       WignerValues values = WignerValues::Stored(),
                       int rounds = 2) {
  using Int = std::ptrdiff_t;
  using Real = typename GridType::Real;
  using Complex = typename GridType::Complex;

  if (count < 1) {
    throw std::invalid_argument("Tuning: the batch count must be positive");
  }
  if (rounds < 1) {
    throw std::invalid_argument("Tuning: rounds must be at least one");
  }

  auto result = TunedKernel{};

  // The three ways the matrix kernel can be unavailable, each named rather
  // than silently collapsing to "use the loop" ([C17], [C19]).
#ifndef GSHTRANS_HAVE_BLAS
  result.skipped = "this build has no BLAS, so the matrix kernel does not "
                   "exist";
  return TuneKernelLoopOnly<GridType>(std::move(result), lMax, nMax, n, count,
                                      policy, flag, chunking, values, rounds);
#else
  if constexpr (!BlasDetails::BlasReal<Real>) {
    result.skipped = "BLAS offers single and double precision only, so the "
                     "matrix kernel is unavailable at this precision";
    return TuneKernelLoopOnly<GridType>(std::move(result), lMax, nMax, n,
                                        count, policy, flag, chunking, values,
                                        rounds);
  } else {
    if (!values.AreStored()) {
      result.skipped =
          "the matrix kernel needs a stored table, and generated values were "
          "asked for";
      return TuneKernelLoopOnly<GridType>(std::move(result), lMax, nMax, n,
                                          count, policy, flag, chunking,
                                          values, rounds);
    }

    result.matrixTried = true;
    auto loopBest = std::numeric_limits<double>::max();
    auto matrixBest = std::numeric_limits<double>::max();

    for (auto round = 0; round < rounds; ++round) {
      const auto loopFirst = (round % 2) == 0;
      for (auto step = 0; step < 2; ++step) {
        const auto wantLoop = (step == 0) == loopFirst;
        const auto kernel = wantLoop ? TransformKernel::Loop()
                                     : TransformKernel::Matrix();
        const auto grid =
            GridType(lMax, nMax, flag, chunking, values, kernel);
        const auto seconds = TuningDetails::TimeRound<GridType, Complex>(
            grid, lMax, n, count, policy);
        auto& best = wantLoop ? loopBest : matrixBest;
        best = std::min(best, seconds);
      }
    }

    result.loopSeconds = loopBest;
    result.matrixSeconds = matrixBest;
    if (matrixBest < loopBest * (1 - TuningMargin)) {
      result.kernel = TransformKernel::Matrix();
      result.conclusive = true;
    }
    return result;
  }
#endif
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_TUNING_GUARD_H
