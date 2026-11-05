#pragma once

/**
 * @file Parallel.hpp
 * @brief Parallel loops and a deterministic reduction, under the library's
 * threading rules.
 *
 * @details Four function templates for the loops a caller would otherwise
 * write as an OpenMP region of their own: pointwise work over a stack of grid
 * values, per-slice or per-degree tasks, and sums over large arrays. Each
 * takes an Execution policy and keeps the rules every threaded operation in
 * the library keeps:
 *
 * - **One level threads.** The team is `Execution::TeamSize()`: one thread
 *   when the policy is sequential or the call is already inside an active
 *   parallel region. A body that calls a threaded GSHTrans transform
 *   therefore gets a sequential transform inside a parallel loop; transforms
 *   belong in the batched calls, outside such loops, and a loop over slices
 *   that transforms one field per iteration gives up the batching that makes
 *   a transform cheap.
 * - **Exceptions leave as they would sequentially.** A body may throw. The
 *   exception rethrown after the loop is the one a sequential loop would
 *   have thrown: the one from the lowest failing index. After a failure at
 *   index `j`, indices above `j` may be skipped but every index below it
 *   still runs, so the lowest failure is always found, under a static or a
 *   dynamic schedule.
 * - **Results do not depend on the thread count.** The loops write only what
 *   each index owns; Reduce sums over blocks fixed independently of the
 *   thread count and combines them in block order, so it gives the same bits
 *   on any team, including one thread.
 * - **Nothing without OpenMP.** Built without it, every form runs the same
 *   blocks in order on the calling thread.
 *
 * They take plain index ranges, so they work on the caller's own storage and
 * layout; nothing here knows about grids or fields. Distribution over
 * processes is the caller's: apply a cross-process combine (an MPI
 * all-reduce, say) to what Reduce returns.
 */

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <limits>
#include <mutex>
#include <utility>
#include <vector>

#include "OpenMP.hpp"
#include "Policies.hpp"

// The pragmas below are for a compiler that was asked for OpenMP. One that was
// not ignores them, which is right, and warns about them under -Wall, which in
// a header is a warning in the caller's build; so the warning is off for this
// file's length, and only then. See OpenMP.hpp.
#ifndef _OPENMP
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

namespace GSHTrans {

namespace Details {

/**
 * @brief Keeps the exception from the lowest failing index of a loop.
 *
 * @details Skip(i) says whether index @p i may be skipped: true once a
 * failure at a lower index is known. Keep(i, e) records @p e if @p i is lower
 * than any failure recorded so far. Every index below the lowest failure is
 * therefore run, and the exception rethrown is the one a sequential loop
 * would have stopped at.
 */
class IndexedCapture {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.

  /** @brief Whether a failure below @p i has been recorded. */
  bool Skip(Int i) const noexcept {
    return i > lowest_.load(std::memory_order_relaxed);
  }

  /** @brief Records @p exception if @p i is the lowest failure so far. */
  void Keep(Int i, std::exception_ptr exception) noexcept {
    const auto lock = std::lock_guard(mutex_);
    if (!exception_ || i < index_) {
      exception_ = std::move(exception);
      index_ = i;
      lowest_.store(i, std::memory_order_relaxed);
    }
  }

  /** @brief Throws what was kept, if anything was. Call after the loop. */
  void Rethrow() const {
    if (exception_) std::rethrow_exception(exception_);
  }

 private:
  std::mutex mutex_;
  std::exception_ptr exception_;
  Int index_ = 0;
  std::atomic<Int> lowest_{std::numeric_limits<Int>::max()};
};

// The contiguous block of [0, n) that member `t` of a team of `threads`
// takes. Shared by every static form, so two loops over the same range hand
// each thread the same block.
inline std::pair<std::ptrdiff_t, std::ptrdiff_t> StaticBlock(
    std::ptrdiff_t n, std::ptrdiff_t t, std::ptrdiff_t threads) {
  return {n * t / threads, n * (t + 1) / threads};
}

// Run `block(begin, end, capture)` once per member of the team over its static
// block of [0, n), then rethrow the lowest failure. The block is responsible
// for catching what its body throws and handing it to `capture` with the
// index it failed at; doing that once per block rather than once per element
// leaves an element-wise body's inner loop free to vectorise. Sequentially the
// same function runs once over the whole range, so both paths keep the same
// exception.
template <typename Block>
void StaticRegion(Execution policy, std::ptrdiff_t n, Block& block) {
  if (n <= 0) return;
  auto capture = IndexedCapture{};
  const auto threads = policy.TeamSize();
  if (threads == 1) {
    block(std::ptrdiff_t{0}, n, capture);
  } else {
#pragma omp parallel num_threads(threads)
    {
      const auto [begin, end] =
          StaticBlock(n, OpenMP::ThreadNumber(), OpenMP::TeamThreads());
      if (begin < end && !capture.Skip(begin)) block(begin, end, capture);
    }
  }
  capture.Rethrow();
}

}  // namespace Details

/**
 * @brief Calls `body(begin, end)` on contiguous blocks covering `[0, n)`.
 *
 * @details One block per member of the team, assigned statically: two calls
 * over the same `n` with the same team hand each thread the same block, which
 * keeps a pass over memory a previous pass touched on the thread that touched
 * it. Sequentially the body is called once, with `[0, n)`. This is the form
 * for element-wise work over large arrays: the body's own inner loop is a
 * plain loop the compiler can vectorise.
 *
 * The body must write only what its indices own. If it throws, it should do
 * so at the first index in its range that fails, processing its range in
 * order, so that the exception rethrown is the sequential one.
 *
 * @param policy Whether to thread, and how many threads.
 * @param n The length of the range; nothing is called when it is not
 * positive.
 * @param body Called as `body(std::ptrdiff_t begin, std::ptrdiff_t end)`.
 */
template <typename Body>
void ForEachBlock(Execution policy, std::ptrdiff_t n, Body&& body) {
  auto block = [&](std::ptrdiff_t begin, std::ptrdiff_t end,
                   Details::IndexedCapture& capture) {
    try {
      body(begin, end);
    } catch (...) {
      capture.Keep(begin, std::current_exception());
    }
  };
  Details::StaticRegion(policy, n, block);
}

/**
 * @brief Calls `body(i)` for every `i` in `[0, n)`, in static blocks.
 *
 * @details ForEachBlock with the per-index loop written for the caller. The
 * index of a throwing call is known exactly, so the exception rethrown is the
 * lowest failing index's whatever the body does.
 *
 * @param policy Whether to thread, and how many threads.
 * @param n The length of the range.
 * @param body Called as `body(std::ptrdiff_t i)`; writes only what @p i owns.
 */
template <typename Body>
void ForEach(Execution policy, std::ptrdiff_t n, Body&& body) {
  auto block = [&](std::ptrdiff_t begin, std::ptrdiff_t end,
                   Details::IndexedCapture& capture) {
    auto i = begin;
    try {
      for (; i < end; ++i) body(i);
    } catch (...) {
      capture.Keep(i, std::current_exception());
    }
  };
  Details::StaticRegion(policy, n, block);
}

/**
 * @brief Calls `body(i)` for every `i` in `[0, n)`, handing indices out
 * dynamically.
 *
 * @details For tasks of uneven cost: per-degree solves, or work sorted
 * largest first so that the last tasks handed out are small. Which thread
 * runs which index varies from run to run, so the body must write only what
 * its index owns.
 *
 * @param policy Whether to thread, and how many threads.
 * @param n The number of tasks.
 * @param body Called as `body(std::ptrdiff_t i)`.
 */
template <typename Body>
void ForEachTask(Execution policy, std::ptrdiff_t n, Body&& body) {
  if (n <= 0) return;
  const auto threads = policy.TeamSize();
  if (threads == 1) {
    for (auto i = std::ptrdiff_t{0}; i < n; ++i) body(i);
    return;
  }
  auto capture = Details::IndexedCapture{};
#pragma omp parallel for schedule(dynamic) num_threads(threads)
  for (std::ptrdiff_t i = 0; i < n; ++i) {
    if (capture.Skip(i)) continue;
    try {
      body(i);
    } catch (...) {
      capture.Keep(i, std::current_exception());
    }
  }
  capture.Rethrow();
}

/**
 * @brief The most blocks Reduce divides its range into.
 * @details Fixed, and independent of the thread count, which is what makes the
 * result independent of it. Enough blocks to keep 64 threads busy several
 * times over; few enough that combining them is negligible.
 */
inline constexpr std::ptrdiff_t ReduceBlocks = 256;

/**
 * @brief `combine` over `body(i)` for `i` in `[0, n)`, with the same bits on
 * any team.
 *
 * @details The range is divided into `min(n, ReduceBlocks)` contiguous blocks
 * fixed by `n` alone. Each block's partial is formed in index order as
 * `p = body(begin); p = combine(p, body(i))` for the rest; the result is
 * `combine(...combine(combine(init, p_0), p_1)..., p_last)`. That grouping is
 * the same whoever computes which block, so the result is bitwise
 * independent of the thread count, sequential included -- provided the
 * build does not license the compiler to reassociate floating-point
 * arithmetic (`-ffast-math` and the like).
 *
 * @param policy Whether to thread, and how many threads.
 * @param n The length of the range; @p init is returned when it is not
 * positive.
 * @param init The starting value, combined once, first.
 * @param body Called as `body(std::ptrdiff_t i)`, returning a `T`.
 * @param combine Called as `combine(T, T)`, returning a `T`.
 * @return The reduction.
 */
template <typename T, typename Body, typename Combine>
T Reduce(Execution policy, std::ptrdiff_t n, T init, Body&& body,
         Combine&& combine) {
  if (n <= 0) return init;
  const auto blocks = std::min(n, ReduceBlocks);
  auto partials = std::vector<T>(static_cast<std::size_t>(blocks), init);
  const auto one = [&](std::ptrdiff_t b) {
    const auto [begin, end] = Details::StaticBlock(n, b, blocks);
    T p = body(begin);
    for (auto i = begin + 1; i < end; ++i) p = combine(p, body(i));
    partials[static_cast<std::size_t>(b)] = p;
  };
  ForEach(policy, blocks, one);
  auto result = init;
  for (const auto& p : partials) result = combine(result, p);
  return result;
}

/** @brief Reduce with `std::plus`: a deterministic sum. */
template <typename T, typename Body>
T Sum(Execution policy, std::ptrdiff_t n, T init, Body&& body) {
  return Reduce(policy, n, init, std::forward<Body>(body),
                [](const T& a, const T& b) { return a + b; });
}

}  // namespace GSHTrans

#ifndef _OPENMP
#pragma GCC diagnostic pop
#endif
