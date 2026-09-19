#pragma once

/**
 * @file Utility.hpp
 * @brief Small numerical helpers used across the library.
 */

#include <atomic>
#include <cstddef>
#include <exception>
#include <initializer_list>
#include <mutex>
#include <utility>

namespace GSHTrans {

/**
 * @brief Returns @f$(-1)^m@f$.
 * @tparam Scalar The type of the result.
 * @param m The exponent; only its parity is used.
 */
template <typename Scalar = std::ptrdiff_t>
constexpr auto MinusOneToPower(std::ptrdiff_t m) -> Scalar {
  return m % 2 ? -1 : 1;
}

/**
 * @brief Whether FFTW has hard-coded codelets for a transform of length @p n.
 *
 * @details The lengths it does are @f$n = 2^a 3^b 5^c 7^d 11^e 13^f@f$ with
 * @f$e + f@f$ either zero or one. Any other length falls back to a general
 * algorithm, which for a large prime factor is much slower — and the grid's
 * longitude count runs into exactly that, since the smallest length resolving
 * orders @f$|m| \le 256@f$ is 513, whose next even neighbour 514 is
 * @f$2 \times 257@f$ with 257 prime.
 *
 * @param n The length to test.
 * @return True if @p n has only the small prime factors above.
 */
constexpr auto IsFastFFTSize(std::ptrdiff_t n) -> bool {
  if (n < 1) return false;
  for (auto p : {2, 3, 5, 7}) {
    while (n % p == 0) n /= p;
  }
  auto large = 0;
  for (auto p : {11, 13}) {
    if (n % p == 0) {
      n /= p;
      large++;
    }
  }
  return n == 1 && large <= 1;
}

/**
 * @brief The smallest length at least @p n for which IsFastFFTSize holds.
 * @param n The lower bound.
 */
constexpr auto FastFFTSize(std::ptrdiff_t n) -> std::ptrdiff_t {
  if (n < 1) return 1;
  while (!IsFastFFTSize(n)) n++;
  return n;
}

namespace Details {

/**
 * @brief Carries an exception out of an OpenMP region.
 *
 * @details An exception that leaves a structured block inside a parallel
 * region does not propagate: the program is terminated. So nothing may leave
 * one, and everything that can throw inside a region -- a per-thread
 * workspace being allocated, a plan being made, and above all a *caller's*
 * radial operator, which the library knows nothing about -- runs through
 * Run(). The first exception is kept and the rest are dropped, which is what
 * a sequential loop would have shown the caller anyway; Rethrow() after the
 * region hands it on, so a call that throws when it runs sequentially throws
 * the same thing when it runs threaded.
 *
 * Two rules for using it, both of which come from OpenMP and not from here.
 *
 * **Every thread must still reach every worksharing construct.** A thread
 * whose set-up failed cannot skip the `omp for` that follows, or the others
 * wait at its barrier for ever. It goes through the loop and does nothing,
 * which is what Failed() is for.
 *
 * **Failed() is a hint and not a lock.** It is read without synchronisation
 * beyond its own atomicity, so an iteration may start after another has
 * failed. That is harmless -- each iteration is guarded by its own Run() --
 * and it is what keeps the cost to one load per iteration of an *outer*
 * loop. It is not for inner loops.
 */
class ExceptionCapture {
 public:
  /**
   * @brief Calls @p f, keeping whatever it throws.
   * @param f Called with no arguments.
   */
  template <typename F>
  void Run(F&& f) noexcept {
    try {
      std::forward<F>(f)();
    } catch (...) {
      Keep(std::current_exception());
    }
  }

  /** @brief Whether anything has been kept, so that remaining work can be
   * skipped. */
  bool Failed() const noexcept {
    return failed_.load(std::memory_order_relaxed);
  }

  /** @brief Throws what was kept, if anything was. Call after the region. */
  void Rethrow() const {
    if (exception_) std::rethrow_exception(exception_);
  }

 private:
  void Keep(std::exception_ptr exception) noexcept {
    {
      const auto lock = std::lock_guard(mutex_);
      if (!exception_) exception_ = std::move(exception);
    }
    failed_.store(true, std::memory_order_relaxed);
  }

  std::mutex mutex_;
  std::exception_ptr exception_;
  std::atomic<bool> failed_{false};
};

}  // namespace Details

}  // namespace GSHTrans
