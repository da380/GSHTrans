#pragma once

/**
 * @file OpenMP.hpp
 * @brief OpenMP where there is one, and nothing where there is not.
 *
 * @details The library threads with OpenMP and does not need it. Built
 * without `-fopenmp` every parallel region is a serial loop, every policy
 * that asks for threads gets one, and the answers are the same -- so copying
 * the `GSHTrans` directory into a project needs a compiler and FFTW and
 * nothing else, and the library builds where OpenMP is awkward to come by.
 * It is also the build to take to a debugger or a race detector.
 *
 * Whether it is there is the compiler's to say, through `_OPENMP`, which is
 * defined exactly when `-fopenmp` or its equivalent is given. There is no
 * macro of this library's to set, and so none to set inconsistently between
 * two translation units.
 *
 * Two things are needed to make that work without noise, and they are all
 * this file is.
 *
 * **The runtime calls go through wrappers.** Not functions named `omp_*`
 * defined here for the serial case: GCC ships `<omp.h>` whether or not
 * OpenMP is enabled, and a caller who includes it would then have two
 * declarations of every one.
 *
 * **The pragmas are left as pragmas**, and the warning about them is
 * silenced where it would otherwise be raised. A `#pragma omp` the compiler
 * was not asked to understand is ignored, which is right, and warned about
 * under `-Wall`, which in a header-only library means warned about in the
 * caller's build. Each header that has one turns `-Wunknown-pragmas` off for
 * its own length, and only when `_OPENMP` is not defined, so a caller's own
 * stray pragma is still reported. (Wrapping them in a macro was tried first.
 * It works, and no formatter can lay out a macro whose argument says `for`
 * and `if`.)
 */

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GSHTrans::OpenMP {

/** @brief Whether this translation unit was compiled with OpenMP. */
#ifdef _OPENMP
inline constexpr bool Available = true;
#else
inline constexpr bool Available = false;
#endif

/** @brief `omp_get_max_threads()`, or one. */
inline int MaxThreads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

/** @brief `omp_get_num_threads()`: the size of the current team, or one. */
inline int TeamThreads() {
#ifdef _OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}

/** @brief `omp_get_thread_num()`: this thread's number in its team, or zero. */
inline int ThreadNumber() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

/** @brief `omp_in_parallel()`: whether an active region encloses this. */
inline bool InParallel() {
#ifdef _OPENMP
  return omp_in_parallel() != 0;
#else
  return false;
#endif
}

/**
 * @brief `omp_set_num_threads()` for the calling task, or nothing.
 * @param threads How many threads a region opened from here may have.
 */
inline void SetThreads([[maybe_unused]] int threads) {
#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif
}

}  // namespace GSHTrans::OpenMP
