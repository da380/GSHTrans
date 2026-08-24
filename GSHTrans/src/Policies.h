#ifndef GSH_TRANS_POLICIES_GUARD_H
#define GSH_TRANS_POLICIES_GUARD_H

/**
 * @file Policies.h
 * @brief How an operation may thread, how a caller's batch is laid out, how
 * much of one the inner loop takes at a time, whether a grid stores its Wigner
 * values, which Legendre kernel it uses, and which interpolation scheme a
 * field offers.
 *
 * @details What these have in common is that each is a decision the caller
 * makes about *the machine* rather than about the mathematics, which is why
 * they live together and apart from the numeric concepts.
 *
 * They are values rather than macros or template parameters, deliberately.
 * The library is header-only, so a preprocessor knob defined differently in
 * two translation units would give a grid two inline bodies and let the linker
 * pick one in silence; and a template parameter would put a machine property
 * into a type that the mathematics has no business carrying. Scheme is the one
 * exception, and says below why.
 */

#include <algorithm>
#include <cstddef>
#include <stdexcept>

namespace GSHTrans {

//-------------------------------------------------------------------------//
//                            Execution policy                              //
//--------------------------------------------------------------------------//

/**
 * @brief Whether an operation may use threads, and how many.
 *
 * @details Sequential by default everywhere: a library should not create
 * threads because it can, only because it was asked to. The rule this exists
 * to make keepable is that exactly one level threads — a caller parallelising
 * over slices, components or realisations calls the transform sequentially,
 * while a caller with one large problem asks the transform to thread. Both at
 * once is worse than either, so an operation asked to run in parallel from
 * inside an existing parallel region runs sequentially instead.
 */
class Execution {
 public:
  /** @brief One thread, which is the default everywhere. */
  static Execution Sequential() { return Execution(1); }

  /**
   * @brief Threads permitted.
   * @param threads How many; zero means whatever OpenMP would choose, which
   * respects `OMP_NUM_THREADS`. On a machine with simultaneous multithreading
   * that is usually the number of hardware threads, and for this library's
   * memory-bound work that measures *slower* than one thread per core, so a
   * caller who cares should say what they want.
   */
  static Execution Parallel(int threads = 0) {
    return Execution(threads > 0 ? threads : 0);
  }

  /** @brief The thread count, zero meaning OpenMP's choice. */
  auto Threads() const { return _threads; }

  /** @brief Whether threading is permitted at all. */
  auto IsParallel() const { return _threads != 1; }

  /** @brief Compares componentwise. */
  bool operator==(const Execution&) const = default;

 private:
  explicit Execution(int threads) : _threads{threads} {}
  int _threads;
};

//-------------------------------------------------------------------------//
//                             Batch descriptor                             //
//--------------------------------------------------------------------------//

/**
 * @brief How a batch of same-spin fields is laid out in the caller's storage.
 *
 * @details The batched transform's primitive is "transform @f$k@f$ fields
 * sharing a grid, a degree and an upper index", with the single field as
 * @f$k = 1@f$. What varies between callers is not how many fields there are
 * but how they are arranged, so a batch is described by
 * `(count, stride, dist)` rather than by contiguity: element @f$j@f$ of field
 * @f$k@f$ lives at `j * stride + k * dist`.
 *
 * This is FFTW's advanced-interface form, and it is used because it costs
 * nothing at the stage that actually reads the caller's layout — the FFT —
 * while covering both arrangements the field layer produces: a radial slab or
 * component-major tensor components, which are `stride = 1`,
 * `dist = FieldSize`; and point-major tensor components, which are
 * `stride = nComponents`, `dist = 1`. The second is why this is a decision
 * rather than a detail: without it the field layer would have to promise
 * component-major storage and repack anything else before transforming.
 *
 * Full FFTW `nembed` generality is deliberately not offered, since the
 * angular layout is fixed by the grid and most of it would be unreachable.
 *
 * A batch shares grid, degree and upper index. The Wigner block for that
 * upper index is precisely what batching amortises, so fields at different
 * upper indices cannot batch together — for a rank-2 tensor that means
 * batching over radii within each @f$n@f$, not across components.
 */
class Batch {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.

  /**
   * @brief Fields laid end to end, each contiguous.
   * @param count How many fields.
   * @param size How many elements each holds.
   */
  static Batch Contiguous(Int count, Int size) {
    return Batch(count, 1, size);
  }

  /**
   * @brief Fields interleaved element by element.
   * @param count How many fields take part; this may be smaller than @p
   * stride, so a batch of two components out of a point-major five is
   * `Interleaved(2, 5)` and the other three are not part of the call.
   * @param stride How many fields lie between successive elements of one.
   */
  static Batch Interleaved(Int count, Int stride) {
    return Batch(count, stride, 1);
  }

  /**
   * @brief The general form, for a layout that is neither.
   * @param count How many fields.
   * @param stride The separation of successive elements of one field.
   * @param dist The separation of successive fields.
   */
  static Batch Strided(Int count, Int stride, Int dist) {
    return Batch(count, stride, dist);
  }

  /**
   * @brief The single field, which is what the unbatched entry points hand
   * down.
   * @param size How many elements it holds.
   */
  static Batch One(Int size) { return Batch(1, 1, size); }

  /** @brief How many fields take part. */
  auto Count() const { return _count; }
  /** @brief The separation of successive elements of one field. */
  auto Stride() const { return _stride; }
  /** @brief The separation of successive fields. */
  auto Dist() const { return _dist; }

  /**
   * @brief Where element @p j of field @p k lives.
   * @param j The element index within a field.
   * @param k The field index within the batch.
   */
  Int Offset(Int j, Int k) const { return j * _stride + k * _dist; }

  /**
   * @brief The smallest range that holds this batch.
   *
   * @details The batched entry points check the caller's range against this
   * rather than for equality: an interleaved batch is a window onto a larger
   * range whose other elements are none of the transform's business, so
   * demanding an exact size would reject exactly the layout this type exists
   * to accept. The unbatched entry points keep the equality check, since their
   * contract is that the range *is* the field.
   *
   * @param size How many elements each field holds.
   */
  Int Span(Int size) const {
    return (size - 1) * _stride + (_count - 1) * _dist + 1;
  }

  /**
   * @brief Whether the fields of this batch are disjoint at that size.
   *
   * @details They must be: the transform writes every element of every field,
   * so overlapping members would silently produce a wrong answer rather than
   * fail. Two elements coincide when `j * stride + k * dist` repeats, which
   * the separation of whichever axis is outer rules out.
   *
   * @param size How many elements each field holds.
   */
  bool Disjoint(Int size) const {
    if (_count <= 1 || size <= 0) return true;
    return _dist >= size * _stride || _stride >= _count * _dist;
  }

  /** @brief Compares componentwise. */
  bool operator==(const Batch&) const = default;

 private:
  Batch(Int count, Int stride, Int dist)
      : _count{count}, _stride{stride}, _dist{dist} {
    if (count < 1) {
      throw std::invalid_argument("Batch count must be at least one");
    }
    if (stride < 1 || dist < 1) {
      throw std::invalid_argument("Batch stride and dist must be positive");
    }
  }

  Int _count;
  Int _stride;
  Int _dist;
};

//-------------------------------------------------------------------------//
//                                Chunking                                  //
//--------------------------------------------------------------------------//

/**
 * @brief How many fields of a batch the inner loop takes at once.
 *
 * @details Not the caller's count. Batching measures a 2.3x gain with an
 * optimum around eight fields and is *worse than no batching at all* beyond
 * it, because the optimum is set by the chunk's coefficient block fitting in
 * last-level cache. A caller who batches a hundred radii must therefore not
 * have a hundred handed to the inner loop: Batch says what the caller has,
 * this says what is done with it.
 *
 * The figure that varies between machines is the cache, so that is what is
 * settable. Automatic() assumes a modest machine and is deliberately
 * conservative — undershooting forgoes some of the gain, overshooting is
 * slower than not batching — and a caller who knows their machine says so
 * with ForCache().
 */
class Chunking {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.

  /**
   * @brief A modest desktop's last-level cache, in bytes.
   * @details At @f$l_{\max} = 256@f$ this gives a chunk of three on one
   * thread, against the eight measured on a 16 MiB laptop — less gain, and no
   * risk of the reversal beyond the optimum.
   */
  static constexpr Int DefaultCacheBytes = Int{8} << 20;

  /**
   * @brief The largest chunk the heuristic may return.
   * @details A guard, not a measurement: at small degrees the formula below
   * grows without bound, and there the limit is per-call overhead rather than
   * cache. It also bounds the scratch a single call can ask for.
   */
  static constexpr Int MaximumCount = 64;

  /** @brief The heuristic, against DefaultCacheBytes. */
  static Chunking Automatic() { return Chunking(DefaultCacheBytes, 0); }

  /**
   * @brief The heuristic, against a stated cache.
   * @param bytes The total last-level cache of the machine. It is divided by
   * the number of copies of the block that will be live at once, since that is
   * what each of them gets.
   */
  static Chunking ForCache(Int bytes) {
    if (bytes < 1) {
      throw std::invalid_argument("Cache size must be positive");
    }
    return Chunking(bytes, 0);
  }

  /**
   * @brief A stated chunk, defeating the heuristic entirely.
   * @details What a benchmark sweeping chunk sizes needs, and what a caller
   * who has measured their own optimum should use.
   * @param count The chunk.
   */
  static Chunking Fixed(Int count) {
    if (count < 1) {
      throw std::invalid_argument("Chunk size must be at least one");
    }
    return Chunking(DefaultCacheBytes, count);
  }

  /**
   * @brief The chunk to use for one call.
   *
   * @details The divisor of two is headroom: the coefficient block is not
   * alone in the cache, the streamed Wigner row and the FFT output are there
   * too. The quotient is rounded to nearest rather than truncated because both
   * anchor points land just under an integer — 7.94 on a 16 MiB cache at
   * @f$l_{\max} = 256@f$, and 1.98 on a 256 MiB cache at 64 threads — so
   * truncating would give seven and *one*, and a chunk of one is the collapse
   * this formula exists to avoid.
   *
   * One thing it cannot see: on a multi-CCD or multi-socket machine a shared
   * block is pulled into each cache domain that touches it, so the inverse
   * transform's single copy is really one per domain. The rule is therefore
   * optimistic there in a way it is not on a single shared cache.
   *
   * @param bytesPerField The size of one field's coefficient block.
   * @param copies How many copies of that block will be live at once. This is
   * **not** the thread count in both directions: the forward transform gives
   * every thread a private accumulator, so its copies are its threads, while
   * the inverse gathers one block, shared and read-only, so it has one copy
   * however many threads read it. Serving both with the thread count starves
   * the inverse, which measures 2.2x slower where the whole batch would fit.
   */
  Int Count(Int bytesPerField, int copies) const {
    if (_fixed > 0) return _fixed;
    if (bytesPerField < 1) return MaximumCount;
    const auto share = _cacheBytes / static_cast<Int>(copies > 0 ? copies : 1);
    const auto count = (share + bytesPerField) / (2 * bytesPerField);
    if (count < 1) return 1;
    return count < MaximumCount ? count : MaximumCount;
  }

  /** @brief Compares componentwise. */
  bool operator==(const Chunking&) const = default;

 private:
  Chunking(Int cacheBytes, Int fixed)
      : _cacheBytes{cacheBytes}, _fixed{fixed} {}

  Int _cacheBytes;
  Int _fixed;  // zero means "use the heuristic"
};

//-------------------------------------------------------------------------//
//                     Stored or generated Wigner values                    //
//--------------------------------------------------------------------------//

/**
 * @brief Whether a grid stores its Wigner values or generates them when it
 * needs them.
 *
 * @details Stored() builds the whole @f$d^l_{nm}(\theta)@f$ table at
 * construction and streams it from every transform. That is 648 MB at
 * @f$l_{\max} = 256@f$ with @f$n_{\max} = 2@f$, and 43 GB at
 * @f$l_{\max} = 1024@f$, which is the first reason for the alternative. The
 * second is that streaming it is a *shared* cost — every thread pulls from the
 * same memory, and that ceiling is reached at around eight threads — whereas a
 * thread generating the values it needs contends with nobody.
 *
 * Generated() pays for that in arithmetic: the same recursion, the same
 * evaluation order and the same values, run inside the transform into
 * per-thread scratch instead of read from a table. Only the storage differs,
 * which is what makes the two paths comparable value by value rather than only
 * to a tolerance.
 *
 * This is a property of the grid and not of the call. A per-call choice would
 * mean carrying the table anyway for the calls that wanted it, which forfeits
 * the whole memory saving, and a template parameter would infect every
 * downstream type for a decision about one object's storage.
 */
class WignerValues {
 public:
  /** @brief Build the whole table at construction. */
  static WignerValues Stored() { return WignerValues(true); }
  /** @brief Run the recursion inside the transform, into per-thread scratch. */
  static WignerValues Generated() { return WignerValues(false); }

  /** @brief Whether the table is held. */
  auto AreStored() const { return _stored; }

  /** @brief Compares componentwise. */
  bool operator==(const WignerValues&) const = default;

 private:
  explicit WignerValues(bool stored) : _stored{stored} {}
  bool _stored;
};

//-------------------------------------------------------------------------//
//                            The Legendre kernel                           //
//--------------------------------------------------------------------------//

/**
 * @brief Which Legendre kernel a grid uses, and therefore how it lays out its
 * Wigner values.
 *
 * @details Loop() works a colatitude at a time, an `axpy` of length @f$k@f$
 * per @f$(l, m)@f$, reading a table contiguous in @f$(l, m)@f$ at fixed
 * @f$(n, \theta)@f$. Matrix() is transform-major: all the FFTs first, then one
 * matrix product per order against a table contiguous in @f$(l, \theta)@f$ at
 * fixed @f$(n, m)@f$. Products at different orders write disjoint outputs, so
 * it needs no accumulator and no reduction in either direction — which is the
 * forward transform's weak point at high thread counts, deleted rather than
 * tuned.
 *
 * **Both are kept, permanently.** This is not a migration with a flag day: the
 * point of the two coexisting is measurement. A GEMM sums in whatever order
 * its kernel chooses, so the matrix path cannot be bit-compared against the
 * loop path — but it can be compared to a tolerance, on identical inputs,
 * which is a check of the layout, the indexing, the FFT ordering and the
 * accumulation that no single-kernel library has available. It also means the
 * right kernel for a machine is a question that machine can answer for itself.
 *
 * A property of the grid rather than of the call, for the reason WignerValues
 * is: the two layouts are mutually exclusive, and 648 MB at
 * @f$l_{\max} = 256@f$ with @f$n_{\max} = 2@f$ is not a thing to hold twice.
 *
 * Matrix() is *absent*, not merely refused, in a build without BLAS, so asking
 * for it there is a compile error at the call site rather than a throw at grid
 * construction. That is the tighter of the two options and the reversible one:
 * offering it later with a runtime throw breaks nobody, while withdrawing it
 * would.
 *
 * Matrix() with WignerValues::Generated() is refused at grid construction, and
 * that is a fact about the recursion rather than an unimplemented case. The
 * recursion's output for one @f$(n, \theta)@f$ spans every order at once, so
 * isolating the single @f$(n, m)@f$ block a per-order product wants means
 * either keeping all of it — which is the table, and not having one is the
 * entire purpose of generating — or re-running the recursion once per order.
 */
class TransformKernel {
 public:
  /** @brief A colatitude at a time, over a table contiguous in @f$(l, m)@f$. */
  static TransformKernel Loop() { return TransformKernel(false); }

#ifdef GSHTRANS_HAVE_BLAS
  /**
   * @brief One matrix product per order, over a table contiguous in
   * @f$(l, \theta)@f$. Present only in a build with BLAS.
   */
  static TransformKernel Matrix() { return TransformKernel(true); }
#endif

  /** @brief Whether the matrix kernel was asked for. */
  auto IsMatrix() const { return _matrix; }

  /** @brief Compares componentwise. */
  bool operator==(const TransformKernel&) const = default;

 private:
  explicit TransformKernel(bool matrix) : _matrix{matrix} {}
  bool _matrix;
};

//-------------------------------------------------------------------------//
//                          The interpolation scheme                        //
//--------------------------------------------------------------------------//

/**
 * @brief Which scheme an interpolant of a field uses.
 *
 * @details Spectral() evaluates the expansion directly, which is exact for a
 * band-limited field and safe at the poles, and costs @f$O(l_{\max}^2)@f$ a
 * point. Bilinear() and Bicubic() are rectilinear schemes over the field's own
 * samples, which cost @f$O(1)@f$ a point and are not exact; they are *absent*,
 * not refused, in a build without the interpolation dependency.
 *
 * This policy differs from the others here in one respect. They are values
 * read at run time by code that is the same either way. These three select
 * between interpolants of *different types* — one holds an expansion, the
 * others hold upstream objects with their own coefficient arrays — so
 * Interpolate must dispatch at compile time and return a scheme-dependent
 * type. The named constructors therefore return distinct tag types rather than
 * a common value type, which keeps the caller's spelling identical to the
 * other policies while making the dispatch overload resolution.
 *
 * What that costs: a scheme cannot be chosen from a configuration file without
 * a switch in the caller's own code.
 */
class Scheme {
 public:
  struct SpectralTag {};  ///< Selects direct evaluation of the expansion.
  struct BilinearTag {};  ///< Selects bilinear interpolation of the samples.
  struct BicubicTag {};   ///< Selects bicubic interpolation of the samples.

  /** @brief Exact for a band-limited field, and safe at the poles. */
  static constexpr SpectralTag Spectral() { return {}; }

#ifdef GSHTRANS_HAVE_INTERPOLATION
  /** @brief Bilinear over the field's samples. Needs the interpolation
   * dependency. */
  static constexpr BilinearTag Bilinear() { return {}; }
  /** @brief Bicubic over the field's samples. Needs the interpolation
   * dependency. */
  static constexpr BicubicTag Bicubic() { return {}; }
#endif
};

}  // namespace GSHTrans

#endif  //  GSH_TRANS_POLICIES_GUARD_H
