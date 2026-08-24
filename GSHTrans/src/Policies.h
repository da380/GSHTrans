#ifndef GSH_TRANS_POLICIES_GUARD_H
#define GSH_TRANS_POLICIES_GUARD_H

// The four policy classes: how an operation may thread, how a caller's batch
// is laid out, how much of one the inner loop takes at a time, and whether a
// grid stores its Wigner values or generates them.
//
// They lived in Concepts.h, which had outgrown its name -- these are not
// concepts, and three of the four carry measured performance reasoning that
// made the concepts hard to find among it. What they have in common is that
// each is a *decision the caller makes about the machine* rather than about
// the mathematics, which is why they belong together and apart.
//
// Values rather than macros or template parameters, deliberately: the library
// is header-only and a macro would have to be defined identically in every
// translation unit, while a template parameter would put a machine property
// into a type that the mathematics has no business carrying.

#include <algorithm>
#include <cstddef>
#include <stdexcept>

namespace GSHTrans {

//-------------------------------------------------------------------------//
//                            Execution policy                              //
//--------------------------------------------------------------------------//

// Whether an operation may use threads, and how many.
//
// Sequential by default everywhere: a library should not create threads
// because it can, only because it was asked to. The rule this exists to make
// keepable is that exactly one level threads -- a caller parallelising over
// slices, components or realisations calls the transform sequentially, while a
// caller with one large problem asks the transform to thread. Both at once is
// worse than either, so an operation asked to run in parallel from inside an
// existing parallel region runs sequentially instead.
//
// A thread count of zero means "whatever OpenMP would choose", which respects
// OMP_NUM_THREADS. On a machine with simultaneous multithreading that is
// usually the number of hardware threads, and for this library's memory-bound
// work that measures *slower* than using one thread per core, so callers who
// care should say what they want.
class Execution {
 public:
  static Execution Sequential() { return Execution(1); }
  static Execution Parallel(int threads = 0) {
    return Execution(threads > 0 ? threads : 0);
  }

  auto Threads() const { return _threads; }
  auto IsParallel() const { return _threads != 1; }

  /** @brief Compares componentwise. */
  bool operator==(const Execution&) const = default;

 private:
  explicit Execution(int threads) : _threads{threads} {}
  int _threads;
};

// How a batch of same-spin fields is laid out in the caller's storage.
//
// The batched transform's primitive is "transform k fields sharing a grid, a
// degree and an upper index", with the single field as k = 1. What varies
// between callers is not how many fields there are but how they are arranged,
// so the batch is described by (count, stride, dist) rather than by
// contiguity: element j of field k lives at j * stride + k * dist.
//
// This is FFTW's advanced-interface form, and it is used because it costs
// nothing at the stage that actually reads the caller's layout -- the FFT --
// while covering the two arrangements the field layer produces:
//
//   radial slab, or tensor components stored component by component
//                                      stride = 1,           dist = FieldSize
//   tensor components stored point by point
//                                      stride = nComponents, dist = 1
//
// The second row is why this is a decision rather than a detail. Without it
// the field layer would have to promise component-major storage and repack
// anything else before transforming; with it the requirement disappears
// (core-plan.md [C9]). Whether the strided path is *faster* than repacking is
// a separate, measurable question. What changed is that the interface no
// longer forces an answer to it.
//
// Full FFTW nembed generality is deliberately not offered: the angular layout
// is fixed by the grid, so most of it would be unreachable.
//
// A batch shares grid, degree and upper index. The Wigner block for that
// upper index is precisely what batching amortises, so fields at different
// upper indices cannot batch together -- for a rank-2 tensor that means
// batching over radii within each n, not across components. It is worth
// saying plainly, because "batch all my fields" is the natural expectation
// and it is the wrong one.
class Batch {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.

  // Fields laid end to end, each one contiguous and `size` elements long.
  static Batch Contiguous(Int count, Int size) {
    return Batch(count, 1, size);
  }

  // Fields interleaved element by element, `stride` of them per element. The
  // count may be smaller than the stride: a batch of two components out of a
  // point-major five is Interleaved(2, 5), and the other three are not part
  // of the call.
  static Batch Interleaved(Int count, Int stride) {
    return Batch(count, stride, 1);
  }

  // The general form, for a caller whose layout is neither.
  static Batch Strided(Int count, Int stride, Int dist) {
    return Batch(count, stride, dist);
  }

  // The single field, which is what the unbatched entry points hand down.
  static Batch One(Int size) { return Batch(1, 1, size); }

  auto Count() const { return _count; }
  auto Stride() const { return _stride; }
  auto Dist() const { return _dist; }

  // Where element j of field k lives.
  Int Offset(Int j, Int k) const { return j * _stride + k * _dist; }

  // The smallest range that holds this batch when each field has `size`
  // elements. The batched entry points check the caller's range against this
  // rather than for equality: an interleaved batch is a window onto a larger
  // range whose other elements are none of the transform's business, so
  // demanding an exact size would reject exactly the layout this type exists
  // to accept. The unbatched entry points keep the equality check, since
  // their contract is that the range *is* the field.
  Int Span(Int size) const {
    return (size - 1) * _stride + (_count - 1) * _dist + 1;
  }

  // Whether the fields of this batch are disjoint at that size.
  //
  // They must be: the transform writes every element of every field, so
  // overlapping members would silently produce a wrong answer rather than
  // fail. Two elements coincide when j * stride + k * dist repeats, which the
  // separation of whichever axis is outer rules out.
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

// How many fields of a batch the inner loop takes at once.
//
// Not the caller's count. Tier-1 batching was measured at 2.3x with an optimum
// around eight fields and *worse than no batching at all* beyond it
// (core-plan.md P2), because the optimum is set by the chunk's coefficient
// block fitting in last-level cache. A caller who batches a hundred radii must
// therefore not have a hundred handed to the inner loop: the descriptor says
// what the caller has, this says what is done with it.
//
// The figure that varies between machines is the cache, so that is what is
// settable. `Automatic` assumes a modest machine and is deliberately
// conservative -- undershooting forgoes some of the 2.3x, overshooting is
// slower than not batching, so the default errs small and a caller who knows
// their machine says so with `ForCache`. `Fixed` defeats the heuristic
// entirely, which is what a benchmark sweeping chunk sizes needs, and what a
// caller who has measured their own optimum should use.
//
// This is a value rather than a macro on purpose. The library is header-only,
// so a preprocessor knob defined differently in two translation units would
// give GaussLegendreGrid two inline bodies and let the linker pick one in
// silence; it could not be swept without a rebuild, which is how measurements
// stop being taken; and it could not read anything about the machine.
class Chunking {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.

  // A modest desktop's last-level cache. At lMax = 256 this gives a chunk of
  // three on one thread, against the eight P2 measured on a 16 MiB laptop --
  // less gain, and no risk of the reversal beyond the optimum.
  static constexpr Int DefaultCacheBytes = Int{8} << 20;

  // A guard, not a measurement: at small degrees the formula below grows
  // without bound, and there the limit is per-call overhead rather than
  // cache. It also bounds the scratch a single call can ask for.
  static constexpr Int MaximumCount = 64;

  static Chunking Automatic() { return Chunking(DefaultCacheBytes, 0); }

  // The total last-level cache of the machine, in bytes. Divided by the number
  // of copies of the block that will be live at once, since that is what each
  // of them gets.
  static Chunking ForCache(Int bytes) {
    if (bytes < 1) {
      throw std::invalid_argument("Cache size must be positive");
    }
    return Chunking(bytes, 0);
  }

  static Chunking Fixed(Int count) {
    if (count < 1) {
      throw std::invalid_argument("Chunk size must be at least one");
    }
    return Chunking(DefaultCacheBytes, count);
  }

  // The chunk to use for a call whose coefficient block is `bytesPerField`
  // bytes and which will hold `copies` of that block at once.
  //
  // The divisor of two is headroom: the coefficient block is not alone in the
  // cache, the streamed Wigner row and the FFT output are there too.
  //
  // `copies` is the count that matters, and it is **not** the thread count in
  // both directions. The forward transform gives every thread a private
  // accumulator of chunk * coefficientSize, so its copies are its threads; the
  // inverse gathers one block, shared and read-only, so it has one copy
  // however many threads read it. Serving both with the thread count starves
  // the inverse: at lMax = 256, k = 8 on eight threads it returns a chunk of
  // one where the whole batch fits, which measured 2.2x slower (section 10).
  //
  // Both of P8's anchors still land, since both were forward-shaped. P2's
  // optimum of eight was measured *sequentially*, so one copy had the whole of
  // a 16 MiB cache: 16 MiB / 1 gives eight, while the 2 MiB per-core share of
  // the same machine would give one. On a 256 MiB, 64-core machine at full
  // width, sixty-four copies give two, which is what P8 predicted.
  //
  // One thing this cannot see, and the target-machine run should. On a
  // multi-CCD or multi-socket machine a *shared* block is pulled into each
  // cache domain that touches it, so the inverse's single copy is really one
  // per domain -- eight CCDs on the deployment target. The rule below is
  // therefore optimistic there in a way it is not on a single shared L3.
  Int Count(Int bytesPerField, int copies) const {
    if (_fixed > 0) return _fixed;
    if (bytesPerField < 1) return MaximumCount;
    const auto share = _cacheBytes / static_cast<Int>(copies > 0 ? copies : 1);
    // Rounded to nearest, not truncated. Both of P8's anchors land just under
    // an integer -- 7.94 on a 16 MiB cache at lMax = 256, and 1.98 on a
    // 256 MiB cache at 64 threads -- so truncation would give seven and *one*,
    // and the second is exactly the collapse to a chunk of one that P8 raised
    // the formula to avoid. Adding half before dividing gives eight and two,
    // which are the numbers that document quotes.
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

// Whether a grid stores its Wigner values or generates them when it needs
// them.
//
// Stored is what the library has always done: the whole d^l_{nm}(theta) table
// is built at construction and streamed by every transform. It is 648 MB at
// lMax = 256 with nMax = 2, and 43 GB at lMax = 1024, which is the first
// reason for the alternative. The second is that streaming it is a shared
// cost -- every thread pulls from the same DRAM, and step H measured that
// ceiling being reached at eight threads -- whereas a thread that generates
// the values it needs contends with nobody, which is the shape the deployment
// target wants (core-plan.md step F', P8).
//
// Generated pays for that in arithmetic: the same recursion, the same
// evaluation order and the same values, run inside the transform into
// per-thread scratch instead of read from a table. Only the storage differs,
// which is what makes the two paths comparable value by value rather than
// only to a tolerance.
//
// This is a property of the grid and not of the call. A per-call choice would
// mean carrying the table anyway for the calls that wanted it, which forfeits
// the whole memory saving, and a template parameter would infect every
// downstream type for a decision about one object's storage ([C10]). It is a
// value rather than a named constructor for the reason Chunking is: ForBand
// is itself a named constructor, so the named form needs a second one to
// reach an oversampled generating grid.
class WignerValues {
 public:
  static WignerValues Stored() { return WignerValues(true); }
  static WignerValues Generated() { return WignerValues(false); }

  auto AreStored() const { return _stored; }

  /** @brief Compares componentwise. */
  bool operator==(const WignerValues&) const = default;

 private:
  explicit WignerValues(bool stored) : _stored{stored} {}
  bool _stored;
};

// Which Legendre kernel a grid uses, and therefore how it lays out its Wigner
// values.
//
// Loop is the kernel of core-plan.md steps F and H, and is what the library
// has always done: a colatitude at a time, an axpy of length k per (l, m),
// reading a table contiguous in (l, m) at fixed (n, theta).
//
// Matrix is the transform-major restructure of section 11: all the FFTs
// first, then one matrix product per order against a table contiguous in
// (l, theta) at fixed (n, m). Products at different orders write disjoint
// outputs, so it needs no accumulator and no reduction in either direction --
// which is the forward transform's weak point at high thread counts, deleted
// rather than tuned.
//
// **Both are kept, permanently** ([C12]). This is not a migration with a flag
// day: the two coexist, and the point of that is measurement. A GEMM sums in
// whatever order its kernel chooses, so the matrix path cannot be
// bit-compared against the loop path -- but it can be compared to a
// tolerance, on identical inputs, which is a check of the layout, the
// indexing, the FFT ordering and the accumulation that no single-kernel
// library has available. It also means the right kernel for a machine is a
// question that machine can answer for itself.
//
// A property of the grid rather than of the call, for the reason WignerValues
// is: the two layouts are mutually exclusive -- 648 MB at lMax = 256 with
// nMax = 2 is not a thing to hold twice -- and a template parameter would put
// a machine decision into every downstream type.
//
// -- Matrix is absent, not merely refused, in a build without BLAS ([C17]).
// The factory below does not exist there, so asking for it is a compile error
// at the call site rather than a throw at grid construction. That is the
// tighter of the two options and the reversible one: offering it later with a
// runtime throw breaks nobody, while withdrawing it would.
//
// -- Matrix with WignerValues::Generated() is refused at grid construction,
// and that is a fact about the recursion rather than an unimplemented case.
// The recursion's output for one (n, theta) spans every order at once, so
// isolating the single (n, m) block a per-order product wants means either
// keeping all of it -- which is the table, and not having one is the entire
// purpose of generating -- or re-running the recursion once per order. There
// is no third way.
class TransformKernel {
 public:
  static TransformKernel Loop() { return TransformKernel(false); }

#ifdef GSHTRANS_HAVE_BLAS
  static TransformKernel Matrix() { return TransformKernel(true); }
#endif

  auto IsMatrix() const { return _matrix; }

  /** @brief Compares componentwise. */
  bool operator==(const TransformKernel&) const = default;

 private:
  explicit TransformKernel(bool matrix) : _matrix{matrix} {}
  bool _matrix;
};

// Which scheme an interpolant of a field uses.
//
// Spectral evaluates the expansion directly, which is exact for a
// band-limited field and safe at the poles, and costs O(lMax^2) a point.
// Bilinear and Bicubic are Interpolation's rectilinear schemes over the
// field's own samples, which cost O(1) a point and are not exact.
//
// -- Bilinear and Bicubic are **absent**, not refused, in a build without
// Interpolation, so asking for one is a compile error at the call site rather
// than a throw. That is [C17]'s argument for the third time in this library,
// after TransformKernel::Matrix() and RadialSplineDerivative.h: the tighter
// option is the reversible one, since offering a factory later with a runtime
// throw breaks nobody while withdrawing one does.
//
// This policy differs from the four above in one respect, and it is worth
// naming rather than leaving to be discovered. The others are values read at
// run time by code that is the same either way. These three select between
// interpolants of *different types* -- one holds an expansion, the others hold
// upstream objects with their own coefficient arrays -- so Interpolate must
// dispatch at compile time and return a scheme-dependent type. The named
// constructors below therefore return distinct tag types rather than a common
// value type, which keeps the caller's spelling identical to the other
// policies while making the dispatch overload resolution.
//
// What that costs: a scheme cannot be chosen from a configuration file
// without a switch in the caller's own code. It is the same limitation
// [C17] records for TransformKernel and it is accepted for the same reason.
class Scheme {
 public:
  struct SpectralTag {};
  struct BilinearTag {};
  struct BicubicTag {};

  static constexpr SpectralTag Spectral() { return {}; }

#ifdef GSHTRANS_HAVE_INTERPOLATION
  static constexpr BilinearTag Bilinear() { return {}; }
  static constexpr BicubicTag Bicubic() { return {}; }
#endif
};

}  // namespace GSHTrans

#endif  //  GSH_TRANS_POLICIES_GUARD_H
