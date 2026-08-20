#ifndef GSH_TRANS_CONCEPTS_GUARD_H
#define GSH_TRANS_CONCEPTS_GUARD_H

#include <concepts>
#include <cstddef>
#include <stdexcept>

#include <NumericConcepts/Numeric.hpp>

namespace GSHTrans {

//-------------------------------------------------------------------------//
//                                 Tag classes                              //
//--------------------------------------------------------------------------//

// Matrix storage options for Wigner class.
struct ColumnMajor {};
struct RowMajor {};

template <typename T>
concept WignerStorage =
    std::same_as<T, ColumnMajor> or std::same_as<T, RowMajor>;

// Index storage options.
struct All {};
struct NonNegative {};
struct Single {};
struct Multiple {};

struct UpperIndexFirst {};
struct AngleFirst {};

template <typename Indices>
concept IndexRange =
    std::same_as<Indices, All> or std::same_as<Indices, NonNegative> or
    std::same_as<Indices, Single>;

template <typename Indices>
concept OrderIndexRange =
    std::same_as<Indices, All> or std::same_as<Indices, NonNegative>;

template <typename Indices>
concept AngleIndexRange =
    std::same_as<Indices, Multiple> or std::same_as<Indices, Single>;

// Value type options.
struct RealValued {};
struct ComplexValued {};

template <typename T>
concept RealOrComplexValued =
    std::same_as<T, RealValued> or std::same_as<T, ComplexValued>;

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
  using Int = std::ptrdiff_t;

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
//                      Numeric concepts, from elsewhere                    //
//--------------------------------------------------------------------------//

// These are NumericConcepts' definitions rather than this library's. The point
// of that dependency is that the projects built on it agree on what "real",
// "complex" and "precision" mean, and a second, subtly different set of
// definitions here would defeat it -- the two did in fact differ, over whether
// a const-qualified std::complex counts as complex.
//
// They are re-spelled rather than imported unqualified because NumericConcepts
// names them Real and Complex, while this library uses Real and Complex as
// member type aliases in nearly every class. Importing those names would have
// them shadowed at exactly the points where the concept is most likely to be
// wanted, and the resulting errors would be obscure. Anything new should use
// the NumericConcepts spelling directly, qualified.
template <typename T>
concept RealFloatingPoint = NumericConcepts::Real<T>;

template <typename T>
concept ComplexFloatingPoint = NumericConcepts::Complex<T>;

template <typename T>
concept RealOrComplexFloatingPoint = NumericConcepts::RealOrComplex<T>;

using NumericConcepts::RemoveComplex;

// Concept for scalar-valued function on S2.
template <typename Function, typename Real, typename Scalar>
concept ScalarFunctionS2 = requires(Function f, Real theta, Real phi) {
  requires RealFloatingPoint<Real>;
  requires RealOrComplexFloatingPoint<Scalar>;
  { f(theta, phi) } -> std::convertible_to<Scalar>;
};

// Concept for scalar-valued function of spherical harmonic indices.
template <typename Function, typename Int, typename Scalar>
concept ScalarFunctionS2Expansion = requires(Function f, Int l, Int m) {
  requires std::integral<Int>;
  requires RealOrComplexFloatingPoint<Scalar>;
  { f(l, m) } -> std::convertible_to<Scalar>;
};

}  // namespace GSHTrans

#endif  //  GSH_TRANS_CONCEPTS_GUARD_H
