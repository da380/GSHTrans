#ifndef GSH_TRANS_CONCEPTS_GUARD_H
#define GSH_TRANS_CONCEPTS_GUARD_H

#include <concepts>

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
