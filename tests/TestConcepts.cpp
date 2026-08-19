#include <gtest/gtest.h>

#include <GSHTrans/All>
#include <complex>
#include <span>
#include <vector>

// Pins the concept surface: what the numeric concepts answer, and exactly
// which argument types the transforms accept.
//
// These are compile-time facts, so the test body is nearly empty; the value is
// that the static_asserts below are checked whenever this file is built. They
// exist because the numeric concepts are defined by NumericConcepts rather
// than here, and because the transforms' requires-clauses are the library's
// only defence against a caller handing over a range of the wrong scalar type.

namespace {

using namespace GSHTrans;

using Real = double;
using Complex = std::complex<Real>;

//--------------------------------------------------------------------------//
//                            The numeric concepts                           //
//--------------------------------------------------------------------------//

static_assert(RealFloatingPoint<float>);
static_assert(RealFloatingPoint<double>);
static_assert(RealFloatingPoint<long double>);
static_assert(!RealFloatingPoint<int>);
static_assert(!RealFloatingPoint<std::complex<double>>);

static_assert(ComplexFloatingPoint<std::complex<float>>);
static_assert(ComplexFloatingPoint<std::complex<double>>);
static_assert(ComplexFloatingPoint<std::complex<long double>>);
static_assert(!ComplexFloatingPoint<std::complex<int>>);
static_assert(!ComplexFloatingPoint<double>);
static_assert(!ComplexFloatingPoint<int>);

static_assert(RealOrComplexFloatingPoint<double>);
static_assert(RealOrComplexFloatingPoint<std::complex<double>>);
static_assert(!RealOrComplexFloatingPoint<int>);

static_assert(std::same_as<RemoveComplex<double>, double>);
static_assert(std::same_as<RemoveComplex<std::complex<double>>, double>);
static_assert(std::same_as<RemoveComplex<long double>, long double>);
static_assert(
    std::same_as<RemoveComplex<std::complex<long double>>, long double>);

// The tag concepts are this library's own and are unaffected by any of the
// above.
static_assert(RealOrComplexValued<RealValued>);
static_assert(RealOrComplexValued<ComplexValued>);
static_assert(!RealOrComplexValued<double>);
static_assert(IndexRange<All> && IndexRange<NonNegative> && IndexRange<Single>);
static_assert(!IndexRange<Multiple>);
static_assert(OrderIndexRange<All> && OrderIndexRange<NonNegative>);
static_assert(!OrderIndexRange<Single>);
static_assert(AngleIndexRange<Single> && AngleIndexRange<Multiple>);
static_assert(!AngleIndexRange<All>);

//--------------------------------------------------------------------------//
//                        What the transforms accept                         //
//--------------------------------------------------------------------------//

template <typename Grid, typename In, typename Out>
concept AdmitsForward = requires(Grid grid, In in, Out out) {
  grid.ForwardTransformation(std::ptrdiff_t{1}, std::ptrdiff_t{0}, in, out);
};

template <typename Grid, typename In, typename Out>
concept AdmitsInverse = requires(Grid grid, In in, Out out) {
  grid.InverseTransformation(std::ptrdiff_t{1}, std::ptrdiff_t{0}, in, out);
};

using Grid = GaussLegendreGrid<Real, All, All>;
using ScalarGrid = GaussLegendreGrid<Real, NonNegative, All>;

template <typename T>
using Vec = std::vector<T>;

// The two supported field scalars, into and out of complex coefficients.
static_assert(AdmitsForward<Grid, Vec<Real>, Vec<Complex>>);
static_assert(AdmitsForward<Grid, Vec<Complex>, Vec<Complex>>);
static_assert(AdmitsInverse<Grid, Vec<Complex>, Vec<Real>>);
static_assert(AdmitsInverse<Grid, Vec<Complex>, Vec<Complex>>);

// Coefficients are always complex, and always of the grid's precision.
static_assert(!AdmitsForward<Grid, Vec<Real>, Vec<Real>>);
static_assert(!AdmitsForward<Grid, Vec<Real>, Vec<std::complex<float>>>);
static_assert(!AdmitsInverse<Grid, Vec<Real>, Vec<Real>>);

// The field must share the grid's precision.
static_assert(!AdmitsForward<Grid, Vec<float>, Vec<Complex>>);
static_assert(!AdmitsForward<Grid, Vec<std::complex<float>>, Vec<Complex>>);
static_assert(!AdmitsInverse<Grid, Vec<Complex>, Vec<long double>>);

// Non-numeric ranges are not fields.
static_assert(!AdmitsForward<Grid, Vec<int>, Vec<Complex>>);

// The output must be writable; the input need not be.
static_assert(AdmitsForward<Grid, std::span<const Real>, Vec<Complex>>);
static_assert(!AdmitsForward<Grid, Vec<Real>, std::span<const Complex>>);
static_assert(!AdmitsInverse<Grid, Vec<Complex>, std::span<const Real>>);

// A grid storing only non-negative orders serves real fields alone: the
// complex path does not compile there, rather than throwing (core-plan step A).
static_assert(AdmitsForward<ScalarGrid, Vec<Real>, Vec<Complex>>);
static_assert(!AdmitsForward<ScalarGrid, Vec<Complex>, Vec<Complex>>);
static_assert(AdmitsInverse<ScalarGrid, Vec<Complex>, Vec<Real>>);
static_assert(!AdmitsInverse<ScalarGrid, Vec<Complex>, Vec<Complex>>);

}  // namespace

TEST(Concepts, CompileTimeSurfaceIsPinned) { SUCCEED(); }
