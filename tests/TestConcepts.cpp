#include <gtest/gtest.h>

#include <GSHTrans/All>
#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
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

using Int = std::ptrdiff_t;
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

//--------------------------------------------------------------------------//
//                          The batch descriptor                             //
//--------------------------------------------------------------------------//
//
// The layouts of core-plan.md [C9], checked as addresses rather than as
// accessors: what matters about (count, stride, dist) is where element j of
// field k actually lands, so the tests enumerate that.

TEST(Batch, ContiguousLaysFieldsEndToEnd) {
  constexpr auto size = Int{4};
  const auto batch = Batch::Contiguous(3, size);

  EXPECT_EQ(batch.Count(), 3);
  EXPECT_EQ(batch.Stride(), 1);
  EXPECT_EQ(batch.Dist(), size);

  auto seen = std::vector<Int>{};
  for (auto k = Int{0}; k < batch.Count(); ++k) {
    for (auto j = Int{0}; j < size; ++j) seen.push_back(batch.Offset(j, k));
  }
  EXPECT_EQ(seen, (std::vector<Int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
}

TEST(Batch, InterleavedPutsOnePointsComponentsTogether) {
  constexpr auto size = Int{3};
  const auto batch = Batch::Interleaved(2, 5);

  EXPECT_EQ(batch.Stride(), 5);
  EXPECT_EQ(batch.Dist(), 1);

  // Two components out of a point-major five: the call touches columns 0 and
  // 1 of each point and nothing else. That is the layout the field plan would
  // otherwise have had to repack.
  auto seen = std::vector<Int>{};
  for (auto k = Int{0}; k < batch.Count(); ++k) {
    for (auto j = Int{0}; j < size; ++j) seen.push_back(batch.Offset(j, k));
  }
  EXPECT_EQ(seen, (std::vector<Int>{0, 5, 10, 1, 6, 11}));
}

TEST(Batch, OneIsTheUnbatchedCase) {
  const auto batch = Batch::One(7);
  EXPECT_EQ(batch.Count(), 1);
  EXPECT_EQ(batch.Span(7), 7);
  for (auto j = Int{0}; j < 7; ++j) EXPECT_EQ(batch.Offset(j, 0), j);
}

TEST(Batch, SpanIsTheRangeTheCallActuallyNeeds) {
  // Contiguous: exactly the whole of it.
  EXPECT_EQ(Batch::Contiguous(3, 4).Span(4), 12);

  // Interleaved: a window onto a larger range. Two components of a five-wide
  // point-major array over three points reach index 11, not 15 -- the last
  // point's remaining three components are not part of the call, and
  // demanding an exact size would reject the layout.
  EXPECT_EQ(Batch::Interleaved(2, 5).Span(3), 12);
}

TEST(Batch, DisjointSeparatesLayoutsThatOverlap) {
  // Both standard layouts are disjoint, by separation of opposite axes.
  EXPECT_TRUE(Batch::Contiguous(3, 4).Disjoint(4));
  EXPECT_TRUE(Batch::Interleaved(5, 5).Disjoint(3));
  EXPECT_TRUE(Batch::One(4).Disjoint(4));

  // A stride that does not clear the count aliases: with stride 2 and dist 1,
  // element 1 of field 0 and element 0 of field 2 are both at index 2.
  const auto overlapping = Batch::Strided(3, 2, 1);
  EXPECT_EQ(overlapping.Offset(1, 0), overlapping.Offset(0, 2));
  EXPECT_FALSE(overlapping.Disjoint(3));

  // Likewise the other way round: fields closer together than they are long.
  EXPECT_FALSE(Batch::Strided(2, 1, 3).Disjoint(4));
  EXPECT_TRUE(Batch::Strided(2, 1, 4).Disjoint(4));
}

TEST(Batch, RejectsDescriptorsThatCannotDescribeAnything) {
  EXPECT_THROW(Batch::Contiguous(0, 4), std::invalid_argument);
  EXPECT_THROW(Batch::Contiguous(-1, 4), std::invalid_argument);
  EXPECT_THROW(Batch::Interleaved(2, 0), std::invalid_argument);
  EXPECT_THROW(Batch::Strided(2, 1, 0), std::invalid_argument);
  EXPECT_THROW(Batch::Strided(2, -1, 1), std::invalid_argument);
}

TEST(Batch, EqualityIsMemberwise) {
  EXPECT_EQ(Batch::Contiguous(3, 4), Batch::Strided(3, 1, 4));
  EXPECT_NE(Batch::Contiguous(3, 4), Batch::Contiguous(3, 5));
}

//--------------------------------------------------------------------------//
//                          The chunking policy                              //
//--------------------------------------------------------------------------//
//
// The coefficient block of one field at (lMax, n) holds (lMax+1)^2 - n^2
// complex values, so at double precision:
//
//   lMax = 128, n = 2   16637 coefficients     266192 bytes
//   lMax = 256, n = 2   66045 coefficients    1056720 bytes
//   lMax = 512, n = 2  263165 coefficients    4210640 bytes
//
// Those are the sizes the numbers below are anchored on.

namespace {
constexpr auto bytesAt128 = Int{266192};
constexpr auto bytesAt256 = Int{1056720};
constexpr auto bytesAt512 = Int{4210640};
constexpr auto mebibyte = Int{1} << 20;
}  // namespace

TEST(Chunking, ReproducesTheTwoMeasuredAnchors) {
  // P2's optimum of eight was measured on a 16 MiB laptop at lMax = 256,
  // running sequentially -- so one thread had the whole cache.
  EXPECT_EQ(Chunking::ForCache(16 * mebibyte).Count(bytesAt256, 1), 8);

  // P8's prediction for a 256 MiB, 64-core machine at full width. The point
  // of dividing by the running threads rather than by a fixed per-core figure
  // is that both of these come out right.
  EXPECT_EQ(Chunking::ForCache(256 * mebibyte).Count(bytesAt256, 64), 2);

  // And the collapse P8 warned about, if the same cache is divided by
  // hardware threads instead of by the cores actually doing the work.
  EXPECT_EQ(Chunking::ForCache(256 * mebibyte).Count(bytesAt256, 128), 1);
}

TEST(Chunking, TheDefaultIsConservative) {
  const auto automatic = Chunking::Automatic();

  // Below the eight measured on a 16 MiB machine, because overshooting the
  // optimum is slower than not batching while undershooting only forgoes
  // gain.
  EXPECT_EQ(automatic.Count(bytesAt256, 1), 4);
  EXPECT_EQ(automatic.Count(bytesAt128, 1), 16);

  // At lMax = 512 one field's coefficients already fill the assumed cache, so
  // there is nothing to batch and the floor holds.
  EXPECT_EQ(automatic.Count(bytesAt512, 1), 1);

  // Never zero, however little cache is assumed to be going spare.
  EXPECT_EQ(automatic.Count(bytesAt512, 64), 1);
  EXPECT_EQ(Chunking::ForCache(1).Count(bytesAt256, 1), 1);
}

TEST(Chunking, FixedDefeatsTheHeuristic) {
  const auto fixed = Chunking::Fixed(3);
  EXPECT_EQ(fixed.Count(bytesAt128, 1), 3);
  EXPECT_EQ(fixed.Count(bytesAt512, 64), 3);
  EXPECT_EQ(fixed.Count(0, 1), 3);
}

TEST(Chunking, IsBoundedAboveAtSmallDegrees) {
  // The formula grows without bound as the degree falls, where the limit is
  // per-call overhead rather than cache. The cap is a guard, not a measured
  // optimum.
  EXPECT_EQ(Chunking::Automatic().Count(16, 1), Chunking::MaximumCount);
}

TEST(Chunking, RejectsPoliciesThatDescribeNothing) {
  EXPECT_THROW(Chunking::ForCache(0), std::invalid_argument);
  EXPECT_THROW(Chunking::ForCache(-1), std::invalid_argument);
  EXPECT_THROW(Chunking::Fixed(0), std::invalid_argument);
  EXPECT_THROW(Chunking::Fixed(-2), std::invalid_argument);
}

//--------------------------------------------------------------------------//
//                              AngularGrid                                  //
//--------------------------------------------------------------------------//
//
// The concept gained the two axes when SphericalGrid took over the transform:
// grids on this library's terms are separable, so Interpolate's local
// refinement asking for them separately was describing every grid there is
// and went away (core-plan.md section 13, G2).
//
// Asserted both ways, because a requirement nothing can fail is decoration. A
// type carrying everything the concept asked for *before* that change must now
// be rejected, and the real grid must still be accepted.

namespace {

using ProbeReal = double;

// Everything AngularGrid required before the axes were added.
struct WithoutAxes {
  using Real = ProbeReal;
  using NRange = GSHTrans::All;
  const void* Identity() const { return this; }
  std::ptrdiff_t FieldSize() const { return 0; }
  std::ptrdiff_t MaxUpperIndex() const { return 0; }
  std::vector<std::ptrdiff_t> UpperIndices() const { return {}; }
  std::vector<std::pair<Real, Real>> Points() const { return {}; }
  std::vector<std::ptrdiff_t> CoLatitudeIndices() const { return {}; }
  std::vector<std::ptrdiff_t> LongitudeIndices() const { return {}; }
  std::vector<Real> CoLatitudeWeights() const { return {}; }
  std::vector<Real> LongitudeWeights() const { return {}; }
};

// The same, plus the axes.
struct WithAxes : WithoutAxes {
  std::vector<Real> CoLatitudes() const { return {}; }
  std::vector<Real> Longitudes() const { return {}; }
};

static_assert(!GSHTrans::AngularGrid<WithoutAxes>,
              "a grid that cannot hand over its two axes separately is no "
              "longer an AngularGrid");
static_assert(GSHTrans::AngularGrid<WithAxes>,
              "and adding them is the only thing that was missing");
static_assert(
    GSHTrans::AngularGrid<
        GSHTrans::GaussLegendreGrid<double, GSHTrans::All, GSHTrans::All>>);

}  // namespace
