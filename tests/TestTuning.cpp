#include <gtest/gtest.h>

#include <GSHTrans/All>

#include <complex>
#include <cstddef>
#include <vector>

namespace {

using namespace GSHTrans;

using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

// Small, and planned with Estimate: this file times transforms, and the point
// of every test below is a structural property of the choice rather than a
// number, so the sizes are the smallest that exercise a batch.
auto TestGrid(Int lMax, Int nMax) {
  return Grid(lMax, nMax, FFTWpp::Estimate);
}

//--------------------------------------------------------------------------//
//                     W2: the timing core's discipline                      //
//--------------------------------------------------------------------------//

TEST(Tuning, BestSecondsTakesTheBestWindowAndWarmsUpFirst) {
  auto calls = 0;
  const auto seconds = TuningDetails::BestSeconds([&] { calls++; }, 4);

  // One warm-up plus one per window. The warm-up is not optional: the first
  // call plans, faults its pages and fills its caches, none of which is what
  // is being compared.
  EXPECT_EQ(calls, 5);
  EXPECT_GE(seconds, 0.0);
}

//--------------------------------------------------------------------------//
//                     W3: what tuning the chunk promises                    //
//--------------------------------------------------------------------------//

// The acceptance criterion of core-plan.md section 12.5, and it is
// deliberately not "it finds the optimum". There may not be one resolvable --
// section 8 records the lMax = 128 peak swapping between two runs of the same
// binary -- so the property a caller needs, and the only one that is
// testable, is that tuning never returns something worse than doing nothing.
//
// It holds by construction, since a candidate displaces the incumbent only by
// beating it by the margin, and these assertions are what say the
// construction is what it claims.
TEST(Tuning, NeverChoosesWorseThanTheDefault) {
  for (auto lMax : {Int{8}, Int{16}}) {
    for (auto count : {Int{1}, Int{4}}) {
      auto grid = TestGrid(lMax, 2);
      const auto tuned = TuneChunking(grid, lMax, 2, count);

      EXPECT_GT(tuned.defaultSeconds, 0.0) << "lMax " << lMax;
      EXPECT_GT(tuned.seconds, 0.0) << "lMax " << lMax;

      // Either the incumbent held, or a candidate beat it by the margin.
      if (tuned.conclusive) {
        EXPECT_LT(tuned.seconds, tuned.defaultSeconds * (1 - TuningMargin))
            << "a conclusive result must have beaten the margin";
        EXPECT_GT(tuned.Speedup(), 1.0);
      } else {
        EXPECT_EQ(tuned.chunking, grid.ChunkingPolicy())
            << "an inconclusive result must return the incumbent";
        EXPECT_DOUBLE_EQ(tuned.seconds, tuned.defaultSeconds);
      }
    }
  }
}

// A tuned grid is the same grid: tuning changes a schedule and must change no
// answer. This is the property that makes substituting one for the other safe,
// and it is asserted through the tuner rather than only through With().
TEST(Tuning, ATunedGridAnswersIdentically) {
  constexpr auto lMax = Int{8};
  constexpr auto n = Int{2};
  constexpr auto count = Int{4};

  auto grid = TestGrid(lMax, n);
  const auto fieldSize = static_cast<Int>(grid.FieldSize());
  const auto coefficientSize =
      static_cast<Int>(grid.CoefficientSize(lMax, n));

  auto fields = FFTWpp::vector<Complex>(
      static_cast<std::size_t>(count * fieldSize));
  for (auto j = std::size_t{0}; j < fields.size(); ++j) {
    fields[j] = Complex(static_cast<Real>(j % 17) / 17,
                        static_cast<Real>(j % 23) / 23);
  }

  const auto Forward = [&](const Grid& g) {
    auto out = FFTWpp::vector<Complex>(
        static_cast<std::size_t>(count * coefficientSize));
    g.ForwardTransformation(lMax, n, fields, Batch::Contiguous(count, fieldSize),
                            out, Batch::Contiguous(count, coefficientSize));
    return out;
  };

  const auto tuned = TuneChunking(grid, lMax, n, count);
  const auto reference = Forward(grid);
  const auto got = Forward(grid.With(tuned.chunking));

  ASSERT_EQ(got.size(), reference.size());
  for (auto j = std::size_t{0}; j < got.size(); ++j) {
    EXPECT_EQ(got[j], reference[j]) << "coefficient " << j;
  }
}

// [C19]: Tune hands back values rather than a configured grid, so nothing is
// substituted behind the caller's back. The grid it was asked about is
// untouched, and using the answer is the caller's own step.
TEST(Tuning, LeavesTheGridItWasAskedAboutAlone) {
  auto grid = TestGrid(8, 2);
  const auto before = grid.ChunkingPolicy();
  const auto tuned = TuneChunking(grid, 8, 2, 4);

  EXPECT_EQ(grid.ChunkingPolicy(), before)
      << "tuning must not configure the grid it measured";
  EXPECT_EQ(grid.With(tuned.chunking).ChunkingPolicy(), tuned.chunking);
  EXPECT_EQ(grid.With(tuned.chunking).Identity(), grid.Identity());
}

// At one field there is nothing to tune: every chunk of one or more takes the
// whole batch in one go, so all the candidates run the same schedule and are
// collapsed to one. Pinned because the first version of the tuner did not
// collapse them, timed seven identical experiments, and duly reported a 2x
// "win" -- which at count == 1 is definitionally impossible and is what said
// the method rather than the result was wrong.
TEST(Tuning, HasNothingToTuneAtOneField) {
  auto grid = TestGrid(16, 2);
  const auto tuned = TuneChunking(grid, 16, 2, 1);

  EXPECT_EQ(tuned.candidates, 1) << "one field admits one schedule";
  EXPECT_FALSE(tuned.conclusive);
  EXPECT_EQ(tuned.chunking, grid.ChunkingPolicy());
  EXPECT_DOUBLE_EQ(tuned.Speedup(), 1.0);
}

// And where there is something to tune, the candidates are still deduplicated
// rather than all seven being timed regardless.
TEST(Tuning, CollapsesCandidatesThatWouldRunTheSameSchedule) {
  auto grid = TestGrid(16, 2);
  const auto tuned = TuneChunking(grid, 16, 2, 8);

  EXPECT_GE(tuned.candidates, 1);
  EXPECT_LE(tuned.candidates,
            static_cast<int>(TuningCacheCandidates().size()) + 1);
}

TEST(Tuning, RefusesAnEmptyBatch) {
  auto grid = TestGrid(8, 2);
  EXPECT_THROW(TuneChunking(grid, 8, 2, 0), std::invalid_argument);
}

// Every candidate is a cache figure rather than a chunk, which is what lets
// one answer serve every batch size (section 12.4).
TEST(Tuning, SweepsCacheFiguresAndTheyAreDistinct) {
  const auto candidates = TuningCacheCandidates();
  ASSERT_GE(candidates.size(), 4u);
  for (std::size_t i = 1; i < candidates.size(); ++i) {
    EXPECT_GT(candidates[i], candidates[i - 1]);
  }
  EXPECT_GE(candidates.front(), Int{1} << 20);
}

}  // namespace
