#include <gtest/gtest.h>

#include <GSHTrans/All>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

// Batches of fields at known offsets.
//
// A Batch was a count, a stride and a dist, which is every layout that is
// affine in the field index. Batch::At is the rest: fields that start wherever
// a table says they do -- a subset of radii, storage with padding between
// elements, the same-N components of a tensor in a caller's own array.
//
// The transform reaches a layout only through Count, Stride, Offset, Span and
// Disjoint, and gathers into its own scratch before FFTW or a BLAS sees
// anything. So nothing in a kernel changed to admit this, and the first test
// below is the statement of that: an offset batch describing an affine layout
// must give the affine answer to the last bit, through every kernel.

using namespace GSHTrans;

namespace {

using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

constexpr Int lMax = 12;
constexpr Int n = 1;

auto Coefficients(const Grid& grid, Int count, Int seed) {
  const auto size = grid.CoefficientSize(lMax, n);
  auto values = std::vector<Complex>(static_cast<std::size_t>(count * size));
  auto j = Int{0};
  for (auto& x : values) {
    x = Complex{std::sin(0.37 * static_cast<Real>(j + seed)),
                std::cos(0.11 * static_cast<Real>(j - seed))};
    j++;
  }
  return values;
}

std::vector<Grid> EveryKernel() {
  auto grids =
      std::vector<Grid>{Grid(lMax, n, FFTWpp::Estimate),
                        Grid(lMax, n, FFTWpp::Estimate, Chunking::Fixed(2),
                             WignerValues::Generated())};
#ifdef GSHTRANS_HAVE_BLAS
  grids.push_back(Grid(lMax, n, FFTWpp::Estimate, Chunking::Fixed(3),
                       WignerValues::Stored(), TransformKernel::Matrix()));
#endif
  return grids;
}

}  // namespace

TEST(BatchAt, AnAffineLayoutGivesTheAffineAnswer) {
  constexpr auto count = Int{5};
  for (const auto& grid : EveryKernel()) {
    const auto fieldSize = grid.FieldSize();
    const auto size = grid.CoefficientSize(lMax, n);
    const auto given = Coefficients(grid, count, 3);

    auto fieldStarts = std::vector<Int>{};
    auto blockStarts = std::vector<Int>{};
    for (auto k = Int{0}; k < count; k++) {
      fieldStarts.push_back(k * fieldSize);
      blockStarts.push_back(k * size);
    }

    for (auto policy : {Execution::Sequential(), Execution::Parallel(4)}) {
      auto affineFields = std::vector<Complex>(count * fieldSize);
      auto offsetFields = std::vector<Complex>(count * fieldSize);
      grid.InverseTransformation(lMax, n, given, Batch::Contiguous(count, size),
                                 affineFields,
                                 Batch::Contiguous(count, fieldSize), policy);
      grid.InverseTransformation(lMax, n, given,
                                 Batch::At(blockStarts, 1, size), offsetFields,
                                 Batch::At(fieldStarts, 1, fieldSize), policy);
      ASSERT_EQ(affineFields, offsetFields);

      auto affineBack = std::vector<Complex>(count * size);
      auto offsetBack = std::vector<Complex>(count * size);
      grid.ForwardTransformation(
          lMax, n, affineFields, Batch::Contiguous(count, fieldSize),
          affineBack, Batch::Contiguous(count, size), policy);
      grid.ForwardTransformation(
          lMax, n, offsetFields, Batch::At(fieldStarts, 1, fieldSize),
          offsetBack, Batch::At(blockStarts, 1, size), policy);
      ASSERT_EQ(affineBack, offsetBack);
    }
  }
}

TEST(BatchAt, AShuffledLayoutWithGapsRoundTripsAndLeavesTheGapsAlone) {
  // Three fields, out of order, with padding between and around them: what a
  // finite-element code's storage looks like. The padding carries a sentinel
  // and must still carry it afterwards.
  constexpr auto count = Int{3};
  const auto sentinel = Complex{-12345.0, 54321.0};
  for (const auto& grid : EveryKernel()) {
    const auto fieldSize = grid.FieldSize();
    const auto size = grid.CoefficientSize(lMax, n);
    const auto given = Coefficients(grid, count, 8);

    const auto pad = Int{7};
    const auto slot = fieldSize + pad;
    const auto fieldStarts =
        std::vector<Int>{pad + 2 * slot, pad + 0 * slot, pad + 1 * slot};
    const auto fieldBatch = Batch::At(fieldStarts, 1, fieldSize);
    auto storage = std::vector<Complex>(
        static_cast<std::size_t>(pad + count * slot), sentinel);
    ASSERT_GE(static_cast<Int>(storage.size()), fieldBatch.Span(fieldSize));

    grid.InverseTransformation(lMax, n, given, Batch::Contiguous(count, size),
                               storage, fieldBatch);

    // Every element that belongs to no field is as it was.
    auto owned = std::vector<bool>(storage.size(), false);
    for (auto start : fieldStarts) {
      for (auto j = Int{0}; j < fieldSize; j++) {
        owned[static_cast<std::size_t>(start + j)] = true;
      }
    }
    for (std::size_t i = 0; i < storage.size(); i++) {
      if (!owned[i]) {
        ASSERT_EQ(storage[i], sentinel) << "element " << i;
      }
    }

    auto back = std::vector<Complex>(count * size);
    grid.ForwardTransformation(lMax, n, storage, fieldBatch, back,
                               Batch::Contiguous(count, size));
    for (std::size_t i = 0; i < back.size(); i++) {
      EXPECT_NEAR(std::abs(back[i] - given[i]), 0, 1e-12)
          << "coefficient " << i;
    }
  }
}

TEST(BatchAt, AnInterleavedLayoutCanBeDescribedByItsOffsets) {
  // Stride other than one: fields interleaved element by element, which is
  // Batch::Interleaved when every one takes part and an offset table when
  // only some do.
  auto grid = Grid(lMax, n, FFTWpp::Estimate);
  const auto fieldSize = grid.FieldSize();
  const auto size = grid.CoefficientSize(lMax, n);
  constexpr auto components = Int{5};
  const auto given = Coefficients(grid, 2, 1);

  auto interleaved = std::vector<Complex>(components * fieldSize);
  auto byOffsets = std::vector<Complex>(components * fieldSize);
  // Components 1 and 3 of a point-major five, and nothing else.
  auto fromComponentOne = std::span<Complex>(interleaved).subspan(1);
  grid.InverseTransformation(lMax, n, given, Batch::Contiguous(2, size),
                             fromComponentOne,
                             Batch::Strided(2, components, 2));
  grid.InverseTransformation(lMax, n, given, Batch::Contiguous(2, size),
                             byOffsets,
                             Batch::At({1, 3}, components, fieldSize));
  EXPECT_EQ(interleaved, byOffsets);
}

TEST(BatchSubset, TransformsTheFieldsNamedAndNoOthers) {
  // "Only the solid regions": some radii of a layered field, chosen by index.
  auto grid = Grid(lMax, 0, FFTWpp::Estimate);
  auto radii = std::vector<Real>{};
  for (auto i = 0; i < 7; i++) radii.push_back(0.4 + 0.1 * i);
  const auto radial = RadialGrid<Real>(radii);

  auto field = LayeredSpinField<0, Grid>(radial, grid);
  auto j = Int{0};
  for (auto& x : field.Data()) {
    x = Complex{std::sin(0.01 * static_cast<Real>(j)), 0.5};
    j++;
  }
  const auto all = Expand(field, lMax);

  const auto solid = std::vector<Int>{1, 2, 5};
  const auto size = grid.CoefficientSize(lMax, 0);
  auto some =
      std::vector<Complex>(solid.size() * static_cast<std::size_t>(size));
  grid.ForwardTransformation(
      lMax, 0, field.Data(), field.Batch().Subset(solid, grid.FieldSize()),
      some, Batch::Contiguous(static_cast<Int>(solid.size()), size));

  for (std::size_t k = 0; k < solid.size(); k++) {
    for (auto l = Int{0}; l <= lMax; l++) {
      for (auto m = -l; m <= l; m++) {
        const auto index = GSHIndices<All>(lMax, lMax, 0).Index(l, m);
        EXPECT_EQ(some[k * static_cast<std::size_t>(size) +
                       static_cast<std::size_t>(index)],
                  (all[solid[k], l, m]))
            << "radius " << solid[k];
      }
    }
  }

  // A subset of a subset is a subset, and an index the batch does not have is
  // refused.
  const auto inner = field.Batch()
                         .Subset(solid, grid.FieldSize())
                         .Subset(std::vector<Int>{2}, grid.FieldSize());
  EXPECT_EQ(inner.Count(), 1);
  EXPECT_EQ(inner.Offset(0, 0), 5 * grid.FieldSize());
  EXPECT_THROW(field.Batch().Subset(std::vector<Int>{7}, grid.FieldSize()),
               std::invalid_argument);
  EXPECT_THROW(field.Batch().Subset(std::vector<Int>{-1}, grid.FieldSize()),
               std::invalid_argument);
}

TEST(BatchAt, RefusesLayoutsInWhichFieldsOverlap) {
  // The transform writes every element of every field, so overlapping fields
  // give a wrong answer and not an error. It is checked once, when the batch
  // is made, which is why the size is given then.
  EXPECT_NO_THROW(Batch::At({0, 10, 20}, 1, 10));
  EXPECT_THROW(Batch::At({0, 9, 20}, 1, 10), std::invalid_argument);
  EXPECT_THROW(Batch::At({20, 0, 9}, 1, 10), std::invalid_argument);
  EXPECT_THROW(Batch::At({4, 4}, 1, 1), std::invalid_argument);

  // Interleaved: offsets in different residue classes of the stride never
  // meet, however close; in the same class they need size * stride between.
  EXPECT_NO_THROW(Batch::At({0, 1, 2}, 3, 100));
  EXPECT_NO_THROW(Batch::At({0, 300}, 3, 100));
  EXPECT_THROW(Batch::At({0, 297}, 3, 100), std::invalid_argument);

  EXPECT_THROW(Batch::At({}, 1, 10), std::invalid_argument);
  EXPECT_THROW(Batch::At({-1, 10}, 1, 10), std::invalid_argument);
  EXPECT_THROW(Batch::At({0, 10}, 0, 10), std::invalid_argument);
  EXPECT_THROW(Batch::At({0, 10}, 1, 0), std::invalid_argument);
}

TEST(BatchAt, IsGoodForFieldsNoLongerThanItWasBuiltFor) {
  const auto batch = Batch::At({0, 10, 20}, 1, 10);
  EXPECT_TRUE(batch.Disjoint(10));
  EXPECT_TRUE(batch.Disjoint(4));
  EXPECT_FALSE(batch.Disjoint(11));
  EXPECT_EQ(batch.Span(10), 30);
  EXPECT_EQ(batch.Span(4), 24);

  // And the transform asks.
  auto grid = Grid(lMax, n, FFTWpp::Estimate);
  const auto size = grid.CoefficientSize(lMax, n);
  const auto given = Coefficients(grid, 2, 0);
  auto fields = std::vector<Complex>(2 * grid.FieldSize());
  EXPECT_THROW(
      grid.InverseTransformation(
          lMax, n, given, Batch::Contiguous(2, size), fields,
          Batch::At({0, grid.FieldSize() - 1}, 1, grid.FieldSize() - 1)),
      std::invalid_argument);
}

TEST(BatchAt, OwnsItsOffsetsAndComparesByValue) {
  auto batch = Batch::Contiguous(1, 1);
  {
    auto offsets = std::vector<Int>{30, 0, 15};
    batch = Batch::At(offsets, 1, 10);
    offsets.assign(3, -999);  // the batch must not have been watching
  }
  EXPECT_EQ(batch.Count(), 3);
  EXPECT_TRUE(batch.HasOffsets());
  EXPECT_EQ(batch.Offset(0, 0), 30);
  EXPECT_EQ(batch.Offset(4, 2), 19);

  // Two batches that describe one layout are equal, whichever table they hold
  // and however it was arrived at.
  EXPECT_EQ(batch, Batch::At({30, 0, 15}, 1, 10));
  EXPECT_NE(batch, Batch::At({30, 0, 16}, 1, 10));
  EXPECT_NE(batch, Batch::At({30, 0, 15}, 2, 10));
  EXPECT_FALSE(Batch::Contiguous(3, 10).HasOffsets());
}
