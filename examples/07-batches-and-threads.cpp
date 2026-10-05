// 07 -- Batches and threads
//
// The transform's primitive is "transform k fields sharing a grid, a degree
// and an upper index", with the single field as k = 1. Batching amortises the
// Wigner values, which is where nearly all of the cost is.
//
// Batched is also the regime where the choice of Legendre kernel matters most,
// and there are two -- see example 20.
//
// What this shows
//   The grid's own transform on raw buffers, one field at a time and as a
//   batch; the batch descriptor; the chunking policy set at construction;
//   and the explicit, per-call threading policy.
//
// Assumes
//   Example 06 (what the transform computes), example 05 (raw buffers in
//   the library's sample order).
//
// Introduced
//   grid.ForwardTransformation(lMax, n, in, [inBatch,] out, [outBatch,]
//   [policy]), grid.CoefficientSize, Batch::Contiguous, Chunking::ForCache,
//   FFTWpp::Measure, Execution::Parallel.
//
// Output
//   Wall-clock times for eight transforms done separately, as one batch, and
//   as one batch on four threads. The numbers depend on the machine; only
//   their order is the point.
//
// See docs/gshtrans-reference.tex, section "Batching, chunking and
// threading".

#include <GSHTrans/GSHTrans.hpp>
#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  using Int = std::ptrdiff_t;

  constexpr auto lMax = Int{64};
  constexpr auto n = Int{2};
  constexpr auto count = Int{8};

  // The chunking policy is a property of the machine, not of the call, so it
  // is given once. Automatic assumes a modest cache; ForCache takes the real
  // figure and does better.
  //
  // A chunk is how many fields of a batch the inner loop takes at once. The
  // gain from batching saturates and then reverses once the chunk's
  // coefficients stop fitting in last-level cache, so a large batch is
  // processed in chunks sized from the cache. ForCache takes the machine's
  // total last-level cache in bytes (16 MiB here); Automatic assumes 8 MiB.
  //
  // The third argument is FFTW's planner flag, how hard FFTW works to plan
  // each transform shape. Plans are made on first use, once per thread and
  // shape, and kept.
  auto grid = Grid(lMax, n, FFTWpp::Measure, Chunking::ForCache(Int{16} << 20));

  // The sizes of one field's samples and one field's coefficients at this
  // degree and upper index.
  const auto fieldSize = grid.FieldSize();
  const auto coefficientSize = static_cast<Int>(grid.CoefficientSize(lMax, n));

  auto fields = std::vector<Complex>(count * fieldSize);
  for (Int i = 0; i < count * fieldSize; i++) {
    fields[i] = Complex{std::cos(0.01 * i), std::sin(0.02 * i)};
  }
  auto coefficients = std::vector<Complex>(count * coefficientSize);

  // A batch is described by (count, stride, dist) rather than by contiguity,
  // so it covers fields laid end to end and fields interleaved with others
  // alike: element j of field k lives at j * stride + k * dist. Here they are
  // contiguous: stride 1, dist FieldSize. Fields and coefficients take
  // separate descriptors, since their sizes, and so their dist, differ.
  // Batch::Interleaved and Batch::Strided name the other affine layouts.
  const auto in = Batch::Contiguous(count, fieldSize);
  const auto out = Batch::Contiguous(count, coefficientSize);

  const auto time = [](auto&& action) {
    const auto start = std::chrono::steady_clock::now();
    action();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         start)
        .count();
  };

  // Each timing below begins with the first call of its kind, so it also
  // pays for FFTW planning of that shape on each thread involved. That is
  // part of what a first call costs; a careful benchmark would warm up first.

  // One field at a time, through the k = 1 overload, which takes one field's
  // samples and one field's coefficients and no descriptors.
  const auto separate = time([&] {
    for (Int k = 0; k < count; k++) {
      auto one = std::span(fields).subspan(k * fieldSize, fieldSize);
      auto block =
          std::span(coefficients).subspan(k * coefficientSize, coefficientSize);
      grid.ForwardTransformation(lMax, n, one, block);
    }
  });

  // All of them together. The Wigner block for this upper index is streamed
  // once for the batch rather than once per field, which is the whole of what
  // batching buys.
  const auto batched = time([&] {
    grid.ForwardTransformation(lMax, n, fields, in, coefficients, out);
  });

  // Threading is an explicit per-call policy and defaults to sequential: the
  // library never creates threads because it can, only because it was asked.
  // Exactly one level threads, so a caller already inside a parallel region
  // gets a sequential transform rather than nested teams.
  //
  // Parallel(4) asks for four threads; Parallel() leaves the count to OpenMP.
  // In a build without OpenMP the same call runs on one thread. On a machine
  // with simultaneous multithreading, ask for cores rather than hardware
  // threads: this work is memory-bound.
  const auto threaded = time([&] {
    grid.ForwardTransformation(lMax, n, fields, in, coefficients, out,
                               Execution::Parallel(4));
  });

  std::cout << "eight fields, one at a time  " << separate * 1e3 << " ms\n"
            << "as one batch                 " << batched * 1e3 << " ms\n"
            << "batched on four threads      " << threaded * 1e3 << " ms\n";

  // A batch shares grid, degree *and upper index*: the Wigner block is what
  // is being amortised, so fields at different upper indices cannot batch
  // together. For a tensor that means batching over radii within each n, not
  // across components -- which is why the tensor fields of example 08 order
  // their buffer by upper index, making each upper index's components one
  // contiguous batch.
}
