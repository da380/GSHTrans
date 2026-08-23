// 07 -- Batches and threads
//
// The transform's primitive is "transform k fields sharing a grid, a degree
// and an upper index", with the single field as k = 1. Batching amortises the
// Wigner values, which is where nearly all of the cost is.
//
// Batched is also the regime where the choice of Legendre kernel matters most,
// and there are two -- see example 20.

#include <GSHTrans/All>
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
  auto grid = Grid(lMax, n, FFTWpp::Measure, Chunking::ForCache(Int{16} << 20));

  const auto fieldSize = static_cast<Int>(grid.FieldSize());
  const auto coefficientSize = static_cast<Int>(grid.CoefficientSize(lMax, n));

  auto fields = std::vector<Complex>(count * fieldSize);
  for (Int i = 0; i < count * fieldSize; i++) {
    fields[i] = Complex{std::cos(0.01 * i), std::sin(0.02 * i)};
  }
  auto coefficients = std::vector<Complex>(count * coefficientSize);

  // A batch is described by (count, stride, dist) rather than by contiguity,
  // so it covers fields laid end to end and fields interleaved with others
  // alike. Here they are contiguous: stride 1, dist FieldSize.
  const auto in = Batch::Contiguous(count, fieldSize);
  const auto out = Batch::Contiguous(count, coefficientSize);

  const auto time = [](auto&& action) {
    const auto start = std::chrono::steady_clock::now();
    action();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         start)
        .count();
  };

  // One field at a time.
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
  const auto batched =
      time([&] { grid.ForwardTransformation(lMax, n, fields, in,
                                            coefficients, out); });

  // Threading is an explicit per-call policy and defaults to sequential: the
  // library never creates threads because it can, only because it was asked.
  // Exactly one level threads, so a caller already inside a parallel region
  // gets a sequential transform rather than nested teams.
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
  // across components -- which example 08 shows the field layer arranging.
}
