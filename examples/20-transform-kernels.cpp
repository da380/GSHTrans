// 20 -- Two Legendre kernels, and choosing between them
//
// The transform's inner stage can be arranged two ways, and the library keeps
// both rather than picking one.
//
//   Loop   -- a colatitude at a time, an axpy per (l, m). What the library has
//             always done, and the default.
//   Matrix -- every FFT first, then one matrix product per order against a
//             table laid out so that each order's block is contiguous. Needs a
//             BLAS.
//
// They compute the same sums in different orders, so they agree to rounding
// rather than exactly. Two reasons both are kept:
//
//   -- Which one suits a machine is a question that machine can answer, and
//      only by being asked. See benchmarks/TransformBenchmark kernels.
//   -- Each is the other's check. A GEMM sums in whatever order its kernel
//      chooses, so the matrix path cannot be verified against itself; it is
//      verified against the loop path on identical inputs.
//
// See core-plan.md section 11.

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

#ifndef GSHTRANS_HAVE_BLAS
  std::cout << "Built without a BLAS, so TransformKernel::Matrix() does not\n"
               "exist and the loop kernel is the only one. Everything the\n"
               "library does is available; nothing is withdrawn.\n"
               "Configure with -DGSHTRANS_WITH_BLAS=ON to require one.\n";
  return 0;
#else

  // The kernel is a property of the grid and not of the call, because the two
  // want different Wigner layouts and holding both would double the table.
  // It comes after the chunking and the Wigner-values policy, and defaults to
  // Loop, so no existing spelling has to move.
  auto loop = Grid(lMax, n, FFTWpp::Measure);
  auto matrix = Grid(lMax, n, FFTWpp::Measure, Chunking::Automatic(),
                     WignerValues::Stored(), TransformKernel::Matrix());

  const auto fieldSize = static_cast<Int>(loop.FieldSize());
  const auto coefficientSize = static_cast<Int>(loop.CoefficientSize(lMax, n));

  auto fields = std::vector<Complex>(count * fieldSize);
  for (Int i = 0; i < count * fieldSize; i++) {
    fields[i] = Complex{std::cos(0.013 * i), std::sin(0.021 * i)};
  }
  const auto fieldBatch = Batch::Contiguous(count, fieldSize);
  const auto coefficientBatch = Batch::Contiguous(count, coefficientSize);

  auto fromLoop = std::vector<Complex>(count * coefficientSize);
  auto fromMatrix = std::vector<Complex>(count * coefficientSize);

  const auto policy = Execution::Parallel();
  loop.ForwardTransformation(lMax, n, fields, fieldBatch, fromLoop,
                             coefficientBatch, policy);
  matrix.ForwardTransformation(lMax, n, fields, fieldBatch, fromMatrix,
                               coefficientBatch, policy);

  // The agreement, checked rather than asserted in prose. Relative to the
  // largest coefficient, because the absolute size depends on the data.
  auto largest = Real{0};
  auto worst = Real{0};
  for (std::size_t i = 0; i < fromLoop.size(); i++) {
    largest = std::max(largest, std::abs(fromLoop[i]));
    worst = std::max(worst, std::abs(fromMatrix[i] - fromLoop[i]));
  }
  std::cout << "kernels agree to " << worst / largest << " relative\n";
  if (!(worst < 1e-12 * largest)) {
    std::cerr << "kernels disagree by more than rounding\n";
    return 1;
  }

  // What it is worth. One size on one machine is an anecdote, not a
  // measurement -- the benchmark's `kernels` section is the measurement, and
  // it reports both regimes separately because the answer differs between
  // them. Batched, which is what the field and layered layers always do, the
  // matrix kernel wins; unbatched at high degree on many threads the loop
  // kernel is already at the memory roof and there is nothing to win.
  auto TimeIt = [&](auto& grid) {
    using Clock = std::chrono::steady_clock;
    auto scratch = std::vector<Complex>(count * coefficientSize);
    grid.ForwardTransformation(lMax, n, fields, fieldBatch, scratch,
                               coefficientBatch, policy);  // warm up
    const auto start = Clock::now();
    for (int rep = 0; rep < 5; rep++) {
      grid.ForwardTransformation(lMax, n, fields, fieldBatch, scratch,
                                 coefficientBatch, policy);
    }
    return std::chrono::duration<double, std::milli>(Clock::now() - start)
               .count() /
           5;
  };
  const auto tLoop = TimeIt(loop);
  const auto tMatrix = TimeIt(matrix);
  std::cout << "forward, lMax = " << lMax << ", k = " << count << ": loop "
            << tLoop << " ms, matrix " << tMatrix << " ms\n";
#ifndef NDEBUG
  // Worth saying rather than letting the reader believe the ratio above.
  // `cmake -S . -B build` leaves CMAKE_BUILD_TYPE empty, and an unoptimised
  // build runs about ten times slow with every figure self-consistent, so
  // nothing looks wrong -- it flatters whichever kernel does less work in the
  // interpreter-like regime rather than in the real one. Configure with
  // -DCMAKE_BUILD_TYPE=Release before believing any of it.
  std::cout << "  (unoptimised build: those two numbers mean nothing --\n"
               "   configure with -DCMAKE_BUILD_TYPE=Release)\n";
#endif

  // The inverse uses the same stored matrices, transposed, so there is no
  // second table and no second decision.
  auto back = std::vector<Complex>(count * fieldSize);
  matrix.InverseTransformation(lMax, n, fromMatrix, coefficientBatch, back,
                               fieldBatch, policy);

  // Two combinations the matrix kernel refuses, and both are facts rather
  // than gaps.
  //
  // Generated Wigner values run the recursion inside the transform, and the
  // recursion produces every order of one colatitude together -- so a single
  // order's block cannot be had without either keeping the whole table, which
  // is what generating exists to avoid, or repeating the recursion for every
  // order.
  try {
    auto bad = Grid(lMax, n, FFTWpp::Measure, Chunking::Automatic(),
                    WignerValues::Generated(), TransformKernel::Matrix());
    std::cerr << "expected the generated combination to be refused\n";
    return 1;
  } catch (const std::invalid_argument&) {
    std::cout << "matrix kernel refuses generated Wigner values, as it must\n";
  }

  // And BLAS offers single and double precision and nothing wider, so a grid
  // over long double keeps the loop kernel and only the loop kernel.
  try {
    auto wide = GaussLegendreGrid<long double, All, All>(
        lMax, n, FFTWpp::Measure, Chunking::Automatic(),
        WignerValues::Stored(), TransformKernel::Matrix());
    std::cerr << "expected long double to be refused\n";
    return 1;
  } catch (const std::invalid_argument&) {
    std::cout << "matrix kernel refuses long double, since BLAS has none\n";
  }

  // One obligation, which no caller can meet without being told. The matrix
  // kernel threads over orders and issues every product from inside an OpenMP
  // region, so a BLAS built on the same OpenMP runtime is nested and runs
  // serial with nothing to configure. A BLAS with a thread pool of its own
  // cannot see that region, and must be held to one thread by the environment
  // -- OPENBLAS_NUM_THREADS=1 or the equivalent. These products are skinny,
  // and threading them loses even when there is no nesting to worry about.
  std::cout << "\nremember: a pthread-pool BLAS must be held to one thread\n";

  return 0;
#endif
}
