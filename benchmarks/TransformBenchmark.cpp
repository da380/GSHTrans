// The benchmark harness of core-plan.md section 5.
//
// Steps E to H are each supposed to be preceded by the measurement that
// justifies them, and this is that measurement. In particular it checks
// finding P1 -- that the transform is Legendre-stage bandwidth-bound by a wide
// margin -- which is currently an arithmetic estimate and which the whole case
// for batching rests on.
//
// The two stages are separated without instrumenting the transform. The
// coefficient loop is filtered on `l <= lMax`, where lMax is the *call's*
// truncation degree, so calling at the smallest legal degree does the full FFT
// work and almost no Legendre work. The difference between that and a
// full-degree call is the Legendre stage, measured on the production code
// rather than on a replica of it.

#include <GSHTrans/All>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

using namespace GSHTrans;

using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Clock = std::chrono::steady_clock;

namespace {

long ResidentMegabytes() {
  auto file = std::ifstream("/proc/self/status");
  auto key = std::string{};
  while (file >> key) {
    if (key == "VmRSS:") {
      long value = 0;
      file >> value;
      return value / 1024;
    }
    file.ignore(1 << 20, '\n');
  }
  return -1;
}

// Run `action` enough times to measure it, and return seconds per call.
//
// Best of several windows, not the mean of one. This machine's clock scales
// under load: repeating one measurement a few seconds apart moves it by tens
// of percent, which is enough to invent a speedup that is not there. Comparing
// two versions of the code therefore needs them interleaved, and any single
// number here should be read as an upper bound rather than a value.
template <typename Action>
double TimePerCall(Action&& action, double target = 0.15, int windows = 5) {
  action();  // warm up: first call plans, faults pages, fills caches
  auto repetitions = 1;
  auto best = std::numeric_limits<double>::max();
  for (auto window = 0; window < windows; ++window) {
    while (true) {
      const auto start = Clock::now();
      for (auto i = 0; i < repetitions; ++i) action();
      const auto elapsed =
          std::chrono::duration<double>(Clock::now() - start).count();
      if (elapsed >= target || repetitions >= (1 << 20)) {
        best = std::min(best, elapsed / repetitions);
        break;
      }
      repetitions *= 2;
    }
  }
  return best;
}

// A STREAM-style triad, to put the transform's achieved bandwidth on a scale.
double TriadBandwidthGBs() {
  constexpr std::size_t n = 40'000'000;  // ~960 MB touched, far beyond L3
  auto a = std::vector<double>(n, 1.0);
  auto b = std::vector<double>(n, 2.0);
  auto c = std::vector<double>(n, 3.0);
  const auto seconds = TimePerCall([&] {
    for (std::size_t i = 0; i < n; ++i) a[i] = b[i] + 3.0 * c[i];
  });
  const auto bytes = 3.0 * n * sizeof(double);  // two read, one written
  return bytes / seconds / 1e9;
}

// Bytes of Wigner values one transform at this degree and upper index streams:
// the (l, m) block for that upper index, once per colatitude.
double WignerBytes(Int lMax, Int n) {
  const auto perColatitude = GSHIndices<All>(lMax, lMax, n).Size();
  const auto nTheta = lMax + 1;
  return static_cast<double>(nTheta) * static_cast<double>(perColatitude) *
         sizeof(Real);
}

struct Row {
  Int lMax;
  Int n;
  const char* scalar;
  const char* direction;
  double total;
  double stage;  // FFT stage plus call overhead
};

void PrintHeader(const char* title) {
  std::printf("\n%s\n", title);
  for (auto i = 0; i < 78; ++i) std::putchar('-');
  std::putchar('\n');
}

}  // namespace

int main() {
  std::printf("GSHTrans transform benchmark (core-plan.md section 5)\n");
  std::printf("double precision, single field per call (k = 1)\n");

  const auto triad = TriadBandwidthGBs();
  std::printf("\nSTREAM-style triad bandwidth: %.1f GB/s\n", triad);

  //------------------------------------------------------------------------//
  //                       Grid construction and size                        //
  //------------------------------------------------------------------------//

  PrintHeader("Grid construction: time and resident size");
  std::printf("%6s %6s %12s %12s %12s\n", "lMax", "nMax", "build (s)",
              "RSS (MB)", "table (MB)");
  for (auto lMax : {Int{32}, Int{64}, Int{128}, Int{256}}) {
    for (auto nMax : {Int{0}, Int{2}}) {
      const auto before = ResidentMegabytes();
      const auto start = Clock::now();
      auto grid = GaussLegendreGrid<Real, All, All>(lMax, nMax,
                                                    FFTWpp::Estimate);
      const auto seconds =
          std::chrono::duration<double>(Clock::now() - start).count();
      const auto after = ResidentMegabytes();

      auto tableBytes = 0.0;
      for (auto n = -nMax; n <= nMax; ++n) tableBytes += WignerBytes(lMax, n);
      std::printf("%6zd %6zd %12.3f %12ld %12.1f\n", lMax, nMax, seconds,
                  after - before, tableBytes / 1e6);
    }
  }

  //------------------------------------------------------------------------//
  //                        Transforms, and the split                        //
  //------------------------------------------------------------------------//

  PrintHeader("Transforms: total, stage split, and Legendre bandwidth");
  std::printf("%6s %4s %8s %9s %10s %10s %9s %9s\n", "lMax", "n", "scalar",
              "direction", "total(ms)", "FFT(ms)", "Leg(ms)", "GB/s");

  for (auto lMax : {Int{32}, Int{64}, Int{128}, Int{256}}) {
    for (auto n : {Int{0}, Int{2}}) {
      auto grid = GaussLegendreGrid<Real, All, All>(lMax, n, FFTWpp::Measure);
      const auto stub = std::abs(n);  // smallest legal call degree

      const auto fieldSize = grid.FieldSize();
      const auto full = grid.CoefficientSize(lMax, n);
      const auto small = grid.CoefficientSize(stub, n);

      auto complexField = FFTWpp::vector<Complex>(fieldSize);
      auto realField = FFTWpp::vector<Real>(fieldSize);
      for (auto i = Int{0}; i < fieldSize; ++i) {
        complexField[i] = Complex{0.5 + 0.001 * i, -0.25 + 0.002 * i};
        realField[i] = 0.5 + 0.001 * i;
      }
      auto fullCoefficients = FFTWpp::vector<Complex>(full);
      auto smallCoefficients = FFTWpp::vector<Complex>(small);
      auto realFull = FFTWpp::vector<Complex>(grid.RealCoefficientSize(lMax));
      auto realSmall = FFTWpp::vector<Complex>(grid.RealCoefficientSize(stub));

      const auto bytes = WignerBytes(lMax, n);

      auto report = [&](const char* scalar, const char* direction,
                        double total, double stage) {
        const auto legendre = total - stage;
        const auto gbs = legendre > 0 ? bytes / legendre / 1e9 : 0.0;
        std::printf("%6zd %4zd %8s %9s %10.3f %10.3f %9.3f %9.1f\n", lMax, n,
                    scalar, direction, total * 1e3, stage * 1e3,
                    legendre * 1e3, gbs);
      };

      report("complex", "forward",
             TimePerCall([&] {
               grid.ForwardTransformation(lMax, n, complexField,
                                          fullCoefficients);
             }),
             TimePerCall([&] {
               grid.ForwardTransformation(stub, n, complexField,
                                          smallCoefficients);
             }));

      report("complex", "inverse",
             TimePerCall([&] {
               grid.InverseTransformation(lMax, n, fullCoefficients,
                                          complexField);
             }),
             TimePerCall([&] {
               grid.InverseTransformation(stub, n, smallCoefficients,
                                          complexField);
             }));

      // Real-valued fields exist only at upper index zero (core step A).
      if (n == 0) {
        report("real", "forward",
               TimePerCall([&] {
                 grid.ForwardTransformation(lMax, 0, realField, realFull);
               }),
               TimePerCall([&] {
                 grid.ForwardTransformation(stub, 0, realField, realSmall);
               }));
        report("real", "inverse",
               TimePerCall([&] {
                 grid.InverseTransformation(lMax, 0, realFull, realField);
               }),
               TimePerCall([&] {
                 grid.InverseTransformation(stub, 0, realSmall, realField);
               }));
      }
    }
  }

  //------------------------------------------------------------------------//
  //                              Threading                                  //
  //------------------------------------------------------------------------//

  PrintHeader("Threading (core-plan.md step H)");
  std::printf("%6s %4s %9s %8s %10s %9s %9s\n", "lMax", "n", "direction",
              "threads", "time(ms)", "speedup", "GB/s");
  for (auto lMax : {Int{128}, Int{256}}) {
    const auto n = Int{2};
    auto grid = GaussLegendreGrid<Real, All, All>(lMax, n, FFTWpp::Measure);
    auto field = FFTWpp::vector<Complex>(grid.FieldSize());
    auto coefficients = FFTWpp::vector<Complex>(grid.CoefficientSize(lMax, n));
    for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
      field[i] = Complex{0.5 + 0.001 * i, -0.25 + 0.002 * i};
    }
    const auto bytes = WignerBytes(lMax, n);

    for (const char* direction : {"forward", "inverse"}) {
      auto base = 0.0;
      for (auto threads : {1, 2, 4, 8, 16}) {
        const auto policy = threads == 1 ? Execution::Sequential()
                                         : Execution::Parallel(threads);
        const auto seconds = TimePerCall([&] {
          if (direction[0] == 'f') {
            grid.ForwardTransformation(lMax, n, field, coefficients, policy);
          } else {
            grid.InverseTransformation(lMax, n, coefficients, field, policy);
          }
        });
        if (threads == 1) base = seconds;
        std::printf("%6zd %4zd %9s %8d %10.3f %8.2fx %9.1f\n", lMax, n,
                    direction, threads, seconds * 1e3, base / seconds,
                    bytes / seconds / 1e9);
      }
    }
  }

  std::printf(
      "\nFFT column is a call truncated to the smallest legal degree: full FFT\n"
      "work, negligible Legendre work. Leg is the difference. GB/s is the\n"
      "Wigner bytes that call streams divided by the Legendre time.\n");

  FFTWpp::CleanUp();
}
