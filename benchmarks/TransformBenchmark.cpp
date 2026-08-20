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
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <utility>
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

//--------------------------------------------------------------------------//
//                             Machine facts                                //
//--------------------------------------------------------------------------//
//
// Printed at the top of every run. The numbers below are unreadable without
// them -- step H already found that cores and hardware threads are different
// answers, and on a multi-socket machine the node count is a third. All of
// this is Linux sysfs; elsewhere the fields come back unknown and the
// benchmark still runs.

std::string FirstLine(const std::string& path) {
  auto file = std::ifstream(path);
  auto line = std::string{};
  std::getline(file, line);
  return line;
}

// "0-7,64-71" -> 16.
int CountCpuList(const std::string& text) {
  auto count = 0;
  auto i = std::size_t{0};
  while (i < text.size()) {
    const auto comma = text.find(',', i);
    const auto part =
        text.substr(i, comma == std::string::npos ? comma : comma - i);
    const auto dash = part.find('-');
    if (dash == std::string::npos) {
      count += 1;
    } else {
      count += std::atoi(part.c_str() + dash + 1) - std::atoi(part.c_str()) + 1;
    }
    if (comma == std::string::npos) break;
    i = comma + 1;
  }
  return count;
}

// Distinct (socket, core) pairs in /proc/cpuinfo. Zero if it does not carry
// them, in which case the caller must say so rather than quietly reporting
// hardware threads as cores.
int PhysicalCores() {
  auto file = std::ifstream("/proc/cpuinfo");
  auto line = std::string{};
  auto seen = std::set<std::pair<int, int>>{};
  auto socket = 0;
  while (std::getline(file, line)) {
    const auto colon = line.find(':');
    if (colon == std::string::npos) continue;
    auto key = line.substr(0, colon);
    while (!key.empty() && (key.back() == ' ' || key.back() == '\t')) {
      key.pop_back();
    }
    const auto value = std::atoi(line.c_str() + colon + 1);
    if (key == "physical id") {
      socket = value;
    } else if (key == "core id") {
      seen.insert({socket, value});
    }
  }
  return static_cast<int>(seen.size());
}

// Size in kilobytes of the cache at this level on cpu0, and how many hardware
// threads share it. The second is what step F's chunk formula needs: P8 sizes
// a chunk against per-core last-level cache, not against the whole of it.
long CacheKilobytes(int level, int& sharedBy) {
  sharedBy = 0;
  for (auto index = 0; index < 10; ++index) {
    const auto dir = "/sys/devices/system/cpu/cpu0/cache/index" +
                     std::to_string(index) + "/";
    const auto levelText = FirstLine(dir + "level");
    if (levelText.empty() || std::atoi(levelText.c_str()) != level) continue;
    if (FirstLine(dir + "type") == "Instruction") continue;
    const auto sizeText = FirstLine(dir + "size");
    if (sizeText.empty()) continue;
    auto kilobytes = std::atol(sizeText.c_str());
    if (sizeText.find('M') != std::string::npos) kilobytes *= 1024;
    sharedBy = CountCpuList(FirstLine(dir + "shared_cpu_list"));
    return kilobytes;
  }
  return 0;
}

int NumaNodes() {
  auto count = 0;
  for (auto node = 0; node < 512; ++node) {
    const auto path = "/sys/devices/system/node/node" + std::to_string(node) +
                      "/cpulist";
    if (!FirstLine(path).empty()) ++count;
  }
  return count;
}

long MemAvailableMegabytes() {
  auto file = std::ifstream("/proc/meminfo");
  auto key = std::string{};
  while (file >> key) {
    if (key == "MemAvailable:") {
      long value = 0;
      file >> value;
      return value / 1024;
    }
    file.ignore(1 << 20, '\n');
  }
  return -1;
}

const char* Environment(const char* name) {
  const auto* value = std::getenv(name);
  return value != nullptr ? value : "(unset)";
}

int HardwareThreads() {
  const auto reported = static_cast<int>(std::thread::hardware_concurrency());
  return reported > 0 ? reported : omp_get_max_threads();
}

void PrintMachineFacts() {
  const auto threads = HardwareThreads();
  const auto cores = PhysicalCores();
  auto l3Shared = 0;
  auto l2Shared = 0;
  const auto l3 = CacheKilobytes(3, l3Shared);
  const auto l2 = CacheKilobytes(2, l2Shared);

  std::printf("\nMachine\n");
  for (auto i = 0; i < 78; ++i) std::putchar('-');
  std::putchar('\n');
  std::printf("  model            %s\n",
              FirstLine("/sys/devices/virtual/dmi/id/product_name").c_str());
  std::printf("  hardware threads %d\n", threads);
  if (cores > 0) {
    std::printf("  physical cores   %d (%.0f threads per core)\n", cores,
                static_cast<double>(threads) / cores);
  } else {
    std::printf("  physical cores   unknown (/proc/cpuinfo carries no core id)\n");
  }
  std::printf("  NUMA nodes       %d\n", NumaNodes());
  if (l2 > 0) std::printf("  L2               %ld KB, shared by %d threads\n", l2, l2Shared);
  if (l3 > 0) {
    std::printf("  L3               %ld KB, shared by %d threads (%.1f MB per core)\n",
                l3, l3Shared, l3Shared > 0 && cores > 0
                    ? l3 / 1024.0 / (l3Shared / (static_cast<double>(threads) / cores))
                    : 0.0);
  }
  std::printf("  MemAvailable     %ld MB\n", MemAvailableMegabytes());
  std::printf("  omp_get_max_threads %d\n", omp_get_max_threads());
  std::printf("  OMP_NUM_THREADS=%s  OMP_PROC_BIND=%s  OMP_PLACES=%s\n",
              Environment("OMP_NUM_THREADS"), Environment("OMP_PROC_BIND"),
              Environment("OMP_PLACES"));
}

// Powers of two up to the hardware thread count, with the physical core count
// inserted -- step H found that the last doubling, from cores to threads, goes
// the wrong way, so the ladder has to contain both to show it.
std::vector<int> ThreadLadder() {
  const auto threads = HardwareThreads();
  const auto cores = PhysicalCores();
  auto ladder = std::set<int>{1};
  for (auto t = 2; t <= threads; t *= 2) ladder.insert(t);
  if (cores > 1 && cores <= threads) ladder.insert(cores);
  ladder.insert(threads);
  return std::vector<int>(ladder.begin(), ladder.end());
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

//--------------------------------------------------------------------------//
//                          Bandwidth roofs                                 //
//--------------------------------------------------------------------------//
//
// Two of them, because the transform's GB/s column needs a ceiling and the
// ceiling is not one number on a multi-socket machine.
//
// The `touch` argument is the point. Pages are placed on the NUMA node of the
// thread that first writes them, and `Wigner::_data` is a std::vector<Real>
// built by its size constructor (Wigner.h:131), so the whole table is
// zero-filled by the single constructing thread and lives on one node however
// many nodes the machine has. Touching with one thread reproduces that;
// touching with the full team is the roof the transform could reach if the
// table were placed deliberately. On a one-node machine the two agree, and
// that agreement is itself the finding.

// Uninitialised, then written by `touch` threads: `new double[n]` does not
// touch, so first touch really is the parallel loop.
std::unique_ptr<double[]> MakeTouched(std::size_t n, int touch, double value) {
  auto data = std::unique_ptr<double[]>(new double[n]);
  auto* raw = data.get();
  const auto count = static_cast<std::ptrdiff_t>(n);
#pragma omp parallel for schedule(static) num_threads(touch)
  for (std::ptrdiff_t i = 0; i < count; ++i) raw[i] = value;
  return data;
}

// STREAM triad: two streams read, one written.
double TriadBandwidthGBs(int threads, int touch, int windows = 5) {
  constexpr std::size_t n = 40'000'000;  // ~960 MB touched, far beyond any L3
  auto a = MakeTouched(n, touch, 1.0);
  auto b = MakeTouched(n, touch, 2.0);
  auto c = MakeTouched(n, touch, 3.0);
  auto* ap = a.get();
  auto* bp = b.get();
  auto* cp = c.get();
  const auto count = static_cast<std::ptrdiff_t>(n);
  const auto seconds = TimePerCall(
      [&] {
#pragma omp parallel for schedule(static) num_threads(threads)
        for (std::ptrdiff_t i = 0; i < count; ++i) ap[i] = bp[i] + 3.0 * cp[i];
      },
      0.15, windows);
  return 3.0 * n * sizeof(double) / seconds / 1e9;
}

// A read-only scan over an array of the given size. This is the closer
// analogue of the Legendre stage, which reads the Wigner table and writes
// something much smaller, so it is the roof the GB/s column should be read
// against.
double ScanBandwidthGBs(double bytes, int threads, int touch, int windows = 5) {
  const auto n = static_cast<std::size_t>(bytes / sizeof(double));
  auto a = MakeTouched(n, touch, 1.0);
  auto* ap = a.get();
  const auto count = static_cast<std::ptrdiff_t>(n);
  static volatile double sink = 0.0;  // keeps the sum from being optimised out
  const auto seconds = TimePerCall(
      [&] {
        auto sum = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : sum) \
    num_threads(threads)
        for (std::ptrdiff_t i = 0; i < count; ++i) sum += ap[i];
        sink = sum;
      },
      0.15, windows);
  return static_cast<double>(n) * sizeof(double) / seconds / 1e9;
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

// Which sections to run. Named on the command line, all of them if none is
// named. A single section takes a fraction of the whole, which is what makes
// an A/B of two builds affordable -- and an A/B run back to back on an idle
// machine is the only kind worth having here: the same binary measured 1.75 ms
// and 2.23 ms on this laptop in different power states, so figures from
// different sessions cannot be compared at all.
std::vector<std::string> sectionsWanted;

bool Want(const std::string& name) {
  return sectionsWanted.empty() ||
         std::find(sectionsWanted.begin(), sectionsWanted.end(), name) !=
             sectionsWanted.end();
}

// For sections that cost too much to run by accident: named or not run.
bool WantNamed(const std::string& name) {
  return !sectionsWanted.empty() && Want(name);
}

double TableMegabytes(Int lMax, Int nMax) {
  auto bytes = 0.0;
  for (auto n = -nMax; n <= nMax; ++n) bytes += WignerBytes(lMax, n);
  return bytes / 1e6;
}

// The thread-scaling table [C11] needs, with both directions side by side.
//
// Side by side is the design, not the presentation. Step H left the lMax = 256
// shortfall with two unseparated candidates -- thread-private accumulators
// competing for last-level cache, and a serialised reduction -- and T10
// removed the second without being able to measure the first. The inverse
// transform's colatitudes write disjoint rows of the field and share only
// read-only input, so it carries no accumulator and does no reduction; the
// forward transform's colatitudes all contribute to every coefficient, so each
// thread accumulates privately and the partials are summed at the end.
// Everything else about the two is the same stream over the same table. The
// fwd/inv column is therefore the accumulator's cost with the rest divided
// out, and its growth with thread count is the answer [C11] is waiting for.
void RunScaling(Int lMax, Int nMax, int windows) {
  const auto n = nMax;
  const auto started = Clock::now();
  auto grid = GaussLegendreGrid<Real, All, All>(lMax, nMax, FFTWpp::Measure);
  const auto buildSeconds =
      std::chrono::duration<double>(Clock::now() - started).count();

  auto field = FFTWpp::vector<Complex>(grid.FieldSize());
  auto coefficients = FFTWpp::vector<Complex>(grid.CoefficientSize(lMax, n));
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    field[i] = Complex{0.5 + 0.001 * i, -0.25 + 0.002 * i};
  }
  const auto bytes = WignerBytes(lMax, n);
  const auto accumulator =
      static_cast<double>(grid.CoefficientSize(lMax, n)) * sizeof(Complex) / 1e6;

  std::printf(
      "\nlMax = %zd, nMax = %zd, complex.  Table %.0f MB built in %.2f s; one\n"
      "transform streams %.0f MB; each thread's forward accumulator is "
      "%.2f MB.\n\n",
      lMax, nMax, TableMegabytes(lMax, nMax), buildSeconds, bytes / 1e6,
      accumulator);
  std::printf("%8s %10s %8s %10s %8s %8s %10s %9s %9s\n", "threads", "fwd(ms)",
              "speedup", "inv(ms)", "speedup", "fwd/inv", "accum(MB)",
              "fwd GB/s", "inv GB/s");

  auto forwardBase = 0.0;
  auto inverseBase = 0.0;
  for (auto threads : ThreadLadder()) {
    const auto policy = threads == 1 ? Execution::Sequential()
                                     : Execution::Parallel(threads);
    const auto forward = TimePerCall(
        [&] {
          grid.ForwardTransformation(lMax, n, field, coefficients, policy);
        },
        0.15, windows);
    const auto inverse = TimePerCall(
        [&] {
          grid.InverseTransformation(lMax, n, coefficients, field, policy);
        },
        0.15, windows);
    if (threads == 1) {
      forwardBase = forward;
      inverseBase = inverse;
    }
    std::printf("%8d %10.3f %7.2fx %10.3f %7.2fx %8.2f %10.1f %9.1f %9.1f\n",
                threads, forward * 1e3, forwardBase / forward, inverse * 1e3,
                inverseBase / inverse, forward / inverse,
                accumulator * threads, bytes / forward / 1e9,
                bytes / inverse / 1e9);
  }
}

// Skip rather than swap: at lMax = 512 the table is 5.4 GB and at 1024 it is
// 43 GB, and a run that starts swapping measures the disk.
bool AffordableAt(Int lMax, Int nMax) {
  const auto needed = TableMegabytes(lMax, nMax) * 1.4;
  const auto available = MemAvailableMegabytes();
  if (available > 0 && needed > available) {
    std::printf(
        "\nlMax = %zd skipped: needs about %.0f MB, MemAvailable is %ld MB.\n",
        lMax, needed, available);
    return false;
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  for (auto i = 1; i < argc; ++i) sectionsWanted.push_back(argv[i]);

  std::printf("GSHTrans transform benchmark (core-plan.md section 5)\n");
  std::printf("double precision, single field per call (k = 1)\n");
  std::printf(
      "sections: stream grid transforms threading server huge "
      "(all, if none named)\n");

  PrintMachineFacts();

  //------------------------------------------------------------------------//
  //                          Bandwidth roofs                                //
  //------------------------------------------------------------------------//

  if (Want("stream") || Want("roof")) {
    PrintHeader("Bandwidth roofs: triad, and a read scan of one transform's table");
    const auto tableBytes = WignerBytes(256, 2);
    std::printf("scan array is %.0f MB, the Wigner bytes one lMax = 256, "
                "n = 2 transform streams\n\n",
                tableBytes / 1e6);
    std::printf("%8s %14s %14s %14s %10s\n", "threads", "triad GB/s",
                "scan GB/s", "scan GB/s", "penalty");
    std::printf("%8s %14s %14s %14s %10s\n", "", "(team touch)", "(team touch)",
                "(1 thread touch)", "");
    for (auto threads : ThreadLadder()) {
      const auto triad = TriadBandwidthGBs(threads, threads, 3);
      const auto scanTeam = ScanBandwidthGBs(tableBytes, threads, threads, 3);
      const auto scanOne = ScanBandwidthGBs(tableBytes, threads, 1, 3);
      std::printf("%8d %14.1f %14.1f %14.1f %9.2fx\n", threads, triad, scanTeam,
                  scanOne, scanOne > 0 ? scanTeam / scanOne : 0.0);
    }
    std::printf(
        "\nThe last two columns differ only in which threads first wrote the\n"
        "pages. The Wigner table is zero-filled by one thread in the grid\n"
        "constructor, so the transform gets the third column, not the second.\n"
        "A penalty well above 1 means the table wants deliberate placement\n"
        "and is not a fact about the transform at all.\n");
  }

  //------------------------------------------------------------------------//
  //                       Grid construction and size                        //
  //------------------------------------------------------------------------//

  if (Want("grid")) {
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

  }

  if (Want("transforms")) {
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

  }

  if (Want("threading")) {
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

  }

  //------------------------------------------------------------------------//
  //             Thread scaling to the full machine ([C11])                  //
  //------------------------------------------------------------------------//

  if (Want("server")) {
    PrintHeader("Thread scaling to the full machine (core-plan.md [C11])");
    std::printf(
        "Step H's decomposition -- colatitudes, with a private accumulator per\n"
        "thread for the forward direction -- was measured to eight threads and\n"
        "is predicted not to survive 64-128: the accumulators become the\n"
        "dominant traffic and the colatitude axis is only lMax + 1 long. These\n"
        "rows are what decides that, and nothing on a laptop can.\n");
    for (auto lMax : {Int{128}, Int{256}, Int{512}}) {
      if (!AffordableAt(lMax, 2)) continue;
      RunScaling(lMax, 2, lMax >= 512 ? 3 : 5);
    }
  }

  // Named explicitly or not run: the table is 43 GB and the grid takes a
  // while to build. Worth one run, because step F' argues its crossover from
  // the table outgrowing last-level cache and this is far past that point.
  if (WantNamed("huge")) {
    PrintHeader("Thread scaling at lMax = 1024");
    if (AffordableAt(1024, 2)) RunScaling(1024, 2, 3);
  }

  if (Want("transforms")) {
    std::printf(
        "\nFFT column is a call truncated to the smallest legal degree: full "
        "FFT\n"
        "work, negligible Legendre work. Leg is the difference. GB/s is the\n"
        "Wigner bytes that call streams divided by the Legendre time.\n");
  }

  FFTWpp::CleanUp();
}
