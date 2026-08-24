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
#include <array>
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
#include <numbers>
#include <optional>
#include <random>
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
  // A 5.4 GB table streamed by every thread is a lot of TLB pressure, and the
  // table is a plain std::vector, so under `madvise` it gets no huge pages.
  std::printf("  huge pages       %s\n",
              FirstLine("/sys/kernel/mm/transparent_hugepage/enabled").c_str());
  std::printf("  cpu governor     %s\n",
              FirstLine("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor")
                  .c_str());
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

// Keep a computed value from being optimised away. One definition rather
// than a volatile sink per section, and it measures nothing itself -- an
// accumulate-into-volatile would add its own arithmetic to the per-point
// columns, which are small enough for that to matter.
template <typename T>
void DoNotOptimise(const T& value) {
#if defined(__GNUC__) || defined(__clang__)
  asm volatile("" : : "m"(value) : "memory");
#else
  static volatile char sink;
  sink = *reinterpret_cast<const volatile char*>(&value);
#endif
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
  // Four accumulators, not one, and this was a defect rather than a
  // refinement. A single `sum += a[i]` is a floating-point dependency chain,
  // and the compiler may not reassociate it without -ffast-math -- so at one
  // thread this measured add latency and not memory at all: 12.7 GB/s against
  // a triad of 36.2 on the same machine, a threefold understatement. It came
  // right at four threads and above, where several chains overlap, which is
  // why the figures section 10 quotes -- all taken at eight threads -- stand.
  // Any single-thread comparison made against the old column was wrong.
  const auto seconds = TimePerCall(
      [&] {
        auto sum = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : sum) \
    num_threads(threads)
        for (std::ptrdiff_t i = 0; i < count; i += 4) {
          auto s0 = ap[i];
          auto s1 = i + 1 < count ? ap[i + 1] : 0.0;
          auto s2 = i + 2 < count ? ap[i + 2] : 0.0;
          auto s3 = i + 3 < count ? ap[i + 3] : 0.0;
          sum += (s0 + s1) + (s2 + s3);
        }
        sink = sum;
      },
      0.15, windows);
  return static_cast<double>(n) * sizeof(double) / seconds / 1e9;
}

// Bytes of Wigner values one transform at this degree and upper index streams:
// the (l, m) block for that upper index, once per colatitude.
// Bytes the *matrix* kernel streams: the same values, but only for m >= 0,
// the rest coming from the reflection (step M6). Order zero is its own
// reflection and is counted once, so this is a shade over half.
double ReflectedWignerBytes(Int lMax, Int n) {
  const auto nTheta = lMax + 1;
  auto values = 0.0;
  for (auto m = Int{0}; m <= lMax; ++m) {
    values += static_cast<double>(lMax - std::max(std::abs(n), m) + 1);
  }
  return values * static_cast<double>(nTheta) * sizeof(Real);
}

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

// The degrees the `server` section walks, chosen against the machine rather
// than hard-coded. 128 to 512 everywhere; 1024 whenever there is room for its
// 43 GB table, because that is the first point far enough past any last-level
// cache to speak to step F'. lMax = 2048 is the `huge` section instead, named
// explicitly: a 344 GB table takes half a minute just to zero, and its
// single-threaded rows are minutes each.
std::vector<Int> ScalingDegrees() {
  auto degrees = std::vector<Int>{128, 256, 512};
  const auto available = MemAvailableMegabytes();
  if (available > 0 && TableMegabytes(1024, 2) * 2 < available) {
    degrees.push_back(1024);
  }
  return degrees;
}

// Fewer repetitions where one call is already long enough to measure.
int Windows(Int lMax) {
  if (lMax >= 1024) return 2;
  if (lMax >= 512) return 3;
  return 5;
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
  // Bumped whenever a section is added or its output changes, so that the
  // wrapper can refuse to run a binary older than the script driving it. Two
  // server runs were lost to exactly that: the source reached the machine with
  // an old timestamp, make saw nothing to do, and the log looked plausible
  // while being produced by the previous harness.
  constexpr auto revision = 9;

  if (argc == 2 && std::string(argv[1]) == "--check") {
    std::printf("harness revision %d\n", revision);
    return 0;
  }

  for (auto i = 1; i < argc; ++i) sectionsWanted.push_back(argv[i]);

  std::printf("GSHTrans transform benchmark (core-plan.md section 5), "
              "harness revision %d\n", revision);
  std::printf("double precision, single field per call (k = 1)\n");
  std::printf(
      "sections: stream grid transforms threading batching generated "
      "kernels kernels-loop kernels-matrix interpolation tuning server huge "
      "(all, if none named)\n");

  // An unoptimised build measures nothing, and the default build directory is
  // one: `cmake -S . -B build` leaves CMAKE_BUILD_TYPE empty, which is fine
  // for the test suite and useless here. The whole harness runs about ten
  // times slow in it -- and every number is self-consistent, so nothing looks
  // wrong. Found the hard way while adding the `kernels` section, whose first
  // run reported speedups of forty.
#ifndef NDEBUG
  std::printf(
      "\n*** WARNING: built without NDEBUG, so probably without "
      "optimisation. ***\n"
      "*** Every figure below is meaningless. Configure with "
      "-DCMAKE_BUILD_TYPE=Release. ***\n");
#endif

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
  //          Generated Wigner values against the stored table (F')          //
  //------------------------------------------------------------------------//

  if (Want("generated")) {
    PrintHeader("Generated Wigner values against the stored table (step F')");
    std::printf(
        "The same recursion, the same order, the same values -- run inside the\n"
        "transform into per-thread scratch instead of read from a table. The\n"
        "two paths agree bit for bit, so the only questions are what it costs\n"
        "and what it saves. Both grids are built in this process and the two\n"
        "are timed alternately, which is the only way a comparison on a\n"
        "clock-scaling machine means anything (core-plan.md section 8).\n\n");

    std::printf("Construction, and what it leaves resident.\n");
    std::printf("%6s %5s %11s %11s %10s %10s\n", "lMax", "nMax", "stored(s)",
                "gen(s)", "storedMB", "genMB");
    for (auto lMax : {Int{64}, Int{128}, Int{256}}) {
      const auto nMax = Int{2};
      if (!AffordableAt(lMax, nMax)) continue;

      // The generating grid first, so that the table's pages are not still
      // resident when its own footprint is measured.
      auto genBefore = ResidentMegabytes();
      auto start = Clock::now();
      {
        auto grid = GaussLegendreGrid<Real, All, All>(
            lMax, nMax, FFTWpp::Estimate, Chunking::Automatic(),
            WignerValues::Generated());
        const auto genSeconds =
            std::chrono::duration<double>(Clock::now() - start).count();
        const auto genMB = ResidentMegabytes() - genBefore;

        auto storedBefore = ResidentMegabytes();
        start = Clock::now();
        auto stored = GaussLegendreGrid<Real, All, All>(lMax, nMax,
                                                        FFTWpp::Estimate);
        const auto storedSeconds =
            std::chrono::duration<double>(Clock::now() - start).count();
        const auto storedMB = ResidentMegabytes() - storedBefore;

        std::printf("%6zd %5zd %11.3f %11.3f %10ld %10ld\n", lMax, nMax,
                    storedSeconds, genSeconds, storedMB, genMB);
      }
    }

    std::printf(
        "\nTransforms, per field. The generated column carries no table\n"
        "traffic at all, so what it competes against is DRAM bandwidth\n"
        "shared between threads rather than arithmetic.\n");
    std::printf(
        "The last column takes the whole batch as one chunk. Generation\n"
        "amortises over a chunk exactly as a table stream does, but the\n"
        "cache-fitting rule that sets the stored optimum is about the\n"
        "coefficient array and says nothing about a table that is not there,\n"
        "so the two paths need not want the same chunk (P2, [C9]).\n");
    std::printf("%6s %4s %4s %8s %10s %11s %11s %9s %11s\n", "lMax", "n", "k",
                "threads", "direction", "stored(ms)", "gen(ms)", "gen/stored",
                "genWhole  storedWhole");
    for (auto lMax : {Int{64}, Int{128}, Int{256}}) {
      const auto n = Int{2};
      if (!AffordableAt(lMax, n)) continue;

      auto stored =
          GaussLegendreGrid<Real, All, All>(lMax, n, FFTWpp::Measure);
      auto generated = GaussLegendreGrid<Real, All, All>(
          lMax, n, FFTWpp::Measure, Chunking::Automatic(),
          WignerValues::Generated());

      const auto fieldSize = static_cast<Int>(stored.FieldSize());
      const auto coefficientSize =
          static_cast<Int>(stored.CoefficientSize(lMax, n));

      for (auto k : {Int{1}, Int{8}}) {
        auto fields = FFTWpp::vector<Complex>(k * fieldSize);
        auto coefficients = FFTWpp::vector<Complex>(k * coefficientSize);
        for (auto i = Int{0}; i < k * fieldSize; ++i) {
          fields[i] = Complex{0.5 + 0.001 * i, -0.25 + 0.002 * i};
        }
        const auto inBatch = Batch::Contiguous(k, fieldSize);
        const auto outBatch = Batch::Contiguous(k, coefficientSize);

        auto wholeChunk = GaussLegendreGrid<Real, All, All>(
            lMax, n, FFTWpp::Measure, Chunking::Fixed(k),
            WignerValues::Generated());

        // The control for it. If taking the whole batch as one chunk costs
        // the forward direction on the stored path too, the cost belongs to
        // the chunk -- each thread's private accumulator is chunk times the
        // coefficient array, which is 8.4 MB per thread at lMax = 256 and
        // k = 8 -- and not to generating anything.
        auto storedWhole = GaussLegendreGrid<Real, All, All>(
            lMax, n, FFTWpp::Measure, Chunking::Fixed(k));

        for (auto threads : {1, 8}) {
          const auto policy = threads == 1 ? Execution::Sequential()
                                           : Execution::Parallel(threads);
          for (const char* direction : {"forward", "inverse"}) {
            auto Time = [&](auto& grid) {
              return TimePerCall([&] {
                if (direction[0] == 'f') {
                  grid.ForwardTransformation(lMax, n, fields, inBatch,
                                             coefficients, outBatch, policy);
                } else {
                  grid.InverseTransformation(lMax, n, coefficients, outBatch,
                                             fields, inBatch, policy);
                }
              });
            };
            const auto a = Time(stored) / static_cast<double>(k);
            const auto b = Time(generated) / static_cast<double>(k);
            const auto c = k > 1 ? Time(wholeChunk) / static_cast<double>(k)
                                 : 0.0;
            std::printf("%6zd %4zd %4zd %8d %10s %11.3f %11.3f %9.2fx", lMax, n,
                        k, threads, direction, a * 1e3, b * 1e3, b / a);
            if (k > 1) {
              const auto e = Time(storedWhole) / static_cast<double>(k);
              std::printf(" %9.3f %11.3f", c * 1e3, e * 1e3);
            }
            std::printf("\n");
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------//
  //             Thread scaling to the full machine ([C11])                  //
  //------------------------------------------------------------------------//

  //------------------------------------------------------------------------//
  //                     Loop kernel against matrix kernel                    //
  //------------------------------------------------------------------------//

#ifdef GSHTRANS_HAVE_BLAS
  // Step M5 of core-plan.md section 11.
  //
  // Three requirements, each guarding a way this measurement can mislead, and
  // section 11.4 states them as part of the step rather than as presentation.
  //
  // -- Two tables, batched and unbatched, and not one with k as a row. The
  // restructure is predicted to buy nothing unbatched at high degree on many
  // threads, where the stage is already at the memory roof. That is the
  // prediction coming true, and a column of 1.0x sitting beside the batched
  // gains reads as failure to anyone scanning it.
  //
  // -- A GB/s column against the roof the `stream` section measures. This is
  // the requirement that does the real work. A caption asserting that no gain
  // is expected is something a reader must take on trust; a row showing 40
  // GB/s against a 42 GB/s roof demonstrates it.
  //
  // -- One kernel per invocation, for a careful A/B. Comparing
  // construction-time policies means two grids, and at lMax = 256 that is
  // 648 MB of table each: both live costs 1.3 GB and the second starts on a
  // cold cache with different page placement. `kernels-loop` and
  // `kernels-matrix` each build one. Plain `kernels` builds both in one
  // process and prints the ratio, which is the convenient form and not the
  // careful one -- it says so in its own output.
  const auto RunKernelSection = [](bool wantLoop, bool wantMatrix,
                                   bool ratios) {
    const auto degrees = std::vector<Int>{64, 128, 256};
    const Int nMax = 2;
    const Int n = 2;

    for (const auto batched : {false, true}) {
      const auto k = batched ? Int{8} : Int{1};
      PrintHeader(batched ? "Kernels, batched (k = 8): the regime phases 2-5 "
                            "run in"
                          : "Kernels, unbatched (k = 1): predicted to gain "
                            "nothing at the roof");
      if (!batched) {
        std::printf(
            "Section 12 predicts no gain here at high degree and many\n"
            "threads, because the loop kernel is already at the memory\n"
            "roof. The GB/s and roof columns are how to check that rather\n"
            "than take it on trust.\n\n");
      }
      std::printf(
          "GB/s is the reported kernel's own table traffic over time --\n"
          "the matrix kernel streams only non-negative orders (M6), a\n"
          "shade over half. Roof is a one-thread-touched read scan. Where\n"
          "the table fits in cache the first exceeds the second, which says\n"
          "the table is not coming from DRAM rather than that anything is\n"
          "wrong.\n\n");
      std::printf("%5s %4s %8s %11s %11s %9s %10s %10s\n", "lMax", "dir",
                  "threads", "loop (ms)", "matrix (ms)", "ratio", "GB/s",
                  "roof");

      for (auto lMax : degrees) {
        const auto bytes = WignerBytes(lMax, n);
        auto loop = std::optional<GaussLegendreGrid<Real, All, All>>{};
        auto matrix = std::optional<GaussLegendreGrid<Real, All, All>>{};
        if (wantLoop) loop.emplace(lMax, nMax, FFTWpp::Measure);
        if (wantMatrix) {
          matrix.emplace(lMax, nMax, FFTWpp::Measure, Chunking::Automatic(),
                         WignerValues::Stored(), TransformKernel::Matrix());
        }
        const auto& any = wantLoop ? *loop : *matrix;
        const auto fieldSize = static_cast<Int>(any.FieldSize());
        const auto coefficientSize =
            static_cast<Int>(any.CoefficientSize(lMax, n));

        auto fields = FFTWpp::vector<Complex>(k * fieldSize);
        for (std::size_t i = 0; i < fields.size(); ++i) {
          fields[i] = Complex{std::cos(0.31 * static_cast<double>(i)),
                              std::sin(0.17 * static_cast<double>(i))};
        }
        auto coefficients = FFTWpp::vector<Complex>(k * coefficientSize);
        const auto fb = Batch::Contiguous(k, fieldSize);
        const auto cb = Batch::Contiguous(k, coefficientSize);

        for (auto threads : ThreadLadder()) {
          const auto policy = threads == 1 ? Execution::Sequential()
                                           : Execution::Parallel(threads);
          const auto roof = ScanBandwidthGBs(bytes, threads, 1, 3);

          for (const auto forward : {true, false}) {
            const auto Run = [&](const auto& grid) {
              return TimePerCall([&] {
                if (forward) {
                  grid.ForwardTransformation(lMax, n, fields, fb, coefficients,
                                             cb, policy);
                } else {
                  grid.InverseTransformation(lMax, n, coefficients, cb, fields,
                                             fb, policy);
                }
              });
            };
            const auto tLoop = wantLoop ? Run(*loop) : 0.0;
            const auto tMatrix = wantMatrix ? Run(*matrix) : 0.0;

            // Table traffic, not work: the Wigner values for this upper
            // index are streamed **once per chunk**, so a batch that fits in
            // one chunk reads them once however many fields it holds. So the
            // figure is bytes over time and not bytes times k -- and read
            // that way it is what the roof column is comparable with. A batch
            // larger than the internal chunk reads the table more than once,
            // which makes this a lower bound there rather than a value.
            // The two kernels stream different amounts: the matrix kernel's
            // table holds only non-negative orders (M6), so reporting the
            // loop kernel's byte count against it would overstate its rate by
            // two and put it above a roof it is nowhere near.
            const auto reported = wantMatrix ? tMatrix : tLoop;
            const auto streamed =
                wantMatrix ? ReflectedWignerBytes(lMax, n) : bytes;
            const auto gbs = reported > 0 ? streamed / reported / 1e9 : 0.0;

            std::printf("%5zd %4s %8d ", lMax, forward ? "fwd" : "inv",
                        threads);
            if (wantLoop) {
              std::printf("%11.3f ", tLoop * 1e3);
            } else {
              std::printf("%11s ", "-");
            }
            if (wantMatrix) {
              std::printf("%11.3f ", tMatrix * 1e3);
            } else {
              std::printf("%11s ", "-");
            }
            if (ratios && tMatrix > 0) {
              std::printf("%8.2fx ", tLoop / tMatrix);
            } else {
              std::printf("%9s ", "-");
            }
            std::printf("%10.1f %10.1f\n", gbs, roof);
          }
        }
      }
    }
    if (ratios) {
      std::printf(
          "\nBoth grids were built in this process, so the second pays a cold\n"
          "cache and whatever placement the first left. For an A/B that turns\n"
          "on tens of percent, run `kernels-loop` and `kernels-matrix` in\n"
          "separate invocations instead.\n");
    }
  };

  if (Want("kernels")) RunKernelSection(true, true, true);
  if (WantNamed("kernels-loop")) RunKernelSection(true, false, false);
  if (WantNamed("kernels-matrix")) RunKernelSection(false, true, false);

  //------------------------------------------------------------------------//
  //          Why the inverse gains less: efficiency against K                //
  //------------------------------------------------------------------------//

  if (Want("kernels")) {
    PrintHeader("Per-order products: Gflop/s against the inner dimension");
    std::printf(
        "M3b attributed the inverse's smaller gain to its inner dimension.\n"
        "The two directions do the same arithmetic on the same matrices:\n"
        "\n"
        "    forward:  (nL x nTheta)(nTheta x 2k)   K = nTheta, a few hundred\n"
        "    inverse:  (nTheta x nL)(nL x 2k)       K = nL, falls to one\n"
        "\n"
        "Equal flops, so Gflop/s compares efficiency directly. If the\n"
        "hypothesis holds, the inverse's rate falls with nL and the\n"
        "forward's does not.\n\n");
    std::printf("%6s %6s %8s %8s %14s %14s\n", "lMax", "order", "nL", "K inv",
                "fwd Gflop/s", "inv Gflop/s");

    const Int lMax = 256;
    const Int n = 2;
    const Int k = 8;
    const auto nTheta = lMax + 1;
    const auto matrices = WignerMatrices<Real, All, All>(
        lMax, lMax, n, std::vector<Real>(nTheta, 0.7));

    auto rhs = std::vector<Complex>(nTheta * k);
    auto coeff = std::vector<Complex>((lMax + 1) * k);
    auto out = std::vector<Complex>(std::max(nTheta, lMax + 1) * k);
    for (std::size_t i = 0; i < rhs.size(); ++i) rhs[i] = Complex{0.5, 0.25};
    for (std::size_t i = 0; i < coeff.size(); ++i) coeff[i] = Complex{0.3, 0.7};

    for (auto m : {Int{0}, Int{64}, Int{128}, Int{192}, Int{240}, Int{254}}) {
      const auto lMin = std::max(std::abs(n), std::abs(m));
      const auto nL = lMax - lMin + 1;
      const auto* a = matrices[n, m].data();
      const auto flops = 2.0 * nL * nTheta * 2 * k;

      const auto fwd = TimePerCall([&] {
        BlasDetails::RowMajorGemm(
            static_cast<int>(nL), static_cast<int>(2 * k),
            static_cast<int>(nTheta), Real{1}, a, static_cast<int>(nTheta),
            reinterpret_cast<const Real*>(rhs.data()), static_cast<int>(2 * k),
            Real{0}, reinterpret_cast<Real*>(out.data()),
            static_cast<int>(2 * k));
      });
      const auto inv = TimePerCall([&] {
        BlasDetails::RowMajorGemmTransposed(
            static_cast<int>(nTheta), static_cast<int>(2 * k),
            static_cast<int>(nL), Real{1}, a, static_cast<int>(nTheta),
            reinterpret_cast<const Real*>(coeff.data()),
            static_cast<int>(2 * k), Real{0},
            reinterpret_cast<Real*>(out.data()), static_cast<int>(2 * k));
      });

      std::printf("%6zd %6zd %8zd %8zd %14.2f %14.2f\n", lMax, m, nL, nL,
                  flops / fwd / 1e9, flops / inv / 1e9);
    }
  }
#endif  // GSHTRANS_HAVE_BLAS

  if (Want("server")) {
    PrintHeader("Thread scaling to the full machine (core-plan.md [C11])");
    std::printf(
        "Step H's decomposition -- colatitudes, with a private accumulator per\n"
        "thread for the forward direction -- was measured to eight threads and\n"
        "is predicted not to survive 64-128: the accumulators become the\n"
        "dominant traffic and the colatitude axis is only lMax + 1 long. These\n"
        "rows are what decides that, and nothing on a laptop can.\n");
    for (auto lMax : ScalingDegrees()) {
      if (!AffordableAt(lMax, 2)) continue;
      RunScaling(lMax, 2, Windows(lMax));
    }
  }

  //------------------------------------------------------------------------//
  //                    Batching, and the chunk heuristic                    //
  //------------------------------------------------------------------------//

  if (Want("batching")) {
    PrintHeader("Batching (core-plan.md step F, tier 1)");
    std::printf(
        "Each row transforms k fields in one call, with the chunk pinned to k\n"
        "so that the row measures one chunk of that width. P2 measured 2.3x at\n"
        "an optimum near k = 8 on a 16 MiB laptop, and *worse than no batching*\n"
        "beyond it. The `auto` column is what Chunking::Automatic would pick\n"
        "here, and the point of these rows is whether it picks near the peak.\n");

    for (auto lMax : {Int{128}, Int{256}}) {
      const auto n = Int{2};
      auto grid = GaussLegendreGrid<Real, All, All>(lMax, n, FFTWpp::Measure);
      const auto fieldSize = static_cast<Int>(grid.FieldSize());
      const auto coefficientSize =
          static_cast<Int>(grid.CoefficientSize(lMax, n));
      const auto bytesPerField =
          coefficientSize * static_cast<Int>(sizeof(Complex));

      std::printf("\nlMax = %zd, n = %zd, complex, sequential.  auto = %zd\n\n",
                  lMax, n, Chunking::Automatic().Count(bytesPerField, 1));
      std::printf("%6s %12s %14s %10s %12s\n", "k", "total(ms)",
                  "per field(ms)", "speedup", "resident(MB)");

      auto base = 0.0;
      for (auto k : {Int{1}, Int{2}, Int{4}, Int{8}, Int{16}, Int{32}}) {
        auto fields = FFTWpp::vector<Complex>(k * fieldSize);
        for (auto i = Int{0}; i < k * fieldSize; i++) {
          fields[i] = Complex{0.5 + 0.001 * i, -0.25 + 0.002 * i};
        }
        auto coefficients = FFTWpp::vector<Complex>(k * coefficientSize);

        // One chunk of width k, whatever the heuristic would have said.
        auto pinned = GaussLegendreGrid<Real, All, All>(
            lMax, n, FFTWpp::Measure, Chunking::Fixed(k));
        const auto inBatch = Batch::Contiguous(k, fieldSize);
        const auto outBatch = Batch::Contiguous(k, coefficientSize);

        auto wholeChunk = GaussLegendreGrid<Real, All, All>(
            lMax, n, FFTWpp::Measure, Chunking::Fixed(k),
            WignerValues::Generated());

        // The control for it. If taking the whole batch as one chunk costs
        // the forward direction on the stored path too, the cost belongs to
        // the chunk -- each thread's private accumulator is chunk times the
        // coefficient array, which is 8.4 MB per thread at lMax = 256 and
        // k = 8 -- and not to generating anything.
        auto storedWhole = GaussLegendreGrid<Real, All, All>(
            lMax, n, FFTWpp::Measure, Chunking::Fixed(k));

        const auto seconds = TimePerCall([&] {
          pinned.ForwardTransformation(lMax, n, fields, inBatch, coefficients,
                                       outBatch);
        });
        const auto perField = seconds / static_cast<double>(k);
        if (k == 1) base = perField;
        std::printf("%6zd %12.3f %14.4f %9.2fx %12.1f\n", k, seconds * 1e3,
                    perField * 1e3, base / perField,
                    static_cast<double>(k) *
                        (fieldSize + coefficientSize) * 16 / 1e6);
      }
    }
    std::printf(
        "\nSpeedup is against k = 1 *per field*, so it is the batching gain\n"
        "and nothing else. A row slower than 1.00x is past the optimum: the\n"
        "chunk's coefficients no longer fit in cache and the Wigner stream is\n"
        "being re-read. If `auto` sits well below the peak, the default cache\n"
        "figure is too conservative for this machine and Chunking::ForCache\n"
        "is what fixes it.\n");
  }

  // Named explicitly or not run. A 344 GB table needs a machine that has it
  // spare and nothing else running, and the single-threaded rows are minutes
  // each. Worth one run on a large server: it is far enough past last-level
  // cache that the stored path has no cache left to lose, which is the regime
  // step F' argues from.
  if (WantNamed("huge")) {
    PrintHeader("Thread scaling at lMax = 2048");
    if (AffordableAt(2048, 2)) RunScaling(2048, 2, Windows(2048));
  }


  //------------------------------------------------------------------------//
  //                    Interpolation (field-algebra-plan.md 22)             //
  //------------------------------------------------------------------------//
  //
  // P5 of section 22, and the step that says whether the local schemes are
  // worth having. Spectral is exact for a band-limited field, so it is the
  // reference the other two are measured against -- which is the whole reason
  // [I2] built it first, and it is what makes the accuracy of a cheap scheme
  // measurable on any field rather than only on one with a closed form.
  //
  // Three columns of error rather than one, because the padding exists for
  // two specific regions and a single number would hide whether it worked:
  // the interior, the two polar cells, and the last longitude cell where the
  // wrap column was added.

  if (Want("interpolation")) {
    PrintHeader("Interpolation: error against the spectral reference");

    // Points off the grid, drawn once and reused at every degree so that the
    // rows are comparable. theta is drawn over the whole of [0, pi], so the
    // polar cells get their share.
    constexpr auto samplePoints = 4000;
    auto engine = std::mt19937_64(20260824);
    auto uniform = std::uniform_real_distribution<Real>(0, 1);
    auto points = std::vector<std::pair<Real, Real>>();
    points.reserve(samplePoints);
    for (auto i = 0; i < samplePoints; ++i) {
      points.emplace_back(
          uniform(engine) * std::numbers::pi_v<Real>,
          uniform(engine) * 2 * std::numbers::pi_v<Real>);
    }

    // The question a caller actually has is not "how good is bicubic at the
    // grid's band limit" -- it is always bad there, since the samples barely
    // resolve the field -- but "how fine a grid do I need". So the band is
    // held fixed and the grid is oversampled, which is the axis that answers
    // it. The first row, oversampling 1, is the band-limit case.
    constexpr Int band = 16;
    constexpr Int N = 2;

    // One set of coefficients, reused at every resolution, so that every row
    // interpolates the *same continuous field* and the rows are comparable.
    // The spectrum falls like 1/(l+1): flat coefficients would put most of
    // the field in its highest degree, which is the least representative case
    // there is.
    auto coefficients = std::vector<Complex>();
    {
      auto normal = std::normal_distribution<Real>(0, 1);
      auto source = std::mt19937_64(11);
      for (auto l = N; l <= band; ++l) {
        const auto scale = 1 / static_cast<Real>(l + 1);
        for (auto m = -l; m <= l; ++m) {
          coefficients.push_back(Complex(normal(source), normal(source)) *
                                 scale);
        }
      }
    }

    std::printf("%6s %7s %10s %14s %14s %14s\n", "grid", "over", "scheme",
                "interior", "polar cells", "last phi cell");

    for (auto factor : {Int{1}, Int{2}, Int{4}, Int{8}}) {
      const auto lMax = band * factor;
      auto grid = GaussLegendreGrid<Real, All, All>(lMax, 2);
      auto expansion = SpinExpansion<N, GaussLegendreGrid<Real, All, All>>(
          grid, band);
      {
        auto next = coefficients.begin();
        for (auto l : expansion.Degrees())
          for (auto m : expansion.Orders(l)) expansion[l, m] = *next++;
      }
      const auto field = Evaluate(expansion);

      // The cell boundaries the padding created.
      auto colatitudes = std::vector<Real>();
      for (auto t : grid.CoLatitudes()) colatitudes.push_back(t);
      auto longitudes = std::vector<Real>();
      for (auto p : grid.Longitudes()) longitudes.push_back(p);
      const auto polar = [&](Real theta) {
        return theta < colatitudes.front() || theta > colatitudes.back();
      };
      const auto lastPhi = [&](Real phi) { return phi > longitudes.back(); };

      // Truncated at the band: above it the coefficients are zero, so this
      // is the same function and a cheaper sum.
      const auto reference = Interpolate(field, Scheme::Spectral(), band);

      // A scale to divide by, so the columns are relative rather than
      // absolute and comparable across degrees.
      auto scale = Real{0};
      for (const auto& [theta, phi] : points) {
        scale = std::max(scale, std::abs(reference(theta, phi)));
      }

      const auto Report = [&](const char* name, auto&& at) {
        auto worst = std::array<Real, 3>{0, 0, 0};
        for (const auto& [theta, phi] : points) {
          const auto error =
              std::abs(at(theta, phi) - reference(theta, phi)) / scale;
          const auto where = polar(theta) ? 1 : (lastPhi(phi) ? 2 : 0);
          worst[static_cast<std::size_t>(where)] =
              std::max(worst[static_cast<std::size_t>(where)], error);
        }
        std::printf("%6td %6tdx %10s %14.2e %14.2e %14.2e\n", lMax, factor,
                    name, worst[0], worst[1], worst[2]);
      };

#ifdef GSHTRANS_HAVE_INTERPOLATION
      Report("bilinear", Interpolate(field, Scheme::Bilinear()));
      Report("bicubic", Interpolate(field, Scheme::Bicubic()));
#else
      std::printf("%6td %6tdx %10s %14s %14s %14s\n", lMax, factor, "(none)",
                  "-", "-", "-");
#endif
    }

    std::printf(
        "\nA band-16 field on grids of 1x to 8x its band, error relative to\n"
        "the largest sampled value. The oversampling column is the answer to\n"
        "\"how fine a grid does a cheap scheme need\"; the first row is the\n"
        "band limit, where the samples barely resolve the field at all.\n"
        "\nThe polar column is what the two extra rows were added for and the\n"
        "last-phi column is what the wrap column was added for; either being\n"
        "far worse than the interior would mean the padding is not doing its\n"
        "job. A bicubic across the wrap is still not a *periodic* spline,\n"
        "which is the one thing the last column can show and an argument\n"
        "cannot.\n");

    //----------------------------------------------------------------------//

    PrintHeader("Interpolation: what it costs, and when to transform instead");

    std::printf("%6s %12s %12s %12s %12s %12s\n", "lMax", "build spec",
                "build bicu", "spectral/pt", "bicubic/pt", "break-even");

    for (auto lMax : {Int{16}, Int{32}, Int{64}, Int{128}}) {
      auto grid = GaussLegendreGrid<Real, All, All>(lMax, 2);
      auto expansion = SpinExpansion<N, GaussLegendreGrid<Real, All, All>>(
          grid, lMax);
      for (auto l : expansion.Degrees())
        for (auto m : expansion.Orders(l)) expansion[l, m] = Complex(1, 0);
      const auto field = Evaluate(expansion);

      const auto buildSpectral =
          TimePerCall([&] {
            auto at = Interpolate(field, Scheme::Spectral());
            DoNotOptimise(at(1.0, 1.0));
          });

      const auto reference = Interpolate(field, Scheme::Spectral());
      auto index = std::size_t{0};
      const auto perPointSpectral = TimePerCall([&] {
        index = (index + 1) & 1023;
        DoNotOptimise(reference(1.0 + 0.0005 * static_cast<Real>(index),
                                0.7));
      });

      // One remesh: expand the field and evaluate it on the same grid, which
      // is what a caller with a whole grid of target points would do instead
      // of asking for a point at a time.
      const auto remesh = TimePerCall([&] {
        auto e = Expand(field, lMax);
        auto f = Evaluate(e);
        DoNotOptimise(f.Data()[0]);
      });

      auto buildBicubic = Real{0};
      auto perPointBicubic = Real{0};
#ifdef GSHTRANS_HAVE_INTERPOLATION
      buildBicubic = TimePerCall([&] {
        auto at = Interpolate(field, Scheme::Bicubic());
        DoNotOptimise(at(1.0, 1.0));
      });
      const auto bicubic = Interpolate(field, Scheme::Bicubic());
      auto j = std::size_t{0};
      perPointBicubic = TimePerCall([&] {
        j = (j + 1) & 1023;
        DoNotOptimise(bicubic(1.0 + 0.0005 * static_cast<Real>(j), 0.7));
      });
#endif

      std::printf("%6td %10.2f ms %10.2f ms %9.2f us %9.3f us %12.0f\n", lMax,
                  buildSpectral * 1e3, buildBicubic * 1e3,
                  perPointSpectral * 1e6, perPointBicubic * 1e6,
                  remesh / perPointSpectral);
    }

    std::printf(
        "\nBreak-even is how many scattered points the spectral interpolant\n"
        "answers in the time one whole-grid remesh takes. Compare it against\n"
        "the grid's own point count: fewer than that and evaluating point by\n"
        "point is the cheaper route, more and it is worth transforming onto a\n"
        "second grid and interpolating there instead.\n"
        "\nBuilding a local interpolant costs a forward transform, because the\n"
        "polar rows are exact ([I3]). Its cheapness is per evaluation, not\n"
        "per interpolant, and these two columns are what says so.\n");
  }


  //------------------------------------------------------------------------//
  //                   Tuning (core-plan.md section 12)                      //
  //------------------------------------------------------------------------//
  //
  // W6: does the mechanism pay? Section 12.5 makes the persistence layer
  // conditional on this answer, so the section reports the two tuners side by
  // side and lets the numbers decide.
  //
  // The answer is machine-dependent by construction -- that is the premise
  // the whole mechanism rests on -- so this exists to be re-run elsewhere
  // rather than to settle anything once.

  if (Want("tuning")) {
    PrintHeader("Tuning: what measuring a policy is worth");

    std::printf("%6s %5s %8s %14s %8s %14s %8s\n", "lMax", "k", "threads",
                "chunk", "cands", "kernel", "chose");

    for (auto lMax : {Int{64}, Int{128}, Int{256}}) {
      for (auto count : {Int{1}, Int{8}, Int{32}}) {
        for (auto threads : {1, 8}) {
          if (threads > omp_get_max_threads()) continue;
          const auto policy = threads == 1 ? Execution::Sequential()
                                           : Execution::Parallel(threads);

          auto grid = GaussLegendreGrid<Real, All, All>(lMax, 2);
          const auto chunk = TuneChunking(grid, lMax, 2, count, policy);
          const auto kernel = TuneKernel<GaussLegendreGrid<Real, All, All>>(
              lMax, 2, 2, count, policy);

          std::printf("%6td %5td %8d %13.2fx %8d %13.2fx %8s\n", lMax, count,
                      threads, chunk.Speedup(), chunk.candidates,
                      kernel.Speedup(),
                      kernel.matrixTried
                          ? (kernel.conclusive ? "matrix" : "loop")
                          : "n/a");
          if (!kernel.skipped.empty()) {
            std::printf("       kernel not tried: %s\n",
                        kernel.skipped.c_str());
          }
        }
      }
    }

    std::printf(
        "\nChunk column is the tuned policy against Chunking::Automatic();\n"
        "cands is how many distinct schedules there were to choose between,\n"
        "and one means the question does not arise -- at k = 1 it never does,\n"
        "since any chunk takes the whole batch in one go.\n"
        "\nKernel column is the matrix kernel against the loop, both built\n"
        "and timed here. Section 12's [C22] margin is ten per cent, so a\n"
        "column reading 1.00x means the incumbent held rather than that the\n"
        "two were identical.\n"
        "\nBoth are the same measurement a caller would make at start-up, on\n"
        "their own problem shape. The point of the section is that the answer\n"
        "differs by machine, so it is worth re-running here rather than\n"
        "reading the figures the plan quotes for a laptop.\n");
  }

  if (Want("transforms")) {
    std::printf(
        "\nFFT column is a call truncated to the smallest legal degree: full "
        "FFT\n"
        "work, negligible Legendre work. Leg is the difference. GB/s is the\n"
        "Wigner bytes that call streams divided by the Legendre time.\n");
  }
}
