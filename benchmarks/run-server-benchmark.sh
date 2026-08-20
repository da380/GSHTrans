#!/usr/bin/env bash
#
# One run of the transform benchmark on a target machine, with the machine's
# own description recorded next to it.
#
#   ./benchmarks/run-server-benchmark.sh [output-directory]
#
# Everything lands in one text file, which is the thing to send back.
#
# What this exists to settle (core-plan.md):
#
#   [C11]  Step H threads the transform over colatitudes, with a private
#          accumulator per thread in the forward direction. That was measured
#          to eight threads on a laptop and is predicted not to survive
#          64-128, because the accumulators become the dominant memory traffic
#          and the colatitude axis is only lMax + 1 long. The alternatives --
#          m-block threading, threading over the batch axis, a two-dimensional
#          split -- cannot be chosen without these numbers, and choosing on
#          laptop numbers would repeat the error section 9 records.
#
#   step H The lMax = 256 shortfall (2.7x against 4.2x at 128) had two
#          unseparated candidates. T10 removed one of them, below the laptop's
#          noise floor. The fwd/inv column separates the other.
#
#   step F P8 sizes a transform's chunk as perCoreL3 / (2 * 16 * nCoefficients).
#          The header records this machine's per-core L3, which is the input.
#
#   NUMA   Not yet in the plan, and testable here for the first time. The
#          Wigner table is a std::vector<Real> built by its size constructor
#          (Wigner.h:131), so it is zero-filled by the single constructing
#          thread and every page first-touches on that thread's NUMA node.
#          On a multi-socket machine all threads then stream one node's
#          memory. Runs 2 and 3 below are the control.

set -u

outdir="${1:-.}"
mkdir -p "$outdir"
stamp="$(date +%Y%m%d-%H%M%S)"
host="$(hostname -s 2>/dev/null || echo unknown)"
log="$outdir/gshtrans-bench-$host-$stamp.txt"

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build="$root/build-bench"
binary="$build/bin/TransformBenchmark"

exec > >(tee "$log") 2>&1

echo "GSHTrans server benchmark"
echo "host    $host"
echo "date    $(date -Is)"
echo "source  $root"
echo "commit  $(git -C "$root" rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "        $(git -C "$root" status --porcelain 2>/dev/null | wc -l) modified files"
echo

echo "=============================================================================="
echo "Machine, as the operating system describes it"
echo "=============================================================================="
command -v lscpu   >/dev/null && lscpu
echo
if command -v numactl >/dev/null; then
  numactl --hardware
else
  echo "numactl not installed: the NUMA control run below will be skipped."
  echo "(apt install numactl / dnf install numactl -- it is worth having.)"
fi
echo
grep -E 'MemTotal|MemAvailable' /proc/meminfo
echo
echo "compiler: $(${CXX:-c++} --version 2>/dev/null | head -1)"
echo

echo "=============================================================================="
echo "Build"
echo "=============================================================================="
#
# Release is -O3 -DNDEBUG. NDEBUG only turns off the point-index asserts; the
# transform's size checks have thrown in every build mode since T4, so nothing
# that guards correctness is being compiled out.
#
# -march=native, because the decisions these numbers feed -- step F's chunk
# size, [C11]'s decomposition -- are decisions about this machine, so it should
# be compiled the way it will be deployed. The one machine where that is not
# obviously right is a Skylake-SP or Cascade Lake Xeon, where heavy AVX-512
# pulls the clock down and buys nothing for load/store-bound work; if lscpu
# above shows one, run again with GSH_CXX_FLAGS="-g -fno-omit-frame-pointer"
# and compare.
#
# -g -fno-omit-frame-pointer costs nothing at -O3 and means perf record gives
# usable stacks if one of the curves below looks strange.
#
# Deliberately absent: -Ofast and -ffast-math. The Wigner recursions are the
# numerically delicate part of this library and the round-trip tolerance is the
# only oracle for them, so reassociation would change the thing being
# validated, for no gain on a kernel that is not flop-bound.
#
: "${GSH_CXX_FLAGS:=-march=native -g -fno-omit-frame-pointer}"
echo "flags: CMAKE_BUILD_TYPE=Release CMAKE_CXX_FLAGS=\"$GSH_CXX_FLAGS\""
echo
if [ ! -d "$build" ]; then
  if ! cmake -S "$root" -B "$build" \
         -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_CXX_FLAGS="$GSH_CXX_FLAGS" \
         -DMY_PROJECT_BUILD_EXAMPLES=OFF; then
    echo
    echo "Configure failed. The usual cause on a fresh machine is that"
    echo "CMakeLists.txt fetches GaussQuad over SSH (git@github.com:...), so it"
    echo "needs a GitHub key this machine may not have. Check with:"
    echo "    ssh -T git@github.com"
    echo "OpenMP and double-precision FFTW3 must also be installed."
    exit 1
  fi
else
  # An existing build directory keeps whatever it was configured with, so
  # report that rather than what was asked for above.
  echo "reusing $build, configured as:"
  grep -E '^(CMAKE_BUILD_TYPE|CMAKE_CXX_FLAGS):' "$build/CMakeCache.txt" 2>/dev/null |
    sed 's/^/  /'
  echo "  (delete it to reconfigure)"
fi
cmake --build "$build" --target TransformBenchmark -j "$(nproc)" || exit 1

# Refuse to measure with a binary older than this script.
#
# Two server runs were lost this way. The source reached the machine carrying
# an older timestamp than the object file already in the build directory, make
# reported "Built target" without compiling anything, and the log that came
# back looked entirely plausible while having been produced by the previous
# harness. Nothing in the output said so, which is the part worth fixing.
expected=4
got="$("$binary" --check 2>/dev/null | awk '/harness revision/ {print $3}')"
if [ "${got:-0}" -lt "$expected" ]; then
  echo
  echo "STALE BINARY: this script needs harness revision $expected, the built"
  echo "binary reports ${got:-none}. Nothing below would mean anything, so stopping."
  echo
  echo "The source did not reach the build. Usually its timestamp is older than"
  echo "the object file, so make sees nothing to do. Either:"
  echo "    touch $root/benchmarks/TransformBenchmark.cpp && $0 $outdir"
  echo "or, to be certain:"
  echo "    rm -rf $build && $0 $outdir"
  echo
  echo "Check the source is actually there first:"
  echo "    grep -c PrintMachineFacts $root/benchmarks/TransformBenchmark.cpp"
  exit 1
fi
echo "harness revision $got"
echo

echo "=============================================================================="
echo "Correctness first: the suite must pass before any timing means anything"
echo "=============================================================================="
cmake --build "$build" -j "$(nproc)" >/dev/null 2>&1
(cd "$build" && ctest --output-on-failure 2>&1 | tail -20)
echo

# Cores, not hardware threads: for this library's memory-bound work the second
# thread on a core measured slower than one thread per core. The ladder inside
# the benchmark walks past this anyway; binding is about placement, not count.
run () {
  local title="$1"; shift
  echo
  echo "=============================================================================="
  echo "$title"
  echo "=============================================================================="
  echo "+ $*"
  echo
  "$@"
}

run "Run 1 of 3 -- threads bound to cores, memory wherever it falls" \
  env OMP_PROC_BIND=spread OMP_PLACES=cores "$binary" stream batching server

run "Run 2 of 3 -- unbound, which is what an unprepared caller gets" \
  env OMP_PROC_BIND=false "$binary" stream batching server

if command -v numactl >/dev/null && [ "$(ls -d /sys/devices/system/node/node* 2>/dev/null | wc -l)" -gt 1 ]; then
  run "Run 3 of 3 -- pages interleaved across nodes, the NUMA control" \
    env OMP_PROC_BIND=spread OMP_PLACES=cores \
    numactl --interleave=all "$binary" stream batching server
else
  echo
  nodes="$(ls -d /sys/devices/system/node/node* 2>/dev/null | wc -l)"
  if [ "$nodes" -le 1 ]; then
    echo "Run 3 skipped: one NUMA node, so there is nothing to control for."
  else
    echo "Run 3 NOT RUN, and it is the one that matters most on this machine."
    echo "This host has $nodes NUMA nodes and numactl is not installed. The"
    echo "Wigner table is zero-filled by the single constructing thread, so"
    echo "every page of it lands on one node and all threads then stream one"
    echo "node's memory. Whether that is what limits the scaling below cannot"
    echo "be told apart from the accumulators without this run."
    echo "    apt install numactl   (or dnf install numactl), then run again."
  fi
fi

echo
echo "=============================================================================="
echo "Done. Send back: $log"
echo "=============================================================================="
echo
echo "Optional, and worth one run if the machine is idle and has the memory:"
echo "  $binary huge        # lMax = 1024, a 43 GB table, ~10 minutes"
echo "It is the only point far enough past last-level cache to speak to step F'."
