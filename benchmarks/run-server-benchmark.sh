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
if [ ! -d "$build" ]; then
  cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release || exit 1
fi
cmake --build "$build" --target TransformBenchmark -j "$(nproc)" || exit 1
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
  env OMP_PROC_BIND=spread OMP_PLACES=cores "$binary" stream server

run "Run 2 of 3 -- unbound, which is what an unprepared caller gets" \
  env OMP_PROC_BIND=false "$binary" stream server

if command -v numactl >/dev/null && [ "$(ls -d /sys/devices/system/node/node* 2>/dev/null | wc -l)" -gt 1 ]; then
  run "Run 3 of 3 -- pages interleaved across nodes, the NUMA control" \
    env OMP_PROC_BIND=spread OMP_PLACES=cores \
    numactl --interleave=all "$binary" stream server
else
  echo
  echo "Run 3 skipped: one NUMA node, or numactl absent. Nothing to control for."
fi

echo
echo "=============================================================================="
echo "Done. Send back: $log"
echo "=============================================================================="
echo
echo "Optional, and worth one run if the machine is idle and has the memory:"
echo "  $binary huge        # lMax = 1024, a 43 GB table, ~10 minutes"
echo "It is the only point far enough past last-level cache to speak to step F'."
