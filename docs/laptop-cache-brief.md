# Brief: the default chunk cache figure, measured on a laptop

Repo: GSHTrans, branch `develop` (pull the latest first). Written 5 November
2026 on the server, for an agent on David's laptop.

## Rules
- **Never commit or push.** David commits himself. Leave any edits in the
  working tree and summarise them at the end.
- This is a measurement task. Don't change library code unless asked.
- Laptop timings are noisy (clock scaling). Before timing: plugged in, other
  heavy apps closed, performance power profile if there is one. Compare
  figures only within one run, or across runs made back to back.

## Background (read first)
- `docs/design.md` §2.3–2.4 and `docs/open-issues.md` §4 ("The default cache
  figure"). Measured 5 Nov 2026 on a 2 × EPYC 9334 server: the default
  `Chunking::DefaultCacheBytes` (8 MiB, `GSHTrans/Policies.hpp`) costs
  14–26% at bands 256–384 against `Chunking::ForCache(32 MiB)`, and nothing
  on that machine argued for keeping it small. The 8 MiB figure was
  originally justified by laptop measurements under *older* kernels: each
  thread kept its own forward accumulator, and the inverse had a serial
  gather. Both are gone now. A floor of four
  (`Chunking::MinimumSharedCount`) already covers the worst case, band 512.
- The question for this run: on a laptop-sized L3 (8–16 MB), does 32 MiB (or
  something between 8 and 32) cost anything with the current kernels? If not,
  the default can be raised.

## Build and test

    cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DGSHTRANS_BUILD_EXAMPLES=OFF
    cmake --build build-bench -j <physical cores>
    cd build-bench && OMP_NUM_THREADS=2 ctest -j <physical cores>

GCC 13+ or Clang 18+. Cap `OMP_NUM_THREADS` for ctest: tests otherwise start
a full OpenMP team each and oversubscribe the machine. All tests should pass
before any timing.

## Measure
Binary: `build-bench/bin/TransformBenchmark`. Use `OMP_PLACES=cores
OMP_PROC_BIND=close`. Set `GSHTRANS_BENCH_THREADS` to the thread counts worth
seeing on this machine, e.g. `4,8` and the physical core count. The output
header records the machine's L3; note it.

1. `TransformBenchmark degree-sweep`: bands 96–512, ms per field, the default
   (`fwd auto` / `inv auto`) against `ForCache(32 MiB)` (`ccd`). Check
   `MemAvailable` first: band 512 needs about 3.6 GB for its table plus about
   2.5 GB of fields. If memory is short, stop after band 384 and say so.
2. `TransformBenchmark chunk-sweep`: Fixed chunks 1–63 at lmax 64 and 128,
   both directions. Shows where this machine's optimum chunk is, and whether
   large chunks get slower. The heuristic's chunk for a cache figure C is
   round((C + b) / (2 b)), with b = one field's coefficient bytes (266192 at
   band 128, 1056720 at band 256), clamped to [4, 63]. So any cache figure
   can be read off this table.
3. Run each twice, back to back, to see the run-to-run noise.

## Decide and report
- For each band and thread count: is `ccd` (32 MiB) slower than `auto`
  (8 MiB) by more than the run-to-run noise? Where, and by how much?
- From chunk-sweep: is there a chunk beyond which either direction gets
  slower on this machine, and what cache figure would that correspond to?
- Recommend one of: raise `DefaultCacheBytes` to 32 MiB; raise it to an
  intermediate figure (say which, and why); or keep 8 MiB (say what it
  protects against).
- Add the measurements, dated and with the machine described (CPU, L3,
  cores, compiler), to `docs/open-issues.md` §4 "The default cache figure",
  in the style of the existing entries. Leave the edit uncommitted.
- Finish with a short summary: machine, the key table rows, the
  recommendation, and the list of files changed.
