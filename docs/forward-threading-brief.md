# Brief: GSHTrans forward transform threading (from the sl3dpp optimisation pass)

Repo: /space/da380/dev/GSHTrans, branch `develop` (fresh, cut from main at
991dd1e, which is sl3dpp's pin). Consumer: /space/da380/dev/sl3dpp.

## House rules
- **David makes all git commits himself. Never commit or push.** Leave changes
  in the working tree and summarise them.
- tunya is shared with students. Builds/tests: around 16-24 cores
  (`cmake --build -j 16`, `ctest -j 16`). Timing runs: up to 32 threads (one
  socket), short and solo, `OMP_PLACES=cores OMP_PROC_BIND=close`; check
  `uptime` first. Never more.
- Machine: 2 x EPYC 9334 (Zen 4, 32 cores/socket, 2 NUMA nodes, 128 logical
  CPUs; logical 0-31 = socket 0 physical cores). GCC 13.2 (/opt/gcc-13.2.0),
  THP is `always`. No `perf` (paranoid 3); gdb works for stack sampling.

## The problem (measured 5 Nov 2026, sl3dpp `step_timing`, lmax 128)
sl3dpp's 3D step is transform-bound. Its grid is
`GaussLegendreGrid<double, All, All>` with `Chunking::Automatic()`,
`WignerValues::Stored()`, the loop kernel, `Execution::Parallel()`; it calls
the batched `ForwardTransformation` / `InverseTransformation` with
`Batch::Contiguous(count, ...)` over hundreds of radial nodes at once, spins
n = 0 (`SynthesiseBatch`) and the spin components (`SynthesiseSpinBatch`, N
from sl3dpp's kSpins), complex fields.

Per pass of the step, one process, seconds:

| threads | analyse (forward) | synthesise (inverse) |
|---|---|---|
| 16 | 1.70 | 1.27 |
| 32 | 1.49 | 0.77 |

Same 32 cores split into MPI ranks over radial nodes (each rank transforms
1/R of the batch on 32/R threads):

| ranks x threads | analyse | synthesise |
|---|---|---|
| 1 x 32 | 1.52 | 0.76 |
| 2 x 16 | 0.87 | 0.61 |
| 4 x 8  | 0.47 | 0.39 |
| 8 x 4  | 0.41 | 0.39 |

So the forward transform scales with the work per process, not with threads
inside a process. The inverse scales.

## Diagnosis (from reading SphericalGrid.hpp; confirm by measurement)
`ForwardLoopKernel` threads over colatitudes with a private accumulator per
thread and a partitioned reduction. `ForwardChunkSize` passes the *thread
count* as `copies` to `Chunking::Count`, so with the 8 MB default cache:
coefficient block per field at lmax 128 ~ 129^2 x 16 B = 266 KB; share = 8 MB /
threads; at 16 or 32 threads `Count` returns **1**. The inverse passes copies
= 1 and gets ~16.

Consequence, per chunk of ONE field: open a parallel region, zero `threads`
accumulators of 266 KB, stream the whole stored Wigner table for that n, and
reduce threads x 266 KB. Total Wigner traffic and reduction work are both
proportional to the thread count, so per-thread work is constant: no scaling.
(docs/open-issues.md §4, first bullet, anticipated this.)

## Suggested design (David agreed to GSHTrans work; the design is yours to
## confirm or improve)
A forward kernel that mirrors the inverse: two stages per chunk,
1. **Fourier stage**, threaded over colatitudes (disjoint writes): pack, FFT,
   copy each row's [m][k] output into a chunk buffer [theta][m][k]
   (or reuse `ForwardFourierStage`'s machinery/layout if it fits).
2. **Legendre stage, threaded over blocks of degrees** (the coefficient
   layout is degree-major, so a degree block is a contiguous run of
   coefficients): each thread loops theta = 0..nTheta-1 in order and, for its
   degrees, does exactly the existing inner axpy into its own disjoint slice
   of the [coefficient][field] scratch. No private accumulators, no
   reduction, each thread reads only its degrees' slice of each theta's
   Wigner row. Balance blocks by coefficient count (or more blocks than
   threads with a dynamic schedule).

Properties to keep/gain:
- Each coefficient is summed over theta in the same order as the sequential
  path, so the parallel result should be **bitwise equal to the sequential
  one, independent of thread count**. Add a test for that. (sl3dpp has
  determinism gates that will benefit.)
- The chunk then needs one live copy of the coefficient block plus the
  Fourier buffer, so `ForwardChunkSize` can use copies = 1 like the inverse
  (revisit the cache budget for the Fourier buffer: c x nTheta x nFourier
  complex).
- Generated Wigner values: a degree block would need the recursion run per
  thread per theta from |n| up, so keep the existing colatitude kernel for
  Generated() (or measure before deciding).
- Respect the existing structure: one level of threading
  (`Execution::TeamSize`, `RunInParallel`), `Details::ExceptionCapture`
  inside regions, thread_local workspaces, real-field (m >= 0) and complex
  paths, `MRange_` variants, the matrix (BLAS) kernel untouched.

Other levers recorded by David for later (not this task unless cheap): a fused
synthesise-weight-analyse transform; a real spin transform producing m >= 0
coefficients directly; a ragged batch for depth-truncated node blocks.

## Measure before and after
- GSHTrans's own `benchmarks/TransformBenchmark.cpp` (and
  `run-server-benchmark.sh`): forward and inverse, lmax 64 and 128, batch of
  a few hundred fields, n = 0 and n = +-1/+-2, threads 1, 8, 16, 32.
- Then sl3dpp: configure a build against the working copy with
  `-DFETCHCONTENT_SOURCE_DIR_GSHTRANS=/space/da380/dev/GSHTrans` (a separate
  build dir, e.g. sl3dpp/build-gsh, so David's `build/` stays at the pin), and
  rerun
  `OMP_NUM_THREADS=32 OMP_PLACES=cores OMP_PROC_BIND=close build-gsh/benchmarks/step_timing tests/data/prem.simp.vm1 128 2 2`
  (baseline above: per pass 2.77 s, analyse 1.49 s, startStep 3.44 s,
  peak RSS 34 GB) plus the 16-thread run. sl3dpp's test suite
  (`ctest -j 24`, ~3 min) must stay green; its threaded determinism
  expectations are documented in sl3dpp README "Reproducibility".
- GSHTrans tests green; add the bitwise sequential-vs-parallel test.
- Record the measurements in GSHTrans docs (open-issues.md §4 is the natural
  place to answer) in the repo's dated, measured style.
