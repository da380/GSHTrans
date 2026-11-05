# Open issues and planning

The one place for what is not done: suspected problems, known limitations,
open design questions and possible work. Everything else about the library —
what it does and why — is in `design.md`, `gshtrans-reference.tex` and
`README.md`. When an item is resolved, delete it here and, if the resolution
is a design decision, record the reason in `design.md`.

Each item says what it is and where; search for the named symbol.

---

## 1. Optional clean-up

- `Wigner.hpp`: the closed-form boundary functions (`WignerMinOrder` and
  siblings) are used only by `tests/CheckWignerBoundary.hpp`; they could live
  in the test tree.

## 2. Documentation to check

Statements that may be wrong and need the owner's judgement or a measurement.

- **Headroom for integration.** `GaussLegendreGrid.hpp` (`ForBand`) and the
  reference ("Quadrature") give integrating |f|² of a band-L field as a reason
  for oversampling. A grid at lMax = L already integrates it exactly
  (L + 1 Gauss nodes are exact to degree 2L + 1; nPhi ≥ 2L + 1). Headroom is
  needed to *represent or transform* a product, which has band 2L.
- **Two timings for one configuration** in the reference: lMax = 256, n = 2,
  one thread, forward is ~19 ms in "The transform" and 12.16 ms in "Batching,
  chunking and threading". Re-measure in one run.
- **A citation with no bibliography.** The reference cites "Woodhouse \&
  Deuss [41]" but has no bibliography; give the reference in full.
- **Radial lines with the loop kernel**: `run-server-benchmark.sh` says the
  direct route lost everywhere with the threaded loop kernel; the benchmark's
  `lines` comment says the routes were close. One is out of date.

## 3. Known limitations (by decision)

Deliberate. Listed so that they are not rediscovered as bugs, and so the
decision can be revisited if the use changes.

**Numerics and transforms**

- **Degree ceiling.** `lMax ≤ SafeDegree<Real>(mMax, nLow, nHigh)` — 1822 in
  double for every order at `n = 0`, higher for restricted orders, 913 at
  worst — set by underflow of the Wigner seed (`design.md` §5.1). Lifting it
  needs an exponent carried per column (Fukushima, *J. Geod.* 86, 2012). Out
  of scope unless a use above lMax ≈ 1800 in double with every order appears;
  long double reaches 30 740, and a table of a few orders goes much higher in
  double already.
- **6-j symbols are not implemented.** Schulten and Gordon give the recursion
  in the same paper, and the 3-j machinery would carry over.
- **The 3-j `Fill…` functions allocate one row of scratch per call.** Calling
  them in a hot threaded loop would want the scratch passed in.
- **ILP64 BLAS is not detected.** `Blas.hpp` assumes 32-bit integers; an ILP64
  BLAS would link and pass wrong arguments. No portable header-side check
  exists.
- **No persistence for tuning results.** Tuning is cheap at start-up up to
  lMax of a few hundred; around lMax = 512 building the tables to measure
  them stops being cheap, and a persisted store would become worth having.
- **Reduced-precision storage of a double table** is not offered; any
  experiment would first need a tighter round-trip oracle.

**Tensors**

- **One slot alphabet per tensor.** Mixed objects such as Dahlen & Tromp's
  T^{rΩ} cannot be one tensor.
- **No trace-free symmetric tensors.** The orbit machinery handles signed slot
  permutations only, so the trace of a tangential symmetric (spin-2) tensor
  cannot be removed from storage; the caller subtracts it.
- **`Materialise<Symmetry, Reality>` does not check** that the value has the
  symmetry or reality stated: only representatives are copied, so the
  difference is discarded. A check would need every derived component and a
  tolerance.
- **Tensor `Materialise` has no contiguous fast path**; it copies element by
  element even into a component-major target where `EvaluateInto` on a
  contiguous span would do.

**Layered fields**

- **Radial operators are conveniences.** They refuse what they cannot do
  correctly (narrow elements, Lagrange across interfaces) rather than gaining
  the capability (`design.md` §6).
- **Spline end conditions apply at every interface.** `SplineDerivative`'s
  `left` and `right` are imposed on every element, so the default `Natural`
  forces zero curvature on both sides of each interface. Interior pieces
  cannot be given different conditions.
- **`Gradient` assumes a radial operator with real weights.** Derived
  components are reconstructed by conjugation, which commutes with the
  operator only if its weights are real. Every ready-made operator's are; a
  caller's complex-weighted one would give wrong derived components silently.
- **The layered `Gradient` refuses an operand too short to hold its result.**
  Its radial half maps stacks to stacks of the same shape, so it cannot change
  the degree as the flat `SurfaceGradient` does.
- **The layered `Gradient` threads only its radial half**, and
  `SurfaceGradient` takes no `Execution` policy.
- **No automatic choice between `ExpandToLines` and `Expand` + `RadialMajor`.**
  The timing difference depends on `nR` (power-of-two strides) and threads in
  a way no simple rule captures.

## 4. Open questions for the 64-core target

These need measurements on the production machine (dual-socket, 64 cores)
rather than the development laptop.

- **The forward loop kernel's decomposition** is answered for stored tables
  on one socket (`ForwardBlockKernel`; measurements in `design.md` §2.3).
  The inverse's limit was a serial gather, now divided between the team
  (same section). Still open: (a) a threaded forward on *generated* values
  still keeps private accumulators — the candidate is to generate each
  colatitude's row once per chunk, threaded over colatitudes, into a
  chunk-local table and run the
  degree-block Legendre stage on it (one n of the table per chunk: 26 MB at
  lMax 128, 8.6 GB at 1024); (b) one process across two sockets is slower
  than on one, in both directions, which is not this kernel's doing — a
  two-dimensional split (groups of fields × blocks of coefficients) would
  only matter if a single process must span sockets, and nothing needs that
  yet.
- **The chunk rule on multi-cache-domain machines.** `Chunking::Count` counts
  the inverse transform's shared block as one copy, but on a multi-CCD or
  multi-socket machine it is pulled into every cache domain that reads it, so
  the rule is optimistic there. On one Zen 4 socket (four 32 MB L3s) this
  measured as nothing: chunks from 16 to 63 are within about 8% in both
  directions at lMax 64 and 128 (`TransformBenchmark chunk-sweep`, 5 November
  2026). Whether it matters across sockets is moot while one process per
  socket is the configuration.
- **The default cache figure.** `Chunking::DefaultCacheBytes` (8 MiB) costs
  14–26% at bands 256–384 against 32 MiB on one Zen 4 socket, and nothing
  measured on this machine argues for keeping it small now that neither
  direction keeps per-thread copies of a stored-table block (`design.md`
  §2.4). The floor of four covers the collapse at band 512. Run
  `TransformBenchmark chunk-sweep` and `degree-sweep` on a laptop (8–16 MB
  L3) with the current kernels; if 32 MiB costs nothing there either, raise
  the default.
- **Generated Wigner values across NUMA domains.** Whether generating per
  thread beats streaming a stored table from a remote node has been argued,
  not measured.
- **`WignerMatrices` construction** scatters at stride `NumberOfAngles()`, up
  to a cache line per value; blocking over colatitudes would remove that if
  table construction time matters at large lMax.

## 5. Open design questions

Each wants a decision rather than an edit.

- **`Scheme`** (`Policies.hpp`) is shaped unlike the other policy values: a
  class with no instances whose named constructors return distinct tag types,
  because the schemes have different interpolant types and dispatch must be by
  overload. Rename it, or accept the exception.
- **`Chunking::Count`** takes `int copies` where everything around it uses
  `std::ptrdiff_t`.
- **`SphericalGrid.hpp`** is by far the largest header (grid, both kernels,
  Fourier stages, thread caches). A split along those lines is possible.
- **The two layered field types** (`LayeredSpinField`, `LayeredTensorField`)
  duplicate a good deal of code. The part that was error-prone (derived
  components) is shared; the rest is not.
- **A transform plan object**, FFTW-style (validate once, allocate outside the
  parallel region, carry tuning results). Not built because it cuts across
  keeping policies on the grid with `With()`. Revisit if `Batch::At` grows.

## 6. Possible work

- **Race-checking the threaded paths.** ThreadSanitizer is unusable with GCC's
  `libgomp`; a Clang build against `libomp` built with `LIBOMP_TSAN_SUPPORT`
  would make it possible.
- **A test that throws inside a transform kernel's or `Resample`'s OpenMP
  region.** Neither has a natural way to throw any more: the kernels since
  `FFTWpp::WisdomOnly` is refused, and `Resample` since it checks its pieces
  before the loop. Their `ExceptionCapture` wrapping is covered by the
  utility's own tests, and by `ApplyRadially` with a throwing operator.
- **Release pinning.** A release commit should set the four sibling
  dependencies (`GSHTRANS_FFTWPP_TAG` and siblings) to fixed revisions;
  development tracks their `main`. googletest is pinned to a release.
- **The field layers at an application's scale** (planned 5 November 2026,
  after the transform threading work). sl3dpp hand-codes its pseudospectral
  products — a five-real deviatoric strain per point, stacked over ~10³
  radial nodes, multiplied by a pointwise weight between batched spin
  transforms — rather than using `SpinField`, the layered fields or
  `TensorField`. Its view of why, and what would change that, is in
  `docs/gshtrans-notes-2026-11-04.md` §3c; the first step is a read-only
  assessment of those headers against sl3dpp's `GridCoupling` that confirms
  or corrects it. The candidates, in that note's order of payoff:
  1. a layered `Expand` that takes a layered expression node (a stack times
     a pointwise weight, say), evaluated slice by slice into the batch
     buffer before one batched transform — removes a pass over a 4 GB grid
     tensor at lMax 128 and sl3dpp's hand-written multiply loops;
  2. caller-supplied buffers (`Into` forms) on the layered transforms;
  3. a ragged batch (a per-slice degree), for depth-truncated node blocks;
  4. `GridCoupling`'s tensor as a `LayeredTensorField<2, Symmetric, Real>`,
     worth it when a rank-4 pointwise weight arrives (the relabelling work),
     not before.
  Beyond item 1 is a fully fused synthesise → weight → analyse kernel that
  never holds more of the grid tensor than one chunk of latitude rows: a
  kernel project, not an API one, whose case is memory traffic and peak
  memory (34 GB in sl3dpp's step at lMax 128), not flops — the pointwise
  multiply itself is 0.03 s of a 1.57 s pass. The assessment should estimate
  what each would save on that step before either is built.
- **A small parallel-loop layer** ("RAJA, but simpler"), to be designed with
  the field-layer work above. The kernels repeat one pattern — the
  `RunInParallel` decision, `num_threads(ThreadCount(policy))`,
  `ExceptionCapture` around set-up, a static split of a range, every thread
  reaching the barrier — and pointwise work over stacked grids (a product at
  every node of every slice, local tensor contractions) will repeat it
  again. Sketch: `ForEach(policy, n, body(i))`; `ForEachBlock(policy, n,
  body(begin, end))`, the form that vectorises; a stacked form that flattens
  radius × point into one range so a 10³-slice stack is one loop, not 10³;
  and `Reduce` with a fixed partition and fixed combine order, so results
  are bitwise independent of the thread count. It keeps the library's rules
  (one level via `Execution::TeamSize`, exceptions captured, static
  contiguous blocks for first touch, nothing when OpenMP is off), lives in
  Core beside `OpenMP.hpp`, and is a few function templates, not a framework.
  `std::execution::par` is not the route: it needs TBB under GCC and knows
  nothing of either rule. A layered `Expand(node)` is the first user: its
  evaluation into the batch buffer is exactly the stacked loop. Start by
  listing the loops in GSHTrans and sl3dpp that would use it.

  The layer is for user code as much as for the library, and serves two
  separate aims. One is *functionality users need not write*: the loops,
  the deterministic `Reduce`, and the threading and exception rules that
  come with them, so an application's pointwise work is a lambda, not a
  hand-written region. The other is *working with the user's own
  structures*: the loops take the application's layouts, strides, buffers
  and decomposition, not only the library's types, so a code with its own
  stacks or its own MPI layout can use them without first converting to
  GSHTrans's. The batched transform's `Batch` descriptor already works this
  way for strides; the loop layer should do the same for ownership.

  **MPI, as an option.** Not in the loop: ranks share no memory, so a loop
  cannot hand iterations to another rank. MPI enters through *ownership*
  and *reduction*.
  - Ownership: a range over a stack knows its global extent and which slices
    this rank holds; `ForEach` runs over the local ones, with threads inside
    (one process per socket, as measured in `design.md` §2.3). The loop is
    written the same way on one rank or many. GSHTrans *accepts* the
    decomposition ("my slices are r0..r1 of N") from the application rather
    than imposing one: sl3dpp already distributes radial segments over
    ranks, and a second owner would fight it.
  - Reduction: local, then across ranks. `MPI_Allreduce` does not promise the
    same bits for different rank counts, so the deterministic form splits the
    global range into fixed blocks independent of ranks and threads, sums
    each block where it lives, gathers the per-block partials and combines
    them in block order: bitwise the same for any number of ranks and
    threads, for one small gather.
  - Packaging: Core stays free of MPI. The ownership descriptor and a
    pluggable cross-rank combine step live in Core with no MPI dependency; an
    optional header behind a CMake option supplies the MPI one. A build
    without MPI sees nothing of it.

## 7. A Python interface

Possible work, not started. The aim is a PyPI package exposing what no Python
library offers directly: the grid, the generalised spherical harmonics
(Wigner d-functions, Y^N_lm) and the fast transforms. Double precision only;
no templated forms.

**Shape.** An array-based toolkit, not a port of the field algebra. Fields are
numpy arrays `(nTheta, nPhi)`, which is the library's colatitude-major,
longitude-fastest layout, so they pass without copying; a batch is a 3-D
array, and `Batch(count, stride, dist)` maps onto numpy strides; coefficients
are 1-D with `(l, m)` index helpers. `n` and `lMax` are runtime arguments, as
they already are in `ForwardTransformation`. The lazy `SpinField` algebra is
not exposed: its fused loops and compile-time index checks do not survive
per-call crossing into Python, and numpy does pointwise arithmetic well. A
field is built from the grid's `theta` and `phi` arrays, not from a callback
per point (slow, and serialised by the GIL under OpenMP). Spectral operators
templated on the upper index, if exposed later, are instantiated over a small
range of `N` and chosen at runtime.

**Tools.**

- Bindings: **nanobind** (zero-copy numpy, small fast modules, C++17+). It
  maps `std::invalid_argument` to `ValueError`, and should release the GIL
  around transforms. It can target Python's stable ABI, so one wheel per
  platform covers Python ≥ 3.12.
- Build backend: **scikit-build-core**, which drives the existing CMake; the
  module is one more target behind an option that is off by default. Not
  Poetry: `poetry-core` is designed for pure-Python packages and its build-
  script hook is experimental. Develop with `pip install -e .` or `uv`.
- Layout: a `python/` directory with the binding source, `pyproject.toml` and
  pytest tests checked against the same analytic fields as the C++ tests.

**Distribution.** Publish pre-built wheels so that users need no compiler; the
sdist (which needs C++23, CMake and FFTW headers) is the fallback. Wheels are
built by **cibuildwheel** in GitHub Actions, which also runs pytest against
each one, and uploaded with PyPI trusted publishing. Order: Linux to TestPyPI,
then macOS and the real PyPI, then Windows and/or conda-forge on demand.
conda-forge already packages FFTW, OpenBLAS and modern compilers on every
platform, so it avoids most of the problems below.

**Decisions and problems, in order of weight.**

- **FFTW is GPL.** A wheel bundles it (`auditwheel` / `delocate` do this), so
  the wheel is GPL though GSHTrans is BSD. Either accept GPL for the Python
  package, or give FFTWpp a BSD backend such as pocketfft (numpy's and
  scipy's), which is a significant piece of work. Settle this before writing
  any code: it fixes the package's licence.
- **C++23 on each platform.** The `manylinux` images carry a modern GCC
  toolset. Apple clang and MSVC support for the C++23 features used here
  (multidimensional `operator[]` among them) has to be found by trying a
  build. Windows is the least certain.
- **BLAS and OpenMP in wheels.** Simplest: loop kernel only, no BLAS; macOS
  either bundles `libomp` or is built serial (`GSHTRANS_WITH_OPENMP=OFF`).
  Users who want the matrix kernel and threading build from source.
- **Fetched dependencies.** CMake fetches NumericConcepts, GaussQuad and FFTWpp
  at configure time. An sdist should build offline: pin the tags (see release
  pinning, §6) and include the sources in the sdist.

**Stages.** (1) Grid (constructor, `ForBand`, angles, weights, sizes),
forward and inverse transforms (single and batched), Wigner d-functions and
Y^N_lm, 3-j symbols. (2) Expansion class, raising and lowering, gradient,
interpolation. (3) Tensor and layered fields, only if asked for.
