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

- **Degree ceiling.** `lMax ≤ MaxSafeDegree<Real>()` — 1827 in double — set by
  underflow of the Wigner seed (`design.md` §5.1). Lifting it needs an exponent
  carried per column (Fukushima, *J. Geod.* 86, 2012). Out of scope unless a
  use above lMax ≈ 1800 in double appears; long double reaches 30 747.
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

- **The forward loop kernel's decomposition.** Each thread keeps a private
  accumulator over its colatitudes, followed by a partitioned reduction
  (`SphericalGrid.hpp`, forward loop kernel). That shape is expected to stop
  scaling well before 128 threads; the alternatives (blocks of m, the batch
  axis, a two-dimensional split) have not been measured.
- **The chunk rule on multi-cache-domain machines.** `Chunking::Count` counts
  the inverse transform's shared block as one copy, but on a multi-CCD or
  multi-socket machine it is pulled into every cache domain that reads it, so
  the rule is optimistic there.
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
