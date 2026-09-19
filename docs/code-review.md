# Code review, 2026-09-19

A file-by-file review of `develop` at `342e718` plus the `_x` → `x_` rename.
Six independent passes (core transforms; Wigner and 3-j; spin fields and
expansions; tensors; layered; build, CI, tests and docs), each reading its
files in full and confirming suspicions with small probe programs. Every
finding marked **reproduced** was re-run or re-read by a second pass before
being written down here. Reproducers are given inline as one-liners because
the probe directory was session-scoped and is gone.

What has since been fixed is listed under Status. Items that `lessons.md`
records as settled or deliberately open are not re-raised, with one exception
(C2), where the settled mechanism turns out not to work.

Severity is for *this* library's use: real rank-2 and rank-4 tensors, layered
Earth models, a 64-core target, Release builds.

## Status

Kept current as `fix-plan.md` is worked through. Anything not listed is open.

- **T1, T2, L3 — fixed** (plan Phase 1). One pair of functions in
  `Tensor/Orbits.h`, `DerivedComponent` and `DerivedCoefficient`, replaces the
  four hand-written copies; `tests/TestTensorOrbitValues.cpp` checks the
  defining relations on every component of every tensor type, spatially and
  spectrally, the layered type against the flat one, and a real elastic tensor
  applied to a strain against the double sum written out. Found on the way and
  fixed with them: `LayeredTensorField::ComponentStack<0,-1>()` and its
  expansion twin did not compile for a symmetric tensor although `Writable`
  said they would, because the slot was looked up for the component rather
  than for its representative.

- **B1, B2, B3, B4, B5, C2 — fixed** (plan Phase 2).
  - *B1.* The package config finds BLAS when the build had it. `tests/package`
    now runs the matrix kernel and a bicubic interpolant when the package says
    it has them, so a forgotten link dependency fails instead of passing.
  - *B2.* Every CI build entry, the sanitizer job and the package job install
    a BLAS and configure with `GSHTRANS_WITH_BLAS=ON`. The `undefined`
    sanitizer entry, a strict subset of `address`, is replaced by a clang one,
    advisory until it has been seen green.
  - *B3.* The option's boolean spellings are folded onto `ON`/`OFF`; anything
    else is a configure error.
  - *B4, C2.* `Details::InSerialisingRegion` sets the nested thread count to
    one inside the team, and is tested for teams of 1, 2 and max — the test
    failed with 16 threads from a team of one before the fix. The one-level
    rule now lives in `Execution::TeamSize`, replacing five copies, and is
    tested directly; the assertion that could not fail is gone. **The
    measured effect is smaller than this review's probe suggested:** OpenBLAS
    threads a GEMM only above a size threshold, and the library's default
    chunking at `lMax = 256` stays under it, so nothing changes there. At
    `lMax = 512` with a chunk of eight, a transform under
    `Execution::Sequential()` went from 4.7 cores and 340 ms to one core and
    299 ms.
  - *B5.* googletest is pinned to `v1.17.0`. The four sibling revisions are
    cache variables (`GSHTRANS_FFTWPP_TAG`, …) defaulting to `main`, for a
    release commit to set.
  - Not done, on reflection: making example 20 exit non-zero "when built `ON`
    without the kernel". That state cannot arise — `ON` without a BLAS fails
    the configure — and under `AUTO` returning zero is the documented
    behaviour.

## Summary

The numerical core is in good shape. The loop and matrix kernels, the batching
and chunking arithmetic, the parallel reduction, the row-major GEMM wrappers,
the Wigner truncation and critical-degree logic, the orbit table, the metric
signs in contraction, the `RadialMajor` transposes and the optional-dependency
guards were all checked line by line and found correct. The spin-weighted
layer's index algebra and aliasing argument hold.

The defects cluster in four places:

1. **Where the orbit table is consumed.** The derived-component logic exists
   in four hand-copied places (flat/layered × spatial/spectral) and they have
   drifted three different ways. One of them gives wrong signs on a **real
   elastic tensor** (T1). No test reads a value out of a real tensor of rank
   three or more.
2. **Where the element partition meets code written before it.** Finite
   differences, Lagrange and `Resample`'s target side ignore interfaces (L1,
   L2).
3. **Large degree and non-double precision.** The Wigner seeds underflow
   silently near lMax ≈ 1900 in double; 3-j throws on valid stretched
   triangles from l ≈ 260; `float` is advertised and effectively untested
   (W1–W3).
4. **BLAS slipped every net at once.** Missing from the installed package
   config, absent from every CI job, and its thread-serialisation mechanism
   does not do what the comment says (B1, B2, C2).

Two cross-cutting themes: preconditions that are `assert`-only vanish in the
Release builds `lessons.md` tells people to use (C1, W5, W6), and no OpenMP
region is exception-safe (C3).

## High

### T1. Real-pinned components read back with the wrong sign — reproduced
`Tensor/TensorField.h:405-407`. `turn = conjugated ? -scale : scale` is applied
to both pinned cases. Right for `Imaginary` (conj(i r) = −i r), wrong for
`Real` (conj(r) = r). It bites when the orbit table reaches a Real-pinned
member through a conjugating step, which never happens at rank 2.

    auto c = TensorField<4, ElasticSymmetry, RealTensor, Grid>(grid);
    // set c.Component<-1,0,0,1>() to 7 everywhere
    std::as_const(c).Component<1,0,0,-1>()[i,j];   // -7, should be +7

Also `<1,-1,1,-1>`; five components of `Symmetric<4>`; `<1,0,-1>` of
`Symmetric<3>`. A brute-force check of every component against the symmetry
generators and the reality condition: real elastic fails 10 permutation
relations and 4 reality relations (identically under `PointMajor`); rank 2 and
complex elastic fail none. `TensorExpr` reads operands through this accessor,
so contractions with a real elastic tensor inherit the signs.

Fix: `scale` in the Real branch (verified: zero failures across every symmetry
tried). Test gap: `ElasticTensorAppliedToAStrain` uses `ComplexTensor` and
compares an expression with itself.

### T2. `Coefficient` has the mirror-image error for Imaginary pinning — reproduced
`Expansion/TensorExpansion.h:225-234` and the copy at
`Layered/LayeredTensorField.h:402-411`. `turn = i` multiplies *outside* the
conjugation, so the result is sign·i·x where sign·conj(i·x) = −sign·i·x is
wanted. `Antisymmetric<3>` real: `Coefficient<1,0,-1>` equals
`Coefficient<-1,0,1>` where antisymmetry wants the negative; propagates
through `SurfaceGradient`. Needs rank ≥ 3 antisymmetry, so it is outside the
ranks that matter — but it is the same drift as T1 and the same fix pattern
(`std::conj(turn)`).

**Recommendation for T1, T2, L3 together:** one shared helper for "derived
component from representative", used by all four accessors, and one property
test comparing every `Component<>` and `Coefficient<>` of real `Symmetric<3>`,
`Antisymmetric<3>`, `ElasticSymmetry` against the `ComplexTensor` widening
(and the layered type against the flat one on a one-radius stack).

### T3. `Contract` and `Symmetrise` turn a partly-vanishing sum into a silent zero — reproduced
`Tensor/TensorExpr.h:465-472`, `571-575`. "Every term must be represented"
makes the whole component unrepresented if any term is, and `Materialise`
leaves unrepresented components at zero.

    Materialise(Contract<1,2>(TensorProduct(A, v)))   // A antisymmetric: all zero
    Materialise(Contract<0,1>(TensorProduct(Embed(u), v)))  // zero; Trace(Embed(T)) does not compile

`Tensor/BundleMaps.h:29-34` claims Contract handles this. Fix: fold only the
represented terms; unrepresented only if none are.

### T4. `Trace(rvalue)` returns views into a field destroyed inside `Trace` — reproduced (ASan)
`Tensor/TensorExpr.h:512-514`. `auto tr = Trace(MakeStrain(grid)); tr[1,1];`
is a heap-use-after-free: the node rightly takes ownership of the rvalue, then
`.Component<>()` hands out spans into it and the node dies. Same for
`Contract<0,1>(f()).Component<>()`, `f().Component<0,0>()`. Fix:
ref-qualify `Component() const&` and delete the `&&` overload on `TensorField`
and the nodes; have `Trace` return the rank-0 node or refuse rvalues.

### L1. Finite differences and Lagrange ignore interfaces and return NaN — reproduced
`Layered/RadialDerivatives.h:181-215`, `286-296`. Neither looks at
`HasElements()` or at repeated radii; the Fornberg and barycentric weights
divide by zero.

    auto g = RadialGrid<double>::WithElements({0.4,0.6,0.8,0.8,1.0,1.2}, {0,3,6});
    FiniteDifferenceDerivative(g, 2)   // 2 2 nan nan -3 -3, no exception
    LagrangeDerivative(g)              // all nan

FD is documented as the default, so `Gradient(e, FiniteDifferenceDerivative(mesh))`
on any model with an interface is NaN. `SplineDerivative` and
`ElementDerivative` show the right pattern. Fix: per-element stencils with the
width checked against each element, or throw.

### L2. `Resample` onto a target with interfaces overwrites the lower-side value — reproduced
`Layered/RadialResample.h:110-128`. Both copies of a coincident target radius
are answered from the upper piece: a field that is 1 below and 2 above 0.8,
resampled onto its own mesh, comes back `1 1 2 2 2 2`. The first of a
coincident pair should be answered from below when the target has elements.

### L3. Layered const `Component(i)` dangles, and drops the factor i — reproduced
`Layered/LayeredTensorField.h:228-237`. `return scale * view;` with `view` a
local lvalue; nodes hold lvalue terminals by reference. The flat type has
`std::move(view)` and a comment explaining exactly this. Reading the radial
component of a real layered vector through a const reference is a
stack-use-after-return. The Imaginary branch also returns `turn * view`
(real-valued) where the flat type returns `Complex{0, turn} * view`.

### W1. 3-j throws on valid stretched triangles at high degree — reproduced
`3j.h:279`, `337-347`. The forward pass stops at the first rise of |c1|, which
can be deep in the lower forbidden region; the backward pass then decays below
`rootSmall`, `ratio*ratio` overflows, the row normalises to zero, and the
residual check throws.

    Wigner3jMatrix<double>(260, 260, 520);   // throws; so do (264,264,527), (302,531,229)
    Wigner3jMatrix<double>(300, 300, 600);   // passes, as does (1000,1000,1999)

8 of 51 tables (l, 2l, l) for l in 250..400 throw. `Wigner3jStack(l1, l3)`
includes l2 = l1 + l3 by default. When the overflow does not happen the same
cause costs accuracy instead: (74,466,392) at m1 = −74 is 0.7 % off against an
exact value, with the residual just under tolerance. All l ≤ 7 match exact
values to 1.3e-15 and everything with l1, l3 ≤ 24 is clean. Tests stop at
l = 128.

### W2. The Wigner seeds underflow silently — reproduced
`Wigner.h:199-219`, `314-326`. The seed at l = |m| is ~sin^m θ; the degree
recursion is homogeneous, so a zero or denormal seed poisons the column.
Independent check, Σ_m |d^l_{m0}(θ)|² = 1 at θ = 0.367: exact to 1e-13 at
lMax = 1800, off by 3e-3 at 2000 and by 0.5 at 3000. In `float`: fine to 200,
O(1) wrong at 300. Nothing signals it and no limit is documented. At minimum
document or enforce a usable lMax per precision; the usual cure is a scaled
seed with a carried exponent.

### B1. The installed package is unusable when BLAS was found — reproduced
`cmake/GSHTransConfig.cmake.in:9-20`. The exported target carries `BLAS::BLAS`
but the config never calls `find_dependency(BLAS)`; a consumer fails at
generate time. Fix as for Interpolation: configure in `@GSHTRANS_BLAS_FOUND@`
and find conditionally. This affects any machine with a BLAS, including the
target.

### B2. No CI job has a BLAS
`.github/workflows/ci.yml`. Every job logs "Could NOT find BLAS". The matrix
kernel, `Blas.h`, the three `#ifdef GSHTRANS_HAVE_BLAS` test blocks and
`TuneKernel`'s BLAS branch are never compiled, run or sanitised in CI, and
example 20 returns 0 vacuously. This is why B1 went unnoticed. Add
`libopenblas-openmp-dev` and `-DGSHTRANS_WITH_BLAS=ON` to a GCC leg, the ASan
leg and the package leg.

## Medium

### C1. Grid preconditions are `assert`-only — reproduced
`SphericalGrid.h:1815-1816`. `GaussLegendreGrid<double,All,All>(2, 3, flag)`
under `-DNDEBUG` dies with `munmap_chunk(): invalid pointer`. The same
arguments with the matrix kernel throw cleanly, because `WignerMatrices`
checks. Make it a throw beside the node checks.

### C2. A team of one does not make a same-runtime BLAS serial — reproduced
`SphericalGrid.h:1346-1371`, and the matching text in `lessons.md` and
`CMakeLists.txt:148-151`. A `parallel num_threads(1)` region is *inactive*:
`omp_in_parallel()` is false inside it and a nested region gets the full team.

    outer team 1: in_parallel=0, nested team got 8 threads
    outer team 2: in_parallel=1, nested team got 1 threads

With OpenMP OpenBLAS, 20 000 dgemms of the lMax = 256, c = 8 shape: 3.4 s from
a team of one (CPU/wall 3.5), 0.85 s from an active team. So
`Execution::Sequential()` on a matrix grid still lets the BLAS take the
machine and runs slower for it. Verified fix: `omp_set_num_threads(1)` inside
the region (1.09 s, CPU/wall 1.0; it sets only that task's ICV). The laptop's
default BLAS is the pthread build, which is probably why this was not seen.

### C3. An exception inside any OpenMP region terminates the process — reproduced
`SphericalGrid.h:667-698`, `747-770`, `1156-1178`, `1285-1292`, `1366-1371`;
`Layered/RadialOperator.h:147-152`, `RadialMajor.h:218-223`,
`RadialResample.h:255-260`. The same call throws catchably when sequential.
Triggers: `bad_alloc` from 64 per-thread accumulators at large lMax;
`FFTWpp::WisdomOnly` with no wisdom (the guards at `:113`, `:163` are asserts);
any throwing user operator through the seam, e.g.
`ApplyRadially(e, FiniteDifferenceDerivative(otherGrid), Execution::Parallel(4))`.
Fix: validate what depends only on the grids before the region, and capture
an `exception_ptr` in the region to rethrow after it.

### L4. `ApplyToLines` promises alias safety it does not provide — reproduced
`Layered/RadialMajor.h:189-216`. The comment says the operator writes to a
scratch line; the code calls `op(in.Line(j), out.Line(j))`. In place on r²:
`1 -2.55 15.95 -75.7 383.5 -8033.5`. Either add the scratch line or refuse
`&in == &out` as `ApplyRadially` does.

### L5. `DifferentiationMatrix` overflows for physically scaled radii — reproduced
`Layered/RadialDerivatives.h:117-129`. Raw barycentric product: 101
Chebyshev–Lobatto nodes on [0, 6371] give all NaN; the same nodes on [0, 1]
are good to 1e-13. Scale the nodes by the interval or accumulate the ratios.

### L6. The operator's grid is never compared with the stack's
`Layered/RadialOperator.h:96-116`, `LayeredGradient.h:224-226`. All four
supplied operators expose `Radial()`; an operator built on grid A applied to a
stack on grid B with the same nR is silently wrong. A guarded
`if constexpr (requires { op.Radial(); })` keeps the seam open to bare
callables.

### E1. `Coefficient<>()` rebuilds a view on every read, and scales negatively
`Expansion/TensorExpansion.h:210-214`, `257-262`. Each call copies the grid
handle (atomic refcount), rebuilds `GSHIndices` and re-sums slot offsets.
15.7 ns per call on one thread against 2.4 for `view[l,m]`; **292 ns on
eight**. Eight rank-2 `SurfaceGradient`s at lMax = 256: 0.43 s serial, 0.82 s
on eight threads. Every spectral operator sits on this path. Hoist the block
span out of the (l, m) loops or precompute slot offsets in the constructor.

### E2. Expression nodes hold *views* by reference — reproduced (ASan)
`SpinField/SpinWeighted.h:239-241`. The documented lvalue hazard is
unavoidable for owning fields, but views are cheap handles that users name as
locals: `auto u = t.Component<1>(); auto w = t.Component<-1>(); return u*w;`
dangles. Storing views by value in `OperandStorage` costs a `shared_ptr` copy
and a span, and would also have prevented L3.

### W3. `float` is advertised and unusable beyond small degree
3-j in float returns −1.08e-2 where the truth is −2.8e-9 and *passes* its
residual check (tolerance 0.094); Wigner in float fails from lMax ≈ 250. No
test, example or benchmark instantiates `float` anywhere. Either test it and
state its limits, or drop it from the concept.

### W4. `GSHIndices::Indices()` captures `this` — reproduced (ASan, GCC 13)
`Indexing.h:155-162`. `for (auto [l,m] : w[0,0].Indices())` is a
stack-use-after-scope until GCC 15 / Clang 19 extend range-for temporaries.
Capture `mMax_` by value.

### W5. `Wigner` validates nothing; square-root tables are undersized for negative `Single` n — reproduced (ASan)
`Wigner.h:517-530`: `Wigner<double,All,All,Single>(2,2,3,1.0)` is a heap
overflow under NDEBUG. `Wigner.h:250`: tables sized from `max(mMax, nMax)`
where `abs(nMax)` is needed; `Wigner<double,All,Single,Single>(10, 0, -2, 1.0)`
reads past the end. Library callers pass mMax = lMax and are unaffected.

### W6. `FillWigner3jMatrix` / `FillCouplingMatrix` never check the buffer size
`3j.h:497-502`, `517-523`. An undersized span is a silent heap overwrite.

### B3. `GSHTRANS_WITH_BLAS` is compared as a string
`CMakeLists.txt:157-162`. `=0`, `FALSE`, `off` fall into AUTO and enable BLAS;
`=1`, `TRUE` lose the promised hard failure. Normalise or reject.

### B4. The nested-parallelism test cannot fail
`tests/TestGaussLegendreGrid.cpp:830`. `omp_get_level() > 1` is evaluated in
the outer body after the inner calls return, where it is always 1. The
"exactly one level threads" guarantee has no effective test — relevant to C2.

### B5. Tagged releases fetch moving `main`, googletest included
`CMakeLists.txt:85,91,97,126`, `tests/CMakeLists.txt:5`. Tracking `main` of
the four sibling libraries on `develop` is a recorded choice and is not
contested. But tag `1.1.0` does not identify a buildable state, and googletest
is third-party. Pin googletest to a release; pin sibling SHAs at release
commits only.

## Low

**Core.** The protected `SphericalGrid` constructor does not check
nPhi ≥ 2 lMax + 1, which both kernels rely on (derived grid, lMax = 8,
nPhi = 6: heap overflow at `:1082`); the `lMax == 0` shortcut at `:369-376`,
`:493-505` assumes the one-point grid (gate on `FieldSize() == 1`); the
reflected forward kernel at `:1486-1493` uses `w[iTheta]` where the formula
wants `w[mirror]` (harmless for Gauss–Legendre, free to fix).
`Chunking::MaximumCount = 64` is one of the power-of-two FFT strides
`ChooseThetaBlock` exists to avoid: 12–20 % slower than 62/63/65/66 at
lMax 32–63 (laptop measurement). `Tuning.h`: the schedule mirror at `:249-252`
models `copies = threads` where the matrix forward kernel uses 1;
`conclusive` is set only when Matrix wins, so a 3× loop win reads as a tie;
`TuneChunking` hard-codes complex fields. Thread-local plan caches only grow
and make `FFTWpp::CleanUp()` throw for the life of an OpenMP worker — a static
`ReleaseThreadCaches()` would help. `CoLatitudes()` and friends return views
into the shared `Impl` and dangle off a temporary grid; ref-qualify.

**Wigner / 3-j.** The seed-row binomial at `Wigner.h:314-324` is formed in
`Real` and overflows for |n| ≳ 520 in double (harmless for nMax ≤ 4, but the
API admits nMax = lMax). `Wigner() = default` leaves four members
indeterminate. `Indexing.h:101,104` assert m against `l_` rather than `mMax_`.
`atRight_` is never true at θ = π because `cos(π/2)` is 6e-17.

**Spin / expansion.** `SpinExpansion.h:86-92` and `TensorExpansion.h:93-101`
run an asserting member initialiser before their own throwing check, so
`SurfaceGradient` of an expansion with lMax == Rank aborts in Debug and throws
in Release. An expression owning a moved-in field is deep-copied each time it
is reused as an lvalue (`SpinWeighted.h:224-232` says copying is cheap). The
spectral free functions (`Raise`, `Lower`, `Evaluate`, `Interpolate`) take
only the owning `SpinExpansion`, so a tensor component cannot be raised or
evaluated. `Map` with a wrongly-typed callable is a hard error inside `Unary`
rather than a constraint failure. Scalars must be exactly `Real` or
`complex<Real>`: `f * 2` and `2.0 * floatField` do not compile. The bicubic
interpolant's φ wrap is one copied column under a not-a-knot spline, so the
seam cells are ~5× worse than interior; φ = NaN passes through silently. Local
`Interpolate` evaluates the field expression twice.

**Tensor.** `Permute`'s prose convention (`TensorExpr.h:231-232`, `270-271`)
is the inverse of the code for non-involutions, and the only test uses
`{2,3,0,1}`, an involution — this matters for rank-4 reorderings. `Permute`
and `GeneratedBy` accept non-bijective images; the documented
`Permute<std::array{1,0}>` deduces `array<int,2>` and fails deep inside.
`TensorField.h:368-372` implies a real tensor can live on a `NonNegative`
grid; representatives are the most *negative* N, so it cannot.
`IsTerminalTrait<TensorField>` is specialised in `TensorExpr.h` rather than
beside the class. `Materialise` defaults to `ComplexTensor` where every alias
defaults to `RealTensor`.

**Layered.** `Gradient` takes an `Execution` policy but threads only the ∂r
half; `SurfaceGradient` takes none. `ElementOf` returns −1 without elements
and the last element for any out-of-range index. A NaN radius passes
validation. `LayeredSpinExpansion::Offset` has no bounds check where
`LayeredSpinField::Offset` throws.

**Build / tests / docs.** `tests/CheckLegendre.h:35-36`:
`if (auto norm = std::abs(x) > tiny)` makes `norm` a bool, so the "relative"
check is absolute, and negative orders are silently skipped because
`std::sph_legendre` takes unsigned. `CheckAdditionTheorem.h:24` and
`CheckLegendre.h:19` seed from `std::random_device` against the policy in
`TestRandom.h`. The header self-sufficiency list omits `LayeredGradient.h`,
`LayeredTensorField.h`, `RadialMajor.h`, `RadialOperator.h` (all four compile
standalone today). `run-server-benchmark.sh`: the "suite must pass" gate at
L158-159 discards status; L109 sets `MY_PROJECT_BUILD_EXAMPLES`, which is not
this project's option. README drift: `examples/FieldExample.cpp` (L63) does
not exist; "CMake 3.20+" (L306) vs 3.24; the clang leg described as advisory;
test counts; `Batch{count, stride, dist}` (L131) does not compile. The
`undefined` sanitiser leg is a strict subset of `address`. No warning flags
anywhere; under `-Wall -Wextra -Wpedantic` the only hit is an unused typedef
at `Tuning.h:408`. The feature macros and `-fopenmp` are undocumented for the
vendoring route. `3j.h` uses `#pragma once` where everything else has guards.

**Stale comments.** Section numbers from the deleted plans survive at
`SphericalGrid.h:601,621,1317`, `Eth.h:47`, `TestConventions.cpp:3`,
`TestInterpolate.cpp:202,286`, `TestSpinField.cpp:217`,
`TestTensorIndices.cpp:403`, `TestTuning.cpp:222`,
`TestGaussLegendreGrid.cpp:1528`, `TransformBenchmark.cpp:345,881,915`.
`SphericalGrid.h:1439-1445` describes one GEMM of width 4c where the code
issues two of 2c, and the scratch at `:1338-1343` is sized for the former.
`:1903-1909` says a non-symmetric grid "does not get the halved table";
`WignerMatrices::Reflected` throws. `RadialDerivatives.h:26-29` describes
thread-local scratch that does not exist; `:35-37` and
`RadialSplineDerivative.h:9-11,41` miscount the operators.
`RadialResample.h:6` names the wrong macro. "The plan" is still cited at
`SpinFieldOverloads.h:65-67` and `LayeredSpinField.h:159`. Garbled sentences
at `SphericalGrid.h:816-819,1338`, `LayeredSpinField.h:32`,
`RadialGrid.h:60-62`, `WignerMatrices.h:67-69`, `TensorField.h:363-366`,
`SpinWeighted.h:224-226`, `Expansion/BundleMaps.h:24-26`. Doxygen block at
`TensorField.h:226-244` is attached to the wrong entity; duplicate `@brief`
at `Tensor/BundleMaps.h:90-91`.

## Coverage gaps, largest first

1. Values in a real tensor of rank ≥ 3 (T1, T2) — and the rank-4 test is
   tautological.
2. The BLAS path in CI (B2).
3. `float` anywhere; `long double` beyond Wigner and the transforms.
4. Layered grids *with interfaces* under FD, Lagrange, and as a `Resample`
   target; a const pinned component; a throwing operator under `Parallel`.
5. 3-j beyond l = 128 and stretched in l2; Wigner beyond lMax = 300.
6. `NRange = NonNegative` grids (a smoke test round-trips at 8e-14; nothing
   guards it).
7. Nested-parallelism suppression (B4).
8. `Permute` with a non-involution.

## Suggested order

1. T1 with the shared helper and the property test (also closes T2, L3).
2. B1 and B2 together, then C2 once CI can see it.
3. L1, L2, L4 — the interface handling of the layered layer.
4. C1, C3, W5, W6 — promote asserts to throws; make the regions exception-safe.
5. W1, W2 — both need a little design (stopping rule; scaled seeds), so a note
   first.
6. T3, T4, E2 — representation semantics and lifetimes.
7. E1 before any serious run on the 64-core target.
8. The comment and README sweep.
