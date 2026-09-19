# Plan: fixing the findings of the 2026-09-19 review

Companion to `code-review.md`; finding labels (T1, L2, …) are the ones used
there. The decisions table was settled on 2026-09-19; the phases are agreed in
outline and each is started only on a go-ahead. Nothing in it has been
started.

## Ground rules

1. **Test first, and watch it fail.** Every correctness fix begins with a test
   that fails on the current tree for the reason the review gives. A fix whose
   test never failed has not been shown to fix anything. Several of these bugs
   survived precisely because a test existed that could not fail (the
   tautological elastic test, the nesting test, `CheckLegendre`'s bool).
2. **One phase, one reviewable series.** Each phase below is independent
   unless a dependency is stated, leaves the tree green, and is sized to be
   read in one sitting. Commits are yours; I stop at the end of each phase.
3. **The gate for every phase:** default build + ctest; the minimal
   configuration (`-DGSHTRANS_WITH_BLAS=OFF -DGSHTRANS_WITH_INTERPOLATION=OFF`);
   `scripts/test_sanitized.sh address`; `clang-format --dry-run -Werror`; the
   `docs` target. From Phase 2 on, CI sees BLAS too.
4. **Property tests over example tests** where the property is cheap to state:
   "every component obeys every generator and the reality condition" catches
   a whole family, where "`<1,0,0,-1>` is +7" catches one member.
5. **No behaviour change rides along.** Comment and README repairs are their
   own phase at the end, so that a diff which claims to fix signs contains
   only signs.

## Decisions (settled 2026-09-19)

| # | Question | Decision |
|---|---|---|
| D1 | lMax ≳ 1900 in double? | **Out of scope.** Enforce the limit; no extended-range recursion (6b stage 2 is dropped) |
| D2 | `float` | **Keep as an option if doable.** It is, cheaply: run the recursions in `double` and narrow on store — see 6b |
| D3 | `Contract` over a partly-unrepresented sum | **Skip the absent terms**; unrepresented only if all are |
| D4 | FD on a grid with interfaces | **Per-element stencils**, throw if an element is narrower than the stencil. The radial operators are a convenience, not a core function — user codes handle this themselves — so Phase 3 makes wrong answers impossible and adds no capability |
| D5 | Views inside expression nodes | **By value** |
| D6 | Dependency pinning | **googletest to a release now; sibling SHAs at release commits** |
| D7 | `f().Component<>()` on an rvalue | **Delete the `&&` overload** |

---

## Phase 1 — derived components (T1, T2, L3)

The one that affects real work: wrong signs on a real elastic tensor.

**Cause.** "Component α from its orbit representative" is written out four
times — `TensorField::Component() const`, `LayeredTensorField::Component(i)
const`, `TensorExpansion::Coefficient`, `LayeredTensorExpansion::Coefficient` —
and each copy is wrong differently. The relation itself is small. With
`(sign, conj, constraint)` from the orbit table and `R` the stored
representative:

| constraint | value |
|---|---|
| None | `sign · R`, or `sign · conj(R)` |
| Real (R = r real) | `sign · r` — conjugation is the identity |
| Imaginary (R = i r) | `sign · i r`, or `sign · conj(i r) = −sign · i r` |

**Steps.**

1. *Tests (fail today).* New `tests/TestTensorOrbitsValues.cpp`:
   - for real `Symmetric<2>`, `Antisymmetric<2>`, `Symmetric<3>`,
     `Antisymmetric<3>`, `ElasticSymmetry`, `Symmetric<4>`, a Riemann-type
     `GeneratedBy`, both layouts: fill storage at random (seeded), read every
     representable component through a const reference, and check every
     generator relation and the reality condition
     `T^{−α} = (−1)^{N(α)} conj(T^α)`. Expected today: real elastic fails 10
     + 4, `Symmetric<4>` 12 + 10, rank 2 none.
   - the same for `Coefficient<>` on `Expand(t)`, checked against the
     coefficients of the `ComplexTensor` widening of the same field.
     Expected today: `Antisymmetric<3>` fails.
   - layered against flat on a one-radius stack, every component, through a
     const reference, under ASan. Expected today: stack-use-after-return, and
     a `RealValued`/`ComplexValued` type mismatch on Imaginary-pinned
     components (a `static_assert` on the value kind catches that half).
   - replace `ElasticTensorAppliedToAStrain`'s oracle with one computed from
     explicit component loops over a `ComplexTensor` widening, and run it on
     a `RealTensor`.
2. *One helper*, in `Tensor/Orbits.h` beside the table it interprets:
   `DerivedField<relation>(View&&)` for the spatial side (always takes the
   view by value/move, so the dangle cannot be re-created) and
   `DerivedCoefficient<relation>(stored, m, repN)` for the spectral side.
   Pure functions of the relation; no class state.
3. Route all four accessors through it. Delete the four local copies of the
   logic and the `turn` variables.
4. Correct the three comments that describe the old behaviour
   (`LayeredTensorField.h:229-231`, "as on the flat type").

**Size.** ~150 lines of tests, ~60 of helper, net deletion in the accessors.
**Risk.** Low: rank-2 behaviour is provably unchanged (the table has no
conjugated pinned members there) and the existing suite pins it.

## Phase 2 — BLAS reaches the package and CI (B1, B2, B3, B4, B5, C2)

Done early so every later phase is checked with the matrix kernel built.

1. **B1.** `GSHTransConfig.cmake.in`: `set(GSHTRANS_BLAS_FOUND @GSHTRANS_BLAS_FOUND@)`
   and a conditional `find_dependency(BLAS)`, mirroring Interpolation.
   Extend `tests/package` to call `TransformKernel::Matrix()` when
   `GSHTRANS_HAVE_BLAS` is defined, and one Interpolation symbol likewise, so
   a missing link dependency fails at link time rather than never.
2. **B3.** Normalise `GSHTRANS_WITH_BLAS`: uppercase; map `1/TRUE/YES` → `ON`,
   `0/FALSE/NO` → `OFF`; anything else is a `FATAL_ERROR` naming the three
   values.
3. **B2.** `ci.yml`: add `libopenblas-openmp-dev` and
   `-DGSHTRANS_WITH_BLAS=ON` to one GCC leg, the address-sanitiser leg and the
   package leg. `ON`, not `AUTO`, so a missing BLAS fails the job instead of
   passing vacuously. Make example 20 exit non-zero when built `ON` without
   the kernel. Replace the `undefined` leg (a strict subset of `address`)
   with ASan+UBSan on clang/libomp + BLAS.
4. **C2.** The serialisation mechanism.
   - Factor the region in `OverOrders` into one helper,
     `Details::InSerialisingRegion(threads, body)`, which opens the team and
     calls `omp_set_num_threads(1)` first thing inside it. That sets the
     nthreads ICV of each implicit task only; the caller's value is restored
     on exit (verified in the review: `omp_get_max_threads()` is unchanged
     after).
   - *Test (fails today for threads = 1):* inside the helper's body, open a
     nested `parallel` and record `omp_get_num_threads()`; require 1 for
     outer team sizes 1, 2 and max. This replaces the assertion at
     `TestGaussLegendreGrid.cpp:830` (**B4**), which cannot fail.
   - Correct the comment at `SphericalGrid.h:1346-1365`, the CMake note at
     `CMakeLists.txt:148-151`, and the `lessons.md` entry: the lesson stands
     (one runtime, one level), the stated mechanism was wrong.
   - Re-run the sequential matrix-kernel rows of the benchmark on a machine
     with OpenMP OpenBLAS, interleaved A/B from `build-bench`, and record the
     figure. Expected ~4× on those rows.
5. **B5 / D6.** Pin googletest to a release tag, and add a
   `GSHTRANS_PIN_DEPENDENCIES` list of SHAs that the release commit fills in.

**Size.** Small. **Risk.** CI churn only; the first BLAS-enabled ASan run may
surface something new, which is the point.

## Phase 3 — interfaces in the layered layer (L1, L2, L4, L5, L6)

Scoped by D4: these operators are conveniences. The aim is that none of them
can return a wrong number silently; none of them gains a feature. The seam
items (L4, L6, and Phase 5's exceptions) matter more than the operator items,
because the seam is what user codes actually run through.

1. **L1 / D4.** `FiniteDifferenceDerivative`:
   - with elements: stencils are clipped to the node's element. The loop at
     `RadialDerivatives.h:199-214` runs per element with `[lo, hi)` in place
     of `[0, nR)`. An element with fewer than `width` nodes throws at
     construction, naming the element and its size.
   - without elements: any coincident pair throws, and the message points at
     `RadialGrid::WithElements`.
   - `LagrangeDerivative` is a single global polynomial, which has no meaning
     across a discontinuity: on a grid with elements or coincident radii it
     throws and names `ElementDerivative`.
   - *Tests (fail today with NaN):* the review's mesh
     `{0.4,0.6,0.8,0.8,1.0,1.2}, {0,3,6}` with f = 2r below and −3r above →
     `2 2 2 -3 -3 -3` at orders 1 and 2; exactness to order on each element
     for polynomials; the three throws.
2. **L2.** `Resample`, target side. Rule: a target node that is the **top
   node of a target element** takes the left limit; every other node keeps
   the present right-continuous rule. For a target without elements nothing
   changes. *Test (fails today):* the step field resampled onto its own mesh
   is the identity, Linear and CubicSpline; plus a target whose interface
   sits strictly inside a source piece, which must stay continuous.
3. **L4.** `ApplyToLines`: do what the comment says. A `thread_local` scratch
   line of nR scalars, `op(in.Line(j), scratch)`, copy to `out.Line(j)`. The
   line is contiguous already so the cost is one copy of nR per line; the
   in-place repeated-application loop is what the class is for. *Test (fails
   today):* in place equals out of place. Tighten the wording at
   `RadialOperator.h:40-41` to say which side guarantees what.
4. **L5.** `DifferentiationMatrix`: map the nodes affinely to [−1, 1] before
   forming barycentric weights and divide the matrix by the half-length.
   *Test (fails today):* 101 Chebyshev–Lobatto nodes on [0, 6371] and on
   [0, 6.371e6] differentiate r³ to 1e-10 relative.
5. **L6.** In `ApplyRadially`, `ApplyToLines`, `Gradient`:
   `if constexpr (requires { op.Radial(); })` compare with the stack's radial
   grid by identity and throw, before any parallel region. Bare callables are
   untouched.

**Size.** Medium; FD is the bulk. **Risk.** Medium — FD results at element
ends change on layered grids (from NaN, so no one can be depending on them).
Plain grids are bit-identical; assert that in a test.

## Phase 4 — preconditions that survive Release (C1, W5, W6, core lows)

Mechanical. `Indexing.h` already states the policy: "the layers above validate
and throw in every build mode". This phase makes it true. Each item gets a
`EXPECT_THROW` that today either corrupts the heap under NDEBUG or aborts in
Debug, so these tests go in a file that is also built with `-DNDEBUG` in the
`build-bench` tree.

- `SphericalGrid::Impl` ctor: `lMax >= 0`, `|nMax| <= lMax`,
  `nPhi >= 2 lMax + 1`; WisdomOnly guards at `:113`, `:163`.
- `Wigner` ctor: the same checks `WignerMatrices` already makes; `{}`
  initialisers on the four members; tables sized from `abs(nMax)`
  (`Wigner.h:250`); constrain `thetaRange` to `sized_range` +
  `random_access_range`, which is what the body uses.
- `FillWigner3jMatrix`, `FillCouplingMatrix`: check the span's size; negative
  degrees throw.
- `SpinExpansionBase` and `TensorExpansion` ctors: validate through a
  `Checked(lMax)` static, as the owning `SpinExpansion` does, so the throw
  precedes the asserting member. `ContravariantDerivative.h:197` and
  `IntrinsicDerivative`: build the result at `max(operand.MaxDegree(),
  Rank + 1)` — the gradient of an lMax = 0 scalar is a valid zero, not an
  error.
- The `lMax == 0` shortcut: gate on `FieldSize() == 1`.
- Reflected forward kernel: `w[mirror]`.
- `RadialGrid`: reject NaN radii; `ElementOf` throws out of range and without
  elements; `LayeredSpinExpansion::Offset` throws as `LayeredSpinField::Offset`
  does.
- `Indexing.h:101,104`: assert against `mMax_`.
- `GSHIndices::Indices()`: capture `mMax_` and `l`-range by value (**W4**);
  test under ASan on a temporary.

**Size.** Small, wide. **Risk.** Low.

## Phase 5 — exceptions and OpenMP regions (C3)

1. A small utility in `Utility.h`:

       class ExceptionCapture {
        public:
         template <typename F> void Run(F&& f) noexcept;  // try/catch, first one wins
         bool Failed() const noexcept;                    // lets loops skip remaining work
         void Rethrow();                                  // after the region
       };

   First-writer-wins under `omp critical`; `Failed()` is an atomic flag read,
   so the uncontended cost is one load per iteration of the *outer* loop
   only.
2. Apply it to the eight regions (five in `SphericalGrid.h`, three in
   `Layered/`). Bodies are wrapped at the granularity of one order / one
   line / one chunk, not per sample.
3. Hoist what depends only on the arguments out of the regions: shape checks,
   L6's grid check, Akima's minimum element size in `Resample`.
4. *Tests (terminate today):* a user operator that throws on line 7 under
   `Execution::Parallel(4)` is catchable and carries its message;
   `WisdomOnly` with no wisdom under Parallel throws `runtime_error`. Run the
   first under TSan as well.

**Size.** Small. **Risk.** Low; one benchmark pass to confirm no measurable
cost on the transform rows.

## Phase 6 — large degree (W1, W2, W3) — *each needs a short note agreed first*

### 6a. 3-j: stretched triangles (W1)

**What is wrong.** Two separate things. (i) The join always scales *up*, so
`ratio * ratio` overflows when forward and backward runs differ by more than
~1e154 — that is the throw. (ii) The stopping rule is "first rise of |c1|",
inherited from SLATEC, where it is a proxy for "the values have stopped
growing". For near-stretched rows the proxy fails: at (260,260,520), m1 = 259
the forward run stops at k = 491 of 521 while still growing, and the backward
run then descends 30 steps in its *unstable* direction. When (i) does not
overflow, (ii) still costs accuracy: 0.7 % at (74,466,392).

**Proposal A (small, try first).**
- Stop on the quantity itself: run forward while |g[k]| grows, stop at the
  first local maximum (which is the edge of the classical region, or the
  single peak of a stretched row). Guard against a spurious first-step
  maximum with a two-step confirmation.
- Overflow-safe join: form `hypot(ratio·√sumForward, √sumBackward)` rather
  than `ratio²·sumForward + sumBackward`, and scale the side that is scaled
  *down*.
- Keep the phase tracking exactly as it is (`lessons.md`: do not relearn).

**Proposal B (if A does not close it).** Luscombe & Luban's scheme
(Phys. Rev. E 57, 7274): two-term *ratio* recursion in the non-classical
regions, three-term in the classical one, joined at the turning points.
Ratios are O(1), so nothing over- or underflows by construction. This is a
fourth algorithm, so it has to earn it on the table below.

**Acceptance table — built before either proposal, and failing today:**
- no throw on any (l, 2l, l), (l, 2l−1, l), (l, 2l−2, l) for l ≤ 1000, nor
  on the (l1, l1+l3, l3) grid 150..400;
- `RacahReference.h` is *exact to rounding on short sums*, which is exactly
  the stretched region: agreement to 100 ε wherever `RacahSumLength ≤ 3`,
  over those same sweeps;
- long double against double to 1e-13 relative-to-row-scale on 300 seeded
  random triples with l ≤ 500 (9 exceed 2000 ε today);
- cyclic-permutation invariance at the existing high-degree points, which is
  the only check that sees a negated row;
- cost within 1.2× of today on the existing benchmark rows.

**W3, 3-j part.** Once A or B is in, run the table in `float`. Then either
tighten `ResidualTolerance<float>` until garbage fails it, or state the
supported range. `ResidualTolerance` at 1000·(n+1)·ε is 0.094 for a long
float row, which is not a check.

### 6b. Wigner: seed underflow (W2), and what `float` means (W3, D2)

**The limit has a closed form.** The seed at l = |m| is ~ (sin θ)^m, and the
column recovers to O(1) only once l sin θ ≳ m. A seed that underflows *and*
matters therefore needs m |ln sin θ| > |ln min| with m < lMax sin θ, i.e.
lMax · s|ln s| > |ln min|, and s|ln s| peaks at 1/e. So the unscaled
recursion is safe iff

    lMax < e · |ln min|  ≈ 1928 (double), 237 (float), 30 870 (long double).

This matches the measurements: double exact at 1800 and wrong at 2000; float
fine at 200 and wrong at 250; worst at θ ≈ asin(1/e) = 0.368. The estimate
drops an algebraic prefactor ((π m)^{-1/4} at n = 0) and is for n = 0, so the
enforced limit should carry a margin — a factor ε⁻¹ on the seed gives 1830,
194 and 30 700 — and be confirmed by the unitarity test below at the limit
itself, for each |n| ≤ 4, before it is written into the code.

**Stage 1 (do regardless).** A `constexpr MaxSafeDegree<Real>()` from that
formula; `Wigner`, `WignerMatrices` and `SphericalGrid::Impl` throw above it,
with a message that names the precision and the limit. Document it in the
README and the reference note. *Test:* the unitarity sum Σ_m |d^l_{m0}|² at
θ = asin(1/e) is 1 to 1e-12 at the limit for double and float (slow test;
label it). This converts silent O(1) error into an exception, which is the
essential fix.

**Stage 2 is dropped (D1).** Extended-exponent recursion after Fukushima
(J. Geod. 86, 2012) is the known cure above the limit; it is recorded here so
that it need not be rediscovered, and not planned.

**`float` (D2): recursion precision is not storage precision.** Enforced
naively, the limit for a `float` grid would be lMax ≈ 190, which is not much
of an option. But `float` earns its place in the *transform* — half the table
traffic, half the field memory — and nowhere in the *recursions*, which run
once at construction. So for `Real = float`: run the Wigner recursion in
`double` into a per-thread scratch block and narrow on store; likewise the
3-j rows. The limit for a `float` grid is then the `double` one, W3's 3-j
garbage cannot occur, and `ResidualTolerance<float>` stops being a question.
Cost: construction time only, plus one block of scratch per thread.
`long double` is untouched. *Tests:* a `float` grid at lMax = 512 round-trips
to ~1e-4 through both kernels (it cannot be built correctly today); the
`float` table equals the narrowed `double` table exactly.

**Also here:** the seed-row binomial (`Wigner.h:314-324`) overflows for
|n| ≳ 520. Form it as a running product interleaved with the powers of s and
c so that no partial product leaves O(1)·(its final size); no lgamma, so the
TSan reason for the current form is respected.

## Phase 7 — representation and lifetimes in the expression layers (T3, T4, E2, Permute)

Follows D3, D5 and D7 as settled above.

1. **T3.** `ContractionNode::RepresentsFn` becomes "any term represented";
   `Sum` folds only the represented terms (a constexpr-filtered index
   sequence). The same in `SymmetriseNode`. Correct
   `Tensor/BundleMaps.h:29-34`. *Tests (fail today):* A·v for antisymmetric A
   against explicit components; `Contract<0,1>(TensorProduct(Embed(u), v))`;
   `Trace(Embed(T))` compiles and equals the tangential trace.
2. **E2 / D5.** `OperandStorage`: a view terminal is held by value; only
   owning terminals stay by reference. One trait
   (`IsViewTerminalTrait`), specialised beside each view class. This is what
   would have prevented L3. *Test:* the review's `Product(T1&)` under ASan.
   Move `IsTerminalTrait<TensorField>` from `TensorExpr.h` to
   `TensorField.h` while there.
3. **T4 / D7.** `Component() const&` with `Component() && = delete` on
   `TensorField` and on every tensor node. `Trace` returns the rank-0 node's
   component through a small owning wrapper node rather than calling
   `.Component<>()` on a local. *Tests:* `Trace(MakeStrain(grid))` under
   ASan; negative compile tests for `f().Component<0,0>()`.
4. **Permute.** Make the prose agree with the code (the code matches the
   formula at the head of the class, so the code stays). `static_assert` that
   the image is a permutation, in `Permute` and `GeneratedBy`; accept
   `std::array<int, N>` by converting to `Int`. *Test:* a 3-cycle, which
   distinguishes a permutation from its inverse; the present `{2,3,0,1}` is
   an involution and cannot.

**Size.** Medium. **Risk.** Medium: D5 changes node sizes and copy costs
slightly; D7 breaks code that was already undefined.

## Phase 8 — spectral accessor performance (E1, plus small ones)

1. **E1.** Slot offsets and the block's `GSHIndices` are computed once in the
   `TensorExpansion` constructor. A `Reader<Alphas...>()` returns a small
   value type holding `span`, indices and the relation; `reader(l, m)` is the
   hot call. `Coefficient<Alphas...>(l, m)` stays as the convenience wrapper.
   Convert `SurfaceGradient`, `IntrinsicDerivative`, `Embed`, `Tangential`
   and the layered gradient to take readers outside their (l, m) loops.
   *Measure:* the review's two benchmarks — per-call cost at 1 and 8 threads
   (15.7 / 292 ns today), and eight rank-2 `SurfaceGradient`s at lMax = 256
   (0.43 s serial, 0.82 s on 8). Target: parallel faster than serial.
   From `build-bench`, interleaved.
2. *(Optional, per D4.)* Layered `Gradient`: thread the angular half and
   `FillRadialComponent` over radii; give `SurfaceGradient` a policy.
3. Local `Interpolate`: expand from the samples already evaluated.
4. `Chunking::MaximumCount` 64 → 63, after confirming the 12–20 % on a quiet
   machine.

**Risk.** Low; results bit-identical, assert it.

## Phase 9 — API edges and tests (lows)

- Spectral free functions (`Raise`, `Lower`, `Evaluate`, `Interpolate`,
  `Coefficient`) take `SpinExpansionBase`-shaped arguments, so a tensor
  component can be raised or evaluated. Fold the three spellings of the
  real-storage rule into one function.
- `Map`: add the result-type requirement to the `requires` clause.
- Scalars: accept arithmetic types and convert to `Real` (**D2**: this is
  what makes `float` usable).
- Bicubic φ seam: three periodic ghost columns a side; reject non-finite φ.
- `Tuning.h`: schedule mirror uses the matrix kernel's `copies = 1`;
  `conclusive` means "beat the margin" for either winner; constrain
  `TuneChunking` on the grid's field type; move `TuneKernelLoopOnly` into
  `Details`.
- `ReleaseThreadCaches()`; ref-qualify `CoLatitudes()` and friends.
- Tests: `CheckLegendre`'s bool and its silently skipped negative orders;
  seed the two `random_device` checks through `TestRandom.h`; header list
  cross-checked by `file(GLOB)`; a `NonNegative`-NRange round trip; one
  `float` smoke test per layer (**D2**); `-Wall -Wextra -Wpedantic` PRIVATE
  on test targets.
- `Materialise`'s default reality: follow the operand's rather than
  `ComplexTensor`.

## Phase 10 — prose

Stale section numbers (list in the review), the matrix-kernel comments that
describe one GEMM of 4c and the scratch sized for it, the miscounted
operators, the garbled sentences, README drift, `run-server-benchmark.sh`
(the gate that does not gate, `MY_PROJECT_BUILD_EXAMPLES`), the feature macros
for the vendoring route, `#pragma once` in `3j.h`. Add to `lessons.md`: the
corrected OpenMP entry, "a relation copied four times drifts three ways", the
`MaxSafeDegree` formula, and "a test that cannot fail is worse than none".

## Phase 11 — general batching (a feature, not a fix)

**What is already there.** More than it may seem. `Batch` in `Policies.h` *is*
FFTW's advanced interface — `(count, stride, dist)`, given independently for
input and output — and it is public on
`ForwardTransformation(lMax, n, in, inBatch, out, outBatch, policy)` and its
inverse. The radial case already sits on it and is one line:
`LayeredSpinField::Batch()` returns `Batch::Contiguous(nR, FieldSize())`. The
flat tensor sits on it too, batching components that share an upper index,
with `Interleaved` for the point-major layout. So any *affine* layout works
today, including a user's `[r][point][component]` or `[point][r][component]`
array, by taking a subspan at the component's base and choosing stride and
dist.

**What is not there, and has merit.**

*(a) Known offsets.* A batch whose members sit at arbitrary offsets:
`Batch::At(offsets, stride, size)`. Uses: a subset of radii (the solid
regions only; a mask), finite-element storage with padding or duplicated
interface nodes, the several same-N components of a tensor across radii in a
user's own array. It is cheap because of how the kernels were written: they
touch a layout *only* through `Count`, `Stride`, `Offset(j, k)`, `Span` and
`Disjoint` — fourteen call sites — and always gather into chunk-local scratch,
so FFTW never sees the caller's layout and no kernel changes. `Offset` becomes
`j * stride + base(k)`. Points to settle in the note: `Batch` is three
integers passed by value today and would carry a non-owning span (FFTW's
guru interface has the same lifetime rule); `Disjoint` becomes a sort, so it
is done once at construction, which is why `At` takes the size; the branch in
`Offset` is hoisted per field, and measured on the existing benchmark rows to
confirm the affine path costs nothing.

*(b) Transforms straight into and out of radial-major.* `[(l, m)][r]` is
`Batch::Interleaved(nR, nR)` — already expressible. Giving `RadialMajor` a
`Batch()` and letting `Expand` write into one and `Evaluate` read from one
removes both transposes from the Expand → lines → Evaluate loop, which is the
production pattern. No new batching machinery at all. It is, as that file
says of itself, "a measurement rather than a principle": the scatter has
stride nR, which is the power-of-two cache hazard documented at the Fourier
stage, and the present tiled transpose was built to dodge exactly that. So:
measure direct against transpose at nR = 64, 100, 128, 200 from
`build-bench`, interleaved, and keep whichever wins — possibly both, chosen
by nR.

**What I would not build.** A two-level `(countA, distA, countB, distB)`
form: (a) subsumes it. Mixed upper indices in one call: the Wigner block and
the GEMM matrix are per-n, so there is nothing for the kernel to share, and
`TensorField` already is that loop. An FFTW-style *plan object*: it has real
attractions — validate once, make workspaces outside the parallel region
(which would help C3), a home for tuning results — but it cuts across the
recorded decision to keep policies on the grid handle with `With()`, and
(a)'s construct-time validation gets most of the benefit. Noted as the thing
to revisit if (a) grows.

**Steps.** A two-page design note for (a) and (b) first, agreed; then (b)'s
measurement, since it may be the larger win for the core application and
needs no API design; then (a). After Phases 4 and 5, which touch the same
entry points. Also fix `README.md:131`, which documents this interface with a
brace syntax that does not compile — part of why it is easy to forget it
exists.

---

## Order and dependencies

    1 ──┐
    2 ──┼── 3 ── 5           (5 uses 3's hoisted checks)
        ├── 4
        ├── 6a, 6b stage 1   (notes first)
        ├── 7                (after 1: shares the accessors)
        ├── 8                (after 1 and 7: same code)
        ├── 9 ── 10
        └── 11               (after 4 and 5; design note first)

Phases 1 and 2 first and in either order; they are small and everything else
benefits. 3, 4 and 6 are independent of one another. 7 and 8 touch the same
files as 1 and should follow it.

## What this plan does not do

It does not touch the items `lessons.md` holds open (`Scheme`,
`Chunking::Count`'s `int`, the size of `SphericalGrid.h`, the layered
duplication) — though Phase 1's helper removes the part of that duplication
that was actually wrong. It does not build Wigner extended range (dropped by D1), and
builds Luscombe–Luban (6a B) only if proposal A fails its acceptance table.
