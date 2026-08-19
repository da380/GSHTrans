# Numerical core: the step −1 plan

Companion to `field-algebra-plan.md`, which is the authority on the field layer,
and to `canonical-components.tex`, which is the authority on the mathematics.
This document covers the part of the library that was previously declared out of
scope: `GaussLegendreGrid`, `GridBase`, `Wigner`, `Indexing`, `Views`.

Nothing here is implemented yet. All decisions are taken: §6 records them, §7
records the three that are deliberately deferred to measurement, and §8 gives
the task order to work through. Where a decision went against the earlier
recommendation, §6 says so.

---

## 0. Why this exists

The field-algebra plan was written on the premise that the numerical core is
frozen. Three of its answered questions broke that premise:

- **Q5** — real-valued fields at `N ≠ 0` are to be cut *from the library*, not
  merely forbidden in the new concept. That is a deletion in
  `GaussLegendreGrid`.
- **Q3** — the grid is to become lightweight and value-semantic. That is a
  restructuring of `GaussLegendreGrid`.
- **Q4** — the field layer promises that a slice target needs no alignment
  beyond `Scalar`'s. The core does not currently honour that promise.

Reading the core closely enough to plan those three turned up further items that
are worth folding into the same pass: one silent-wrong-answer bug, one grid-size
defect that the code currently works around, and a set of performance questions
that decide what the phase-5 transform layer can be built on. Doing them as one
deliberate step is much cheaper than discovering each of them from underneath a
half-built field algebra.

**The core's algorithms are not being rewritten.** The Wigner recursions, the
Gauss–Legendre quadrature, the GSH index arithmetic and `3j.h` are all left
alone. What changes is ownership, the transform's call interface, one grid
dimension, and — in the later steps, and only if benchmarks justify it — the
*storage order* of the Wigner table.

---

## 1. Findings

Each is confirmed against the code at the cited location. Severity is about
consequence, not effort.

### Correctness

**F1 — `ForwardTransformation` accumulates into `out` without zeroing it.**
*Severity: high. Silent wrong answers.*

The coefficient loop is `*outIter++ += ...`, run once per colatitude, and
nothing initialises `out` on entry (`GaussLegendreGrid.h:196–208`). The `+=` is
correct *within* the routine — the θ-loop is the quadrature sum — but the
precondition "`out` arrives zeroed" is nowhere stated and nowhere checked.
Calling the transform twice on the same buffer doubles the answer.

Every current caller satisfies it by accident: `FFTWpp::vector` value-initialises,
and `TestRealFieldSymmetry.cpp:60–61` fills explicitly. The routine is also
inconsistent with itself — the one-point-grid early exit *assigns*
(`GaussLegendreGrid.h:141`) while the general path accumulates.

*Fix:* zero `out` on entry. If accumulation is ever wanted, it is a separate,
named entry point.

**F2 — `nPhi = 2·lMax` is one sample short of resolving `|m| ≤ lMax`.**
*Severity: high. Was masked by two workarounds. Resolved in T3.*

The longitude quadrature is the trapezoid rule on `nPhi` equally spaced points,
which is exact for `e^{i(m−m′)φ}` only when `|m − m′| < nPhi`. Analysis of a
band-`lMax` field needs `|m|, |m′| ≤ lMax`, hence `nPhi ≥ 2·lMax + 1`. The grid
uses `nPhi = max(1, 2·lMax)` (`GaussLegendreGrid.h:96`), so `m = +lMax` and
`m = −lMax` are the same discrete mode.

The code works around this rather than fixing it, in two places:

- the complex forward transform explicitly zeroes the `(lMax, lMax)` coefficient
  (`GaussLegendreGrid.h:217`);
- `RandomComplexCoefficient` zeroes the same coefficient before it is ever
  transformed (`GridBase.h:138`), and `RandomRealCoefficient` zeroes its
  imaginary part (`GridBase.h:166`, `:170`), so the round-trip tests never present the
  unresolvable mode.

That is a workaround built into the *test data generator*, which is exactly
where it is least visible. It is also the reason the field plan's test family 4
carries a caveat about phantom failures.

*Fix:* choose `nPhi` as the smallest FFT-friendly integer `≥ 2·lMax + 1`. Not
`2·lMax + 2` naively: at `lMax = 256` that is `514 = 2 × 257` with 257 prime,
which is a bad FFT length. A small helper returning the least 2-3-5-7-11-13
smooth integer `≥ 2·lMax + 1` gives `520 = 2³·5·13` there, at a cost of ~1.5% more
FFT work — and the FFT is not the expensive stage (§2). Both workarounds then
delete, and so does the field plan's caveat.

*As landed:* `FastFFTSize` in `Utility.h`, restricted to `2^a 3^b 5^c 7^d 11^e
13^f` with `e + f ≤ 1`, which is the set FFTW has codelets for. It gives
exactly the predicted `520` at `lMax = 256`, and `3, 5, 11, 13, 18, 33, 70,
130, 260` at `lMax = 1, 2, 5, 6, 8, 16, 33, 64, 128`.

**F3 — FFTW plans are executed on caller storage.**
*Severity: medium. Latent, environment-dependent.*

`Plan::Execute(in, out)` is FFTW's new-array execute; FFTWpp documents it as
valid only for buffers with "the same layout and **alignment** characteristics"
as the planning buffers (`FFTWpp/src/Plan.h:287`), and neither FFTWpp nor FFTW
checks. Plans are created on `FFTWpp::vector` work buffers, i.e. `fftw_malloc`'d
and SIMD-aligned.

- The forward transform takes the new-array path on the caller's input whenever
  it is a writable range (`GaussLegendreGrid.h:171–178`), at row offset
  `iTheta·nPhi`.
- The inverse transform takes it **unconditionally** on the caller's output
  (`GaussLegendreGrid.h:322`), which is the worse case: there is no fallback at
  all.

Even with `FFTWpp::vector` storage, the *row* offset `iTheta·nPhi·sizeof(Real)`
is a multiple of 16 but not always of 32 or 64. I have not reproduced a failure,
and FFTW's own dispatch may mask it, but the point of the fix is not the bug —
it is that the constraint propagates outward: every stride that phase 2 and the
layered layer choose would inherit an alignment obligation.

*Fix:* always copy the row through the plan's own work buffers, in both
directions. The copy is negligible beside the FFT, and it is the path a lazy
expression takes anyway. This is Q4's answer, and it is what makes the field
plan's §3.2 slice contract honest.

**F4 — `GaussLegendreGrid() = default` leaves `_lMax`, `_nMax`, `_flag`
uninitialised.** *Severity: low, easy.* Reading `MaxDegree()` on a
default-constructed grid is undefined. Once the grid is a handle (step B) the
default constructor should simply not exist.

**F5 — size checks are `assert`-only.** `assert(in.size() == FieldSize())` and
the coefficient-size asserts (`GaussLegendreGrid.h:132–136`, and the inverse's
equivalents) vanish under `NDEBUG`, so a short `out` in a release build is a
silent heap overflow. The field plan §3.10 wants this class of check in all
build modes; the core is where it belongs.

**F6 — `std::advance(wigIter, l)` should advance by `min(l, mMax)`.**
*Severity: low, latent.* In the real-valued branch with `MRange = All`, the
negative orders are skipped by advancing `l` entries
(`GaussLegendreGrid.h:205`, and the inverse's counterpart). But a degree-`l` row
holds `min(l, mMax)` negative orders, not `l`. The grid always builds its Wigner
table with `mMax = lMax`, so `l ≤ mMax` always and the bug is unreachable today.
It becomes reachable the moment a truncated-order table is used.

**F7 — `ConstGSHView::operator[](Int l, Int m)` is not `const`-qualified**
(`Views.h:84`), unlike every sibling accessor. `Wigner::operator[]` returns one
by value, so binding the result to a `const auto&` and indexing it fails to
compile. *Resolved in T1.*

**F8 — `Wigner::ComputeAll` parallelises a non-canonical loop.**
`#pragma omp parallel for` is applied to
`for (auto [n, iTheta] : Indices())` over a `cartesian_product` view
(`Wigner.h:310–312`). OpenMP's canonical loop form wants an integer induction
variable or, from 5.0, a random-access iterator loop — a structured binding over
a `cartesian_product_view` is neither obviously conforming nor portable. It is
also the nested-parallelism hazard the field plan flags for the layered layer.
*Fix:* flatten to an integer loop over `[0, nUpperIndices·nAngles)` and decode
inside.

### Interface and hygiene

**F9 — grid copies copy the whole Wigner table.** `_quad` and `_wigner` are
held by value with defaulted copy (`GaussLegendreGrid.h:343–344`). At
`lMax = 256, nMax = 2, NRange = All` the table is
`5 × 257 × ((257)² − n²) ≈ 8.5 × 10⁷` doubles ≈ **679 MB**. A single complex
field on the same grid is 2.1 MB. Copying a grid by accident is not a
performance wart, it is an out-of-memory event. This is the concrete argument
behind Q3.

**F10 — `GridBase` compiles only through a transitive include.** It uses
`std::normal_distribution`, `std::random_device`, `std::mt19937_64`,
`std::ranges::generate` and `assert` while including only `<concepts>`,
`<ranges>`, `Concepts.h` and `Indexing.h`; it works because `Indexing.h` happens
to include `<random>`, `<algorithm>` and `<cassert>`. `Indexing.h` also includes
`<iostream>` for nothing.

The same defect is in `Wigner.h`: it calls `std::make_shared` (`Wigner.h:301`)
while including no `<memory>`. It compiles today only because every current
consumer reaches it through a translation unit that has already included
`<memory>` transitively; a standalone `#include "Wigner.h"` fails on GCC 14.
Found while pinning the value convention (§6, [C8]), which is exactly the kind
of small standalone use that the library should support.

*Resolved in T1*, which also added the guard: `tests/CMakeLists.txt` generates
one translation unit per core header, including it twice, so a header that
stops standing on its own fails the build rather than waiting for the first
consumer that does not already pull the missing include. Two further defects
turned up and were fixed the same way — `Utility.h` used `std::ptrdiff_t` while
including nothing at all, and `Indexing.h` carried `<random>`, `<array>` and
`<numeric>` unused. `<random>` was the one that mattered: it was the transitive
supplier that made `GridBase`'s random generators (F12) compile.

**F11 — the coefficient-size API is triplicated.** *Resolved in T2.* `RealCoefficientSize`,
`ComplexCoefficientSize`, `CoefficientSize`, `CoefficientSizeNonNegative`, each
in `(lMax, n)` and `(n)` forms: eight functions for two concepts, with
`CoefficientSize` and `ComplexCoefficientSize` identical. After step A the real
form exists only at `n = 0`, which shrinks it further.

Collapsed to four: `CoefficientSize(lMax, n)` / `CoefficientSize(n)` for the
full storage, and `RealCoefficientSize(lMax)` / `RealCoefficientSize()` for the
reduced one. The real form takes **no upper index**, so step A's rule is
visible in the signature rather than only in a runtime check.

**F12 — random-coefficient generators live on the production grid.**
`RandomComplexCoefficient` and `RandomRealCoefficient` (`GridBase.h:124–171`)
are test scaffolding on a library class, and each seeds a fresh
`std::random_device`/`mt19937_64` per call, so no test using them is
reproducible. They belong in the test tree, taking a seed.

*Resolved in T1.* `tests/TestRandom.h` holds the seeded replacements, and the
call sites thread a seed through rather than merely accepting one: the
round-trip tests draw one seed from `TestSeed()`, report it in the failure
message, and reproduce a failing run under `GSHTRANS_TEST_SEED`. Setting that
variable to `random` restores unseeded exploration.

### Performance

These are the findings that decide what phase 5 can be built on. They are stated
with numbers so that the later steps can be argued rather than assumed. All
figures are for `lMax = 256`, `nMax = 2`, `NRange = All`, double precision.

**P1 — the transform is Legendre-stage bandwidth-bound, by a wide margin.**
At fixed `n`, the θ-loop touches the entire `(l, m)` Wigner block once per
colatitude: `257 × 66049 ≈ 1.7 × 10⁷` entries, **136 MB streamed**, for about
`3.4 × 10⁷` flops — an arithmetic intensity of ≈ **0.25 flop/byte**. The FFT
stage over the same data is `257` transforms of length 512, ≈ `6 × 10⁶` flops on
2.1 MB. So the Legendre stage costs ~6× the flops and ~65× the memory traffic of
the FFT stage. Any optimisation effort that does not address it is misdirected.

**P2 — batching is the lever, and it is nearly free.** With `k` same-spin
slices transformed together (radii, tensor components, time levels), the Wigner
table is streamed once for all `k`, so the intensity becomes ≈ `0.25 k`
flop/byte. This is the largest single factor available and it needs no change to
the Wigner data at all — only a loop restructure and a batched FFT.

**P3 — the `(m, n) → (−m, −n)` involution halves the table, and it is one
relation, not two.** Theory note §4.1 gives
`d^l_{−m,−n} = (−1)^{m−n} d^l_{mn}`. It can be realised as "store `n ≥ 0`, all
`m`" or as "store all `n`, `m ≥ 0`", and these are the *same* saving: the orbit
of `(m, n)` under the involution has size two, so the maximum reduction from it
is exactly **2×**, however it is spent. Worth recording because "use the ±n
symmetry *and* reduce to `m ≥ 0`" sounds like 4× and is not.

**P4 — the θ ↔ π−θ reflection is an independent 2×.**
`d^l_{mn}(π−θ) = (−1)^{l+m} d^l_{m,−n}(θ)`, and the Gauss–Legendre nodes are
symmetric about `π/2` (the quadrature is built on symmetric `x` and mapped by
`θ = acos(−x)`, `GaussLegendreGrid.h:48–49`), so the reflection maps the node set
to itself exactly. Composed with P3 the group acting on `(m, n, θ)` has order
four and generic orbits of size four: **679 MB → 170 MB**. The cost is a sign
flip and a reversed index in the inner loop, i.e. trading bandwidth for a less
regular access pattern — which is why it must be benchmarked against P2 rather
than assumed to compose with it.

**P5 — the GEMM restructure, if it is wanted, is a Wigner *layout* change.**
The Legendre stage at fixed `n` is, per order `m`,
`f^n_{lm} = Σ_θ D^{(n,m)}_{lθ} · (w_θ F_m(θ))`, which over `k` slices is a
`(nL × nTheta) × (nTheta × k)` matrix–matrix product. To hand that to BLAS,
`D^{(n,m)}` must be contiguous in `(l, θ)` for fixed `(n, m)`. The present table
is contiguous in `(l, m)` for fixed `(n, θ)` — the transpose of what is wanted.
`Wigner`'s existing `Storage` tag orders the `(n, θ)` axes only; this is a
finer, different axis order.

**P6 — the batched FFT can produce the layout the batched Legendre stage
wants, for free.** FFTW's advanced interface (`fftw_plan_many_dft_r2c` /
`_dft`) takes arbitrary input and output strides, so a single plan over
`howmany = nTheta·k` rows can write its output directly in `[θ][m][κ]` order.
The transpose that a batched or GEMM Legendre stage needs then costs nothing.

**P7 — plans and work buffers are created per call.** Two `FFTWpp::vector`
allocations and one plan per `ForwardTransformation`
(`GaussLegendreGrid.h:152–164`), and the same in the inverse. FFTW's planner is
not thread-safe. Separately: the constructor generates wisdom and then sets
`_flag = WisdomOnly` (`GaussLegendreGrid.h:70`), so any plan shape the
constructor did not anticipate — every batched shape, in particular — will fail
to plan rather than fall back. Whatever the plan cache looks like, that flag
policy has to change with it.

---

## 2. What the findings imply, in one paragraph

The core is sound where it is hard — the recursions, the quadrature, the index
arithmetic — and weak where it is easy: ownership, preconditions, one grid
dimension, and a call interface that allocates and plans on every use. That is a
good position to be in. It means step −1 is a bounded, mostly mechanical piece
of work whose risk is concentrated in exactly one place: changing `nPhi` (F2)
changes `FieldSize()` and therefore every stored layout in the library, so it
must happen before anything is built on top, not after.

---

## 3. The steps

Steps **A–D** precede phase 1 of the field-algebra plan. Steps **E–H** are
sequenced against phase 5 and gate nothing earlier; each should be justified by
a benchmark (§5) before it is done.

### Step A — deletions

*Gates: field plan §3.2. Implements: Q5, [C7].*

Two removals of working code, each because the library should not be able to
express the thing being removed. Grouped because they are the same kind of
change and touch overlapping files.

#### A1 — real-valued transforms at `n ≠ 0`

- Constrain both transforms so that a real-valued field scalar requires `n == 0`.
  `n` is a runtime argument, so this is a `throw` in `ValidateTransformRequest`,
  not a `static_assert` (§6, [C2]).
- Reject at construction a grid whose parameters can no longer serve any
  transform: after this change an `MRange = NonNegative` grid can only do real
  transforms at `n = 0`, so `nMax != 0` on such a grid is a configuration error.
  `MRange` itself survives, documented as "real scalar grid" (§6, [C1]).
- Rewrite or retire `tests/TestRealFieldSymmetry.cpp`, which exercises exactly
  the removed path at `n = 2`.

**Keep the mathematics, delete only the storage scheme.** The relation that test
exercises — theory note `eq:basiclevel`,
`(conj f)^{−N}_{l,−m} = (−1)^{m−N} conj(f^N_{lm})` — is true, and it is the
engine of phase 4's reality reduction. What the reduced `m ≥ 0` storage at
`n ≠ 0` assumes is the *self*-relation `f^N_{l,−m} = (−1)^{m−N} conj(f^N_{lm})`,
which holds only when `f` is its own conjugate, i.e. only at `N = 0`. So: keep
the relation and reuse the test as a phase-4 oracle at `n = 0`; delete the
storage scheme.

*Risk:* low. *Effort:* small. Deletes more than it adds.

#### A2 — the `FourPi` normalisation

Delete the `Normalisation` template axis from `Wigner`, along with the `Ortho`
and `FourPi` tags and the `Normalisation` concept in `Concepts.h`. The library
offers orthonormalised functions and nothing else.

Two things make this cheaper than it looks. First, `FourPi` is a misnomer: it
does not select 4π-normalised harmonics, it selects the *unnormalised*
`d^l_{Nm}` — the raw rotation-matrix element, whose associated harmonic
`d^l_{Nm}(θ) e^{imφ}` has `∫|·|² dΩ = 4π/(2l+1)`, neither 1 nor 4π. The `Ortho` branch (`Wigner.h:474–486`)
is the one that applies `sqrt((2l+1)/(4π))`; deleting `FourPi` means deleting
the `if constexpr` and running that scaling unconditionally. Second,
`GaussLegendreGrid` already hard-codes `Ortho` (`GaussLegendreGrid.h:53`,
`:344`), so no production path changes.

The only real cost is two test headers:

- `tests/CheckAdditionTheorem.h` sums `Σ_m d^l_{Nm} d^l_{N'm}` and compares
  against `δ_{NN'}`. Under the orthonormal scaling the same sum is
  `δ_{NN'} · (2l+1)/(4π)`, so the fix is the comparison constant, one line.
  The oracle itself is unaffected — it is the strongest check on the recursion
  and it stays.
- `tests/TestWigner.cpp` names `FourPi` in two instantiations
  (`TestWigner.cpp:17`, `:27`) that are about *access patterns*, not values;
  drop the argument.

*Risk:* low. *Effort:* small.

### Step B — the grid becomes a value-semantic handle

*Gates: field plan §3.7, and therefore phase-1 step 2. Implements: Q3.*

`GaussLegendreGrid` holds a single `std::shared_ptr<const Impl>`; `Impl` owns
`_lMax`, `_nMax`, the quadrature, the Wigner table and (from step E) the plan
cache. Copy and move are pointer operations. `Identity()` returns `_impl.get()`,
which is what the field layer's grid-equality check compares. No default
constructor (F4). The two commented-out members at `GaussLegendreGrid.h:346–347`
show this was already the intended direction.

Doing the indirection *inside* the grid rather than in the field layer is the
substantive choice: otherwise every consumer — fields, expansions, the layered
wrapper, application code — decides independently how to hold a grid, and they
will not all decide the same way. It also gives step E's plan cache somewhere to
live that is shared by construction rather than by convention.

*Risk:* low, but it touches every existing use site. *Effort:* small.

### Step C — transform I/O contract

*Gates: field plan §3.2 and §3.8. Implements: Q4, F1, F3, F5.*

- Zero `out` on entry to `ForwardTransformation` (F1).
- Copy rows through the plan's own aligned work buffers in **both** directions;
  no new-array execute on caller storage (F3). The field layer's "no alignment
  beyond `Scalar`'s" contract becomes true at this point.
- Promote the size checks from `assert` to `throw std::invalid_argument`, in all
  build modes (F5), matching what `ValidateTransformRequest` already does for
  `lMax` and `n`.
- Fix F6 while in the file.

Whether the copy is later avoided — `FFTW_UNALIGNED` plans, a second plan
variant chosen on the caller's actual alignment, or folding it into step F's
pack stage where it is free — is a question for step F and does not reach the
field layer. This is the "make the layer agnostic now, optimise later" that Q4's
answer asked for.

*Risk:* low. *Effort:* small. *Measure:* the copy should be invisible against
the Legendre stage (P1); confirm rather than assume.

### Step D — grid sizing

*Gates: field plan test family 4, and every stored layout. Implements: F2.*

- `nPhi` = least FFT-friendly integer `≥ 2·lMax + 1`.
- Delete the `(lMax, lMax)` zeroing in the forward transform and the two zeroing
  hacks in the random-coefficient generators.
- Add a named constructor expressing the distinction the library currently
  cannot: the **grid** resolution versus the **band** being worked with, e.g.
  `GaussLegendreGrid::ForBand(lBand, nMax, oversampling)`. The transform already
  takes `lMax` per call, so an oversampled grid with truncated transforms works
  today; what is missing is a way to *say* it. This is the grid-level answer to
  the dealiasing question the field plan defers (§8, "Deliberately deferred"):
  a product of two band-`L` fields has band `2L`, and `Integrate(abs2(f))`
  integrates a degree-`2l` quantity, so both want an oversampled grid rather
  than a silent refinement inside the field layer.

**This is the riskiest of A–D** and the reason the whole step exists before
phase 1 rather than after it: `nPhi` changes `FieldSize()`, hence every field
buffer, every test's expected sizes, and the flat-index contract the field plan
states in §3.2. It is cheap now and expensive later.

*Risk:* medium (broad, shallow). *Effort:* small–medium. *Test:* a round-trip at
`lMax` including the `m = ±lMax` modes, which is currently impossible.

### Step E — plan and buffer cache

*Gates: nothing. Prerequisite for F, G, H. Implements: P7.*

Move plan creation and work-buffer allocation into `Impl`: plans cached by
(direction, scalar kind, length, `howmany`) under a mutex, work buffers
per-thread. Revisit the `WisdomOnly` flag policy so that an uncached plan shape
falls back to planning instead of failing.

*Risk:* low–medium (concurrency). *Effort:* medium. *Measure:* per-call
overhead at small `lMax`, where it is proportionally largest.

### Step F — the batched transform primitive

*Gates: nothing. Enables phase 5. Implements: P2, P6.*

Make the primitive "transform `k` contiguous same-spin slices", with the single
field as `k = 1` — not a scalar primitive looped from outside. Two tiers, to be
taken in order:

1. **Tier 1, no layout change.** Keep the θ-outer loop; the inner operation
   becomes an axpy of length `k` per `(l, m)`. The Wigner table is streamed once
   per *batch* rather than once per transform, which is the whole of P2's
   benefit. Batch the FFT stage with FFTW's advanced interface, using output
   strides to land the data in `[θ][m][κ]` order (P6) so the axpy is unit-stride
   in `κ`.
2. **Tier 2, GEMM.** Per-`m` matrix–matrix products against BLAS, requiring the
   Wigner layout change of P5 — hence step G.

Tier 1 is expected to capture most of the available gain for `k ≳ 8` at a small
fraction of tier 2's cost. Measure before committing to tier 2. The single-field
signature survives as a thin `k = 1` wrapper over the batched primitive
(§6, [C3]).

*Risk:* medium. *Effort:* medium (tier 1), large (tier 2).

### Step G — Wigner storage

*Gates: nothing. Benchmark-driven. Implements: P3, P4, P5.*

Three independent options, to be decided by measurement, not by inspection
(§7, [C4]):

- **symmetry reduction** (P3 and P4): up to 4× less storage and traffic, at the
  price of sign flips and reversed inner-loop indexing;
- **transform-major layout** (P5): required by tier 2 of step F, useless
  otherwise;
- **reduced precision** for the stored `d` values (§7, [C6]): halves the
  dominant traffic directly; accumulation stays in `Real`. The `d` values are
  `O(1)` and enter linearly, so the error is a relative perturbation of order
  `10⁻⁷`, but that has to be measured against the existing round-trip
  tolerances, not argued.

These interact: symmetry reduction and a transform-major layout both change the
same indexing code, and both change what the batched inner loop looks like. Do
them together or not at all.

*Risk:* medium–high (this is the only step that touches numerically delicate
code). *Effort:* large.

### Step H — threading

*Gates: nothing. Implements: P7's other half, and the field plan's layered
layer.*

Thread the batched transform over `m`-blocks or colatitudes, on per-thread
buffers from step E. Fix `Wigner::ComputeAll`'s loop form (F8) and settle the
nesting rule: one parallel region at a time, so that the layered layer's
parallel-over-slices and `ComputeAll`'s internal loop cannot nest.

### Cleanup

Not a step of its own. F7 (`Views.h` const), F10 (missing includes in
`GridBase` and `Wigner`, `Indexing.h`'s `<iostream>`) and F12 (random
generators to the test tree, with a seed) go in **T1**, where they are the
whole of the task. F11 (coefficient-size API) goes in **T2**, because step A
is what makes it collapsible. See §8.

---

## 4. Sequencing

```
A  deletions          ─┐
D  grid sizing        ─┼─→  field-algebra phase 1  ─→ phases 2–4
C  transform I/O      ─┤
B  grid handle        ─┘

E  plan cache  ─→  F  batching  ─→  G  Wigner storage  ─→  H  threading
                                   ↘  field-algebra phase 5
```

A–D are independent of each other in the sense that none needs another's result.
The order above is nevertheless the one to work in, for two reasons: **A** first
because it deletes code that **C** and **D** would otherwise have to carry
through their changes, and **D** early because it changes `FieldSize()`, so
every test written after it is written against the right sizes and every test
written before it has to be revised. **C** and **B** are pure interface work and
come last of the four.

E–H are strictly ordered and each should be preceded by the measurement that
justifies it (§5).

## 5. Benchmarks to establish before step E

There is currently no benchmark harness, and steps F–G are unarguable without
one. What is needed is small:

- forward and inverse transform, `lMax ∈ {32, 64, 128, 256}`, `n ∈ {0, 2}`,
  real and complex, `k ∈ {1, 8, 64}`;
- grid construction time and resident size as a function of `(lMax, nMax)`;
- the Legendre stage and the FFT stage timed separately, so that P1's claim is
  checked rather than inherited from this document;
- achieved bandwidth against the machine's STREAM figure, which is the number
  that says whether P2 has been realised.

The same harness answers the smaller questions raised above: whether step C's
copy is measurable (it should not be), and whether step D's larger `nPhi` costs
what §1 predicts.

---

## 6. Decisions taken

Seven questions were put here; all are answered, and an eighth ([C8]) arrived
from the field plan and is answered too. Six accepted the recommendation — of
those, [C4], [C5] and [C6] accepted a recommendation to *defer to measurement*,
which is a decision in its own right and is recorded separately in §7. One,
[C7], went against the recommendation in both halves.

### [C1] `MRange` survives, renamed in the documentation only

*Recommendation accepted.* The parameter stays. After step A,
`MRange = NonNegative` means "this grid can do exactly one thing: real
transforms at `n = 0`" — a scalar grid — and it remains a genuine 2× saving on
the Wigner table for scalar-only users, so it is not dead weight. What changes
is that step A adds a constructor-time rejection of `MRange = NonNegative` with
`nMax != 0`, and the documentation describes it as "real scalar grid" rather
than as an order range. Revisit at step G, when the Wigner storage is being
reconsidered anyway and a better factorisation of the parameter will be
obvious; deciding it now would be deciding it with the least information.

### [C2] `n` stays a runtime argument at the transform boundary

*Recommendation accepted.* The Wigner table is indexed at runtime regardless,
step F's batched primitive wants to loop over `n`, and a template parameter
would multiply the instantiations of an already heavily-templated header.
Step A's real-at-`n≠0` check is therefore a `throw` in
`ValidateTransformRequest`, not a `static_assert`.

### [C3] The batched primitive gets a thin `k = 1` wrapper

*Recommendation accepted:* one implementation, one wrapper preserving today's
single-field signature, no divergence. Two separately maintained paths would
drift.

### [C7] `FourPi` goes; `RowMajor` stays for now

*Both halves went against the recommendation*, which was to keep `FourPi` for
the addition-theorem oracle and delete `RowMajor` as an unused axis.

- **`FourPi` is deleted** (step A2). The reasoning given was that everyone
  should use fully normalised functions and a library may reasonably require
  it. That is stronger than the recommendation and it is right: the oracle
  argument for keeping `FourPi` was an argument about a *test's* convenience,
  and the test costs one changed constant (§A2). The name is also actively
  misleading — `FourPi` selects unnormalised `d`, not 4π-normalised harmonics.
- **`RowMajor` stays** until step G, where it is decided alongside the other
  storage questions ([C4], [C5]). It is unused today, so it costs an axis in
  the instantiation surface and nothing else; step G is where the whole
  storage question is opened and it is cheaper to answer it once, there, than
  to delete and possibly reinstate.

### [C8] The Wigner value convention is Dahlen & Tromp, and it is already correct

Raised in the field plan (§10 item 1) rather than here, but it is a fact about
`Wigner` so it is recorded here. **Settled and verified against the code**, not
merely expected:

> `Wigner`'s stored value at upper index `N`, degree `l`, order `m` is
> `sqrt((2l+1)/(4π)) · d^l_{Nm}(θ)`, where `d^l_{Nm} = P^N_{lm}(cos θ)` is the
> generalised Legendre function of Dahlen & Tromp (1998) eq. (C.115) — **upper
> index first**.

Checked at `lMax = 1`, `θ = 0.7` against all nine values of D&T (C.115) quoted
in the field plan. The check discriminates: the four entries with `N − m` odd
distinguish `d^l_{Nm}` from `d^l_{mN}`, which differ by `(−1)^{N−m}`, and the
code matches `d^l_{Nm}`. Two independent facts corroborate it — the seed value
`WignerMinOrder(l, n)` (`Wigner.h:52–80`) evaluates to `d^l_{n,−l}`, and the
existing `CheckLegendre` oracle pins the `N = 0` row against
`std::sph_legendre`, fixing both the orthonormal scaling and the Condon–Shortley
phase.

*Consequence:* the field plan's step 6 and its transform-based oracles are
unblocked. The one-line action is to turn the check into a permanent test
(task T1, §8), because the fact is now load-bearing for phase 4's stored
components.

---

## 7. Deferred to measurement

These three are not unanswered; they are answered by "the benchmark decides",
and none of them can be argued honestly before the §5 harness exists. Recorded
here so that steps E and F do not accidentally foreclose them — in particular,
step F must reach the Wigner data through an accessor rather than through raw
offsets.

- **[C4] Which of step G's three options, and in what order.** Symmetry
  reduction (P3, P4), transform-major layout (P5) and reduced precision ([C6])
  are independent in principle and entangled in the code. Decide after step F,
  with numbers.
- **[C5] Whether `NRange` becomes a runtime set of upper indices.** A rank-2
  tensor application wants `n ∈ {0, ±1, ±2}` and currently pays for `All`. It
  is a storage question and belongs with the other storage questions at step G.
  It interacts with P3: "store `n ≥ 0`" and "store an arbitrary set of `n`" are
  different reductions of the same axis.
- **[C6] Precision of the stored Wigner values.** Single-precision storage with
  `Real` accumulation halves the dominant traffic (P1). The working assumption
  is doubles throughout; the case for anything else has to be made by
  measurement, and *the measurement needs a tighter oracle than exists* —
  `CheckCoeff2Coeff.h` uses `50000 · ε`, which is loose enough to hide a
  regression. So this one carries a prerequisite of its own: a round-trip test
  with a defensible tolerance, before the experiment is worth running.

Also deferred to step G, from [C7]: whether `RowMajor` survives.

---

## 8. Task order

The work to start on, in order. Each task is one commit or a short series, and
each compiles and passes the suite before the next begins. T1–T5 are this
document; T6 hands over to `field-algebra-plan.md` §7.

**T1 — pin the conventions and clear the hygiene backlog.** *Done.* No
behaviour change; this is the task that makes the rest safe to review.

- `tests/CheckWignerConvention.h`: the D&T (C.115) `l = 1` table at several
  `θ`, asserting `stored == sqrt((2l+1)/(4π)) · P^N_{lm}` ([C8]). This is the
  regression that stops a future recursion change silently transposing the
  convention.
- F7 (`ConstGSHView::operator[](l, m)` missing `const`).
- F10 (`GridBase` and `Wigner` missing includes, `<iostream>` in `Indexing.h`),
  plus a compile check that each public header is self-sufficient.
- F12 (random-coefficient generators move to the test tree and take a seed).

**T2 — step A, the deletions.** *Done.* Real-valued transforms at `n ≠ 0` (A1)
and the `FourPi` normalisation (A2). `TestRealFieldSymmetry.cpp` is rewritten at
`n = 0` and kept as the phase-4 oracle. F11 (the triplicated coefficient-size
API) collapses here, because A1 is what makes it collapsible.

Two notes on how it landed. The real-at-`n≠0` check went into
`ValidateTransformRequest`, which became a template on the field scalar so that
it can see the real/complex distinction — the plan said "a throw in
`ValidateTransformRequest`" without saying how it would learn the scalar kind.
And `TestRealFieldSymmetry.cpp` was rewritten against raw coefficient buffers
rather than the `CanonicalComponentExpansion` classes it used before: those are
superseded and are deleted at step 7 of the field-algebra plan, so an oracle
that phase 4 is meant to reuse could not keep depending on them.

**T3 — step D, grid sizing.** *Done.* The riskiest of the four and the one that
must not be deferred: `nPhi` changes `FieldSize()`. Smooth-integer helper, the
`ForBand` named constructor, deletion of both `(lMax, lMax)` workarounds, and
the round-trip test at `m = ±lMax` that was impossible before.

It went more smoothly than its risk rating suggested, for a reason worth
recording: no test carried a hard-coded `FieldSize()`. They compute sizes from
the grid, so the only ones that failed were the two that asserted the aliasing
behaviour itself — including the test T2 deliberately left as a tripwire, which
fired exactly as intended. Four tests now guard the sizing, each checked
against the old `nPhi` to confirm it fails there.

`ForBand(lBand, nMax, oversampling)` takes a real oversampling factor rather
than an integer one, so the 3/2 dealiasing rule is expressible; it rounds the
resulting degree **up**, since rounding down would silently remove the headroom
the caller asked for. The grid constructor's `lMax`/`nMax` parameters widened
from `int` to `Int` so that `ForBand`'s computed degree cannot narrow.

**T4 — step C, the transform I/O contract.** F1 (zero `out`), F3 (no new-array
execute on caller storage), F5 (size checks throw in all build modes), F6.
After this the field plan's slice contract is honest.

**T5 — step B, the grid handle.** `shared_ptr<const Impl>`, `Identity()`, no
default constructor (F4). Mechanical but broad; last of the four so that it
rebases over settled code.

**T6 — hand over to phase 1.** `field-algebra-plan.md` §7 steps 1–7. Its step 6
is unblocked by T1.

**T7 — the benchmark harness of §5**, then steps E–H, each gated on the
measurement that justifies it.
