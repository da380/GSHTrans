# Numerical core: the step −1 plan

Companion to `field-algebra-plan.md`, which is the authority on the field layer,
and to `canonical-components.tex`, which is the authority on the mathematics.
This document covers the part of the library that was previously declared out of
scope: `GaussLegendreGrid`, `GridBase`, `Wigner`, `Indexing`, `Views`.

Steps A–E and H are implemented; F, F′ and G remain, and §8 records what each
task did. All decisions are taken: §6 records them, §7 records the three that
are deliberately deferred to measurement, and §8 gives the task order to work
through. Where a decision went against the earlier recommendation, §6 says so.

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
*Severity: high. Silent wrong answers. Resolved in T4.*

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
*Severity: medium. Latent, environment-dependent. Resolved in T4.*

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
uninitialised.** *Severity: low, easy. Resolved in T5:* the default
constructor is deleted, and nothing needed it. Reading `MaxDegree()` on a
default-constructed grid is undefined. Once the grid is a handle (step B) the
default constructor should simply not exist.

**F5 — size checks are `assert`-only.** *Resolved in T4.* `assert(in.size() == FieldSize())` and
the coefficient-size asserts (`GaussLegendreGrid.h:132–136`, and the inverse's
equivalents) vanish under `NDEBUG`, so a short `out` in a release build is a
silent heap overflow. The field plan §3.10 wants this class of check in all
build modes; the core is where it belongs.

**F6 — `std::advance(wigIter, l)` should advance by `min(l, mMax)`.**
*Severity: low, latent. Resolved in T4, at four sites rather than two:* the
same `l`-for-`min(l, mMax)` substitution appears in
`std::prev(outWork.end(), l)` in the complex branch of both transforms, where
it names the count of negative orders in the FFT output. Still unreachable
through `GaussLegendreGrid`, which always builds `mMax = lMax`, and therefore
still untestable from outside; it becomes reachable with a truncated-order
table. In the real-valued branch with `MRange = All`, the
negative orders are skipped by advancing `l` entries
(`GaussLegendreGrid.h:205`, and the inverse's counterpart). But a degree-`l` row
holds `min(l, mMax)` negative orders, not `l`. The grid always builds its Wigner
table with `mMax = lMax`, so `l ≤ mMax` always and the bug is unreachable today.
It becomes reachable the moment a truncated-order table is used.

**F7 — `ConstGSHView::operator[](Int l, Int m)` is not `const`-qualified**
(`Views.h:84`), unlike every sibling accessor. `Wigner::operator[]` returns one
by value, so binding the result to a `const auto&` and indexing it fails to
compile. *Resolved in T1.*

**F8 — `Wigner::ComputeAll` parallelises a non-canonical loop.** *Resolved in
T9,* flattened to an integer loop with the indices decoded inside, and with the
region suppressed when one is already open.
`#pragma omp parallel for` is applied to
`for (auto [n, iTheta] : Indices())` over a `cartesian_product` view
(`Wigner.h:310–312`). OpenMP's canonical loop form wants an integer induction
variable or, from 5.0, a random-access iterator loop — a structured binding over
a `cartesian_product_view` is neither obviously conforming nor portable. It is
also the nested-parallelism hazard the field plan flags for the layered layer.
*Fix:* flatten to an integer loop over `[0, nUpperIndices·nAngles)` and decode
inside.

### Interface and hygiene

**F9 — grid copies copy the whole Wigner table.** *Resolved in T5.* `_quad` and `_wigner` are
held by value with defaulted copy (`GaussLegendreGrid.h:343–344`). At
`lMax = 256, nMax = 2, NRange = All` the table is
`5 × 257 × ((257)² − n²) ≈ 8.5 × 10⁷` doubles ≈ **679 MB**. A single complex
field on the same grid is 2.1 MB. Copying a grid by accident is not a
performance wart, it is an out-of-memory event. This is the concrete argument
behind Q3.

*Measured, before and after.* At `lMax = 128, nMax = 2, NRange = All` one grid
costs 82 MB resident. Twenty copies of it cost **1638 MB more** with the old
value members and **nothing at all** with the handle; `sizeof(Grid)` goes from
152 bytes to 16. At `lMax = 256` the same twenty copies would have been about
14 GB.

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

**P1 — the Legendre stage dominates. *Measured: conclusion confirmed,
reasoning corrected.***
At fixed `n`, the θ-loop touches the entire `(l, m)` Wigner block once per
colatitude: `257 × 66049 ≈ 1.7 × 10⁷` entries, **136 MB streamed**, for about
`3.4 × 10⁷` flops — an arithmetic intensity of ≈ **0.25 flop/byte**. The FFT
stage over the same data is `257` transforms of length 512, ≈ `6 × 10⁶` flops on
2.1 MB. So the Legendre stage costs ~6× the flops and ~65× the memory traffic of
the FFT stage. Any optimisation effort that does not address it is misdirected.

*Measured (T7 harness, one 8-core Zen 5 machine, `lMax = 256`, `n = 2`,
double):* forward transform 19.0 ms, of which the FFT stage is **0.43 ms** and
the Legendre stage **18.6 ms** — the Legendre stage is **44× the FFT stage** in
wall time. The conclusion holds and is if anything understated.

**The reasoning does not hold, and it matters.** "Bandwidth-bound" is wrong.
The Legendre stage achieves **7.3 GB/s** on its 136 MB of Wigner values, while
a single thread on the same machine reads the same array at **37.8 GB/s** once
the accumulator dependency chain is broken. The stage is running at **19% of
what one core can pull from memory**, so DRAM bandwidth is not the binding
constraint.

What binds it is load/store throughput. Per Wigner value the inner loop does
one 8-byte DRAM read and then, from cache, a read of the FFT bin and a
read-modify-write of the coefficient — about 40 bytes of cache traffic per 8
bytes of DRAM traffic, and four memory operations per element. Rewriting the
kernel three ways (weight hoisted out of the inner loop, real and imaginary
parts split into separate streams, raw pointers instead of iterators) moved it
between 5.4 and 8.9 GB/s: none of them is the problem, and none of them is the
fix.

**P2 — batching helps, but by ~2×, not by `k`. *Measured; this finding is a
substantial correction.*** With `k` same-spin slices transformed together
(radii, tensor components, time levels), the Wigner table is streamed once for
all `k`, so the intensity becomes ≈ `0.25 k` flop/byte. That was expected to be
the largest single factor available, needing no change to the Wigner data at
all — only a loop restructure and a batched FFT.

*Measured*, on the true `lMax = 256` geometry with the batch index innermost
and unit-stride (the layout P6 says the batched FFT can produce for free), time
**per field**:

| `k` | per field (ms) | speedup | coefficient array |
|---|---|---|---|
| 1 | 21.6 | 1.00× | 1.1 MB |
| 2 | 13.0 | 1.66× | 2.1 MB |
| 4 | 10.9 | 1.97× | 4.2 MB |
| 8 | **9.6** | **2.26×** | 8.5 MB |
| 16 | 14.0 | 1.54× | 16.9 MB |
| 32 | 16.6 | 1.30× | 33.8 MB |

Two things to take from this. **The gain saturates at about 2.3×**, not at `k`,
because the intensity argument assumes DRAM bandwidth is what binds, and P1 now
says it is not: batching removes `(k-1)/k` of the *table* reads but leaves the
per-element load/store work untouched, which is what the stage is actually
limited by. **And it reverses beyond a cache-determined optimum**: the
coefficient array is `k` times larger, and once it stops fitting in L3 — 16 MB
on this machine, crossed between `k = 8` and `k = 16` — it is streamed from
DRAM too, and the batch becomes worse than no batch. The optimum is
`k ≈ L3 / (16 · nCoefficients)` and therefore depends on the machine and on
`lMax`; it is not a constant to hard-code.

The standalone kernel these numbers come from is a faithful reproduction of the
geometry and access pattern rather than the production loop itself; it agrees
with the production Legendre time to within 15% at `k = 1` (21.6 ms against
18.6 ms), which is the check that it is measuring the right thing.

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

**P8 — the deployment target is 64–128 cores and 256+ GB, and three of the
above are laptop-shaped.** Development and every measurement so far are on an
8-core Zen 5 laptop with a 16 MB shared L3; time-critical runs are on servers
of 64–128 cores and 256 GB or more, which are dual-socket or multi-CCD or both.
Three consequences, none of which the earlier findings account for:

- **Cache is per-core, not per-machine.** EPYC gives 32 MB of L3 per CCD of
  eight cores, so a 64-core part has 256 MB aggregate and each core sees 32 MB;
  Intel's server parts reach a similar place die-wide at about 2 MB per core.
  Across the laptop and the target the **per-core L3 share is roughly constant
  at 2–4 MB** rather than shrinking with core count. Any cache-fitting formula
  must be written per core; written as `L3 / threads` it is right on the
  laptop by coincidence and wrong on the target.
- **Colatitude parallelism runs out.** Step H threads over colatitudes, of
  which there are `lMax + 1` — 129 at `lMax = 128`, 257 at 256. At 64–128 cores
  that is one to four per thread. The forward reduction, a `#pragma omp
  critical` in which each thread adds a full coefficient array into `out`, is
  serialised: 128 MB of serialised adds per transform at 128 threads and
  `lMax = 256`. Partitioning the reduction (step H names it) removes the
  serialisation but not the traffic — `threads × chunk × coeffSize` is 512 MB
  at 128 threads and a chunk of 4, against the table's 136 MB. Thread-private
  accumulators are the wrong shape above roughly sixteen threads, and which
  decomposition replaces them is [C11].
- **NUMA is unaccounted for.** The table is 648 MB at `lMax = 256, nMax = 2`
  and 5.4 GB at 512, and first touch decides which socket holds each page.
  `ComputeAll` does fill it in parallel, but it flattens `(n, iTheta)` and
  schedules statically over the whole range while a transform threads over
  `iTheta` at fixed `n`; the chunk boundaries do not line up, so locality is
  accidental and partial and a good fraction of the stream crosses the
  interconnect.

*This is why step F′ matters more than the laptop numbers suggest.* A streaming
grid has no table, hence no NUMA question at all; its ceiling is cores rather
than the shared DRAM roof that step H already reached at eight threads; and it
makes the **batch axis a usable parallel axis**, which is the decomposition
that answers the second bullet. Threading over fields needs no reduction and no
accumulator — but on the stored path each thread would stream the whole table
independently, multiplying table traffic by the thread count, whereas when each
thread generates its own values that cost is zero. Against that, the *memory*
argument for F′ weakens on a 256 GB machine: 5.4 GB is nothing there, unless
many jobs share a node.

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

**The flag policy is a hard prerequisite for F, not a tidiness item.** The
constructor generates wisdom for exactly two shapes — `(nPhi → nPhi/2+1)` and
`(nPhi → nPhi)`, both at `howmany = 1` — and then sets `WisdomOnly`. Every
batched shape is by definition absent from that wisdom, so it would fail to
plan rather than fall back. Batching cannot be attempted until this changes.

*Risk:* low–medium (concurrency). *Effort:* medium. *Measure:* per-call
overhead at small `lMax`, where it is proportionally largest.

**Done (T8).** Implemented as a **per-thread** cache of workspaces — the
aligned buffers together with the plan bound to them — rather than a shared
plan cache with per-thread buffers. The buffers have to be per-thread whatever
else happens, since they are scratch; keeping the plan with them means
execution uses the plan's own buffers rather than FFTW's new-array form, which
is the alignment obligation F3 was about, and it leaves `Impl` immutable, which
is what makes a shared grid safe to use concurrently without a lock (step B).
The price is planning once per thread per shape rather than once per shape,
which after the first is a wisdom lookup. Plan *creation* is serialised on a
process-wide mutex, because FFTW's planner is not re-entrant; execution is not
serialised, and needs no lock since each thread has its own workspace.

The `WisdomOnly` policy is gone rather than patched. Shapes are planned on
first use with the flag the caller gave, so there is nothing for the
constructor to anticipate and nothing to fall back from. Wisdom pre-generation
went with it, which is why grid construction got faster.

*Measured, interleaving the two versions to control for clock drift:*
**1.22× at `lMax = 32`**, about 1.06× at `lMax = 128`, and **nothing
measurable at `lMax = 256`** — which is what this step predicted, the overhead
being proportionally largest where the transform is smallest. The more useful
effect is on *variance*: repeated calls at `lMax = 32` spanned 0.020–0.049 ms
before and 0.016–0.022 ms after. Per-call allocation and plan lookup were not
just a cost but an erratic one.

**A caveat that outlives this step.** A first, non-interleaved comparison
appeared to show 1.33× at `lMax = 256` and similar gains everywhere. It was an
artefact: this machine's clock scales under load, and the same measurement
taken minutes apart moves by tens of percent. Comparing two versions requires
running them alternately in one session. The §5 harness now reports the best of
several windows rather than the mean of one, but that only reduces the problem.
Within-run *ratios* — the P1 stage split, the P2 batching curve — are not
affected, because their terms are measured moments apart in one process.

### Step F — the batched transform primitive

*Gates: nothing. Enables phase 5. Implements: P2, P6, [C9].*

Make the primitive "transform a batch of `k` same-spin fields", with the single
field as `k = 1` — not a scalar primitive looped from outside. The batch is
**public API** and is described by `(count, stride, dist)` rather than by
contiguity; the threading policy is an explicit per-call argument defaulting to
sequential. Both are settled in [C9], which also explains why the batch is a
first-class facility rather than an internal lever for phase 5.

**Prerequisites, and how big this actually is.** Step E's flag policy must
land first (above). Beyond that, the FFT half is smaller than this document
assumed: FFTWpp already wraps the whole advanced interface —
`Layout(rank, n, howMany, embed, stride, dist)` and
`plan_many_dft`/`_dft_r2c`/`_dft_c2r` — so [C9]'s `(count, stride, dist)`
descriptor maps onto it directly and P6's "for free" claim is confirmed rather
than hoped for. The work is concentrated in the Legendre stage and in the
validation and error reporting around the batch.

**When it must land.** Not before field-algebra phase 1, which never batches;
**before field-algebra phase 2**, which hands out tensor-component views and
whose `Layout` policy is designed around what the transform can consume ([C9]
removed the `PointMajor` repack on the strength of this step existing). A
public interface with no implementation behind it is fine while nothing
consumes it and a liability once something does.

Two tiers, to be taken in order:

1. **Tier 1, no layout change.** Keep the θ-outer loop; the inner operation
   becomes an axpy of length `k` per `(l, m)`. The Wigner table is streamed once
   per *batch* rather than once per transform, which is the whole of P2's
   benefit. Batch the FFT stage with FFTW's advanced interface, using output
   strides to land the data in `[θ][m][κ]` order (P6) so the axpy is unit-stride
   in `κ`.
2. **Tier 2, GEMM.** Per-`m` matrix–matrix products against BLAS, requiring the
   Wigner layout change of P5 — hence step G.

**What the measurements say about this step, and it is not what was expected.**
Tier 1 was expected to capture most of the available gain for `k ≳ 8`. Measured,
it captures **2.3×**, at an optimum `k` of about 8 on this machine, and gets
*worse* beyond it (P2). That is worth having and is not the transformative lever
this document assumed.

Three consequences for how tier 1 should be built, none of which were in the
original sketch:

- **`k` must be chosen, not taken.** The optimum is set by the coefficient
  array fitting in last-level cache, so it depends on the machine and on `lMax`.
  A batched call with a caller-chosen `k` far above the optimum will be slower
  than no batching at all. The primitive should therefore process a large batch
  in *chunks* of an internally chosen size, rather than handing the caller's `k`
  straight to the inner loop. This does not change [C9]'s public interface —
  `(count, stride, dist)` still describes what the caller has — but it does
  change what the implementation does with it.
- **The reason to batch is no longer principally speed.** It remains the right
  primitive for phase 5 and the layered layer, and 2.3× is real. But it should
  not be sequenced ahead of cheaper or larger wins on the strength of P2's
  original arithmetic.
- **Tier 2 is where the larger win probably is, and for a different reason than
  P5 gives.** A GEMM blocks the output in registers so that the per-element
  load/store count falls, which is exactly the constraint P1 now identifies;
  the intensity argument was never the point. That makes step G's layout change
  a prerequisite for the real gain rather than an optimisation on top of it.

**The `d` values reach the inner loop through a supplier, not through the
table** ([C10]). The consumer only ever needs *a contiguous run of `d` values
in `(l, m)` order for one `(n, iTheta)`* — that is the whole of what
`_impl->wigner[n, iTheta]` is used for today. Step F therefore writes the
batched inner loop against a supplier that hands back, per degree, a
contiguous row pointer. The stored supplier returns a pointer into the table
and costs nothing over today's code; step F′ adds a second supplier that
generates the row instead. Writing the loop this way now is what stops step F′
rewriting it.

Two constraints the seam carries, both satisfied by the present loop and
recorded so they are not broken later: degrees are visited in **ascending
contiguous order from `|n|`**, and each `(n, iTheta)` is visited once per pass.
The first forecloses any future "evaluate selected degrees only" optimisation
on the generated path, which is a price worth paying.

Measure before committing to tier 2. The single-field
signature survives as a thin `k = 1` wrapper over the batched primitive
(§6, [C3]).

*Risk:* medium. *Effort:* medium (tier 1), large (tier 2).

### Step F′ — Wigner values on the fly

*Gates: nothing. Prerequisite: step F's supplier seam. Implements: [C10].
May retire most of step G.*

A **second path**, not a replacement: the same recursion with the storage
elided. For a given `(n, iTheta)` the two-term recursion is run up the degrees
inside the transform, into per-thread scratch, and the table is never built.
Selected at grid construction by a policy value, so a streaming grid carries no
table at all; per-call selection would forfeit the memory saving, and a
template parameter would infect every downstream type ([C10]).

*How much of the recursion is fused with the consumer is a separate question
from whether the table exists*, and §8's T11 settles it by measurement rather
than in advance: the first version generates a whole `(n, iTheta)` block into
scratch, which reuses the recursion unchanged and is bit-comparable against the
stored path, and the three-rows-in-L1 form is the second rung of a ladder whose
steps are worth what the numbers say they are worth.

**It is a bandwidth-for-arithmetic trade, and the measurements say the trade
is now worth making.** Single-threaded, P1 says the stage is bound by
load/store throughput, and generating a value costs more memory operations than
loading one — so on those numbers alone this would lose. What changed is step
H: at `lMax = 128` on eight threads the stage reaches 39–40 GB/s against a
machine roof near 42, so the stored path is *at* its ceiling and further cores
buy nothing. A generated path's ceiling is set by cores, which scale. The
crossover is therefore not "very high `lMax`" but roughly *the table for one
`n` no longer fitting in last-level cache* — 17 MB per transform at
`lMax = 128` against a 16 MB L3, so it is crossed inside the range that
matters in practice.

**Memory is the unambiguous win.** 648 MB at `lMax = 256, nMax = 2`, and
5.4 GB at `lMax = 512`, against roughly 12 KB of scratch per thread. The
scratch is the same per-thread workspace step E introduced and step H's
accumulator extended.

Two things must be true for it to be competitive, and both are part of the
step rather than optimisations on top:

1. **The boundary terms stop calling `lgamma`/`exp`.** *Done, in T11's first
   commit; §8 records what it measured.* `Compute` evaluates
   `WignerMinOrder`/`WignerMaxOrder` at `m = ±l` for every `(l, θ)` — three
   `lgamma` and an `exp` apiece. Paid once inside construction that is part of
   the 0.23 s; paid per call it is of order 130k transcendental evaluations per
   transform at `lMax = 256`, comparable to the whole transform. The boundary
   value obeys a one-term recursion in `l` — the ratio to `l-1` is
   `sqrt(2l(2l-1)/((l-n)(l+n))) · sin(θ/2) cos(θ/2)` — so it rides up the
   degree loop with a square root and three multiplies. This is independent of
   everything else here, speeds up construction on the *stored* path too, and
   is more robust than the closed form rather than less, since it never forms
   `lgamma(2l+1)`. It also removes a real data race: glibc's `lgamma` writes
   the global `signgam`, and `ComputeAll` calls it from every thread, so grid
   construction is concurrent-undefined by the letter of the standard today.
   Nothing reads `signgam`, so the effect is benign and no result is wrong —
   but it is the only genuine race ThreadSanitizer finds in this library, and
   the boundary recursion is what deletes it. *Deleting it took the seed row
   as well as the degree loop; §8 says why, since this document had the seed
   row down as not worth touching.*
2. **The colatitudes are blocked.** Of the recursion's coefficients only
   `alpha` carries `cos θ`; `beta`, `gamma`, `sqrtIntInv[l ± m]` and
   `sqrtInt[l-1 ± m]` are θ-independent. Over a block of `B` colatitudes the
   four table loads and the `denom` and `f2` products amortise, leaving about
   five flops and two L1 loads per value — comparable to the stored path's
   per-value cost with no DRAM traffic at all. Unblocked, the path is roughly a
   wash. Blocking needs the FFT run for `B` rows before the Legendre stage,
   which is the batched-FFT machinery step F builds anyway.

**The numerics are not a new question.** This is the same recursion, the same
seeds, the same evaluation order and the same underflow behaviour; only the
storage differs. That is the substantive difference from step G, which
re-derives values under symmetry.

*A prediction the benchmark can settle.* Step H left the `lMax = 256` forward
shortfall — 2.7× against 4.2× at `lMax = 128` — with two unseparated
candidates: the critical-section reduction, and eight 1 MB private accumulators
competing for a 16 MB L3 *with the streaming table*. A streaming grid deletes
the table stream and leaves the accumulators owning L3. If the cache candidate
is the real one, this path fixes it and says so.

*Risk:* medium. *Effort:* medium.

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

**Step F′ changes what this step is for, and may retire most of it.** All three
options above are about the *table*. On a streaming grid there is no table to
re-lay-out and none to store in reduced precision, so two of the three become
irrelevant on that path. The third inverts in F′'s favour: P4's reflection,
composed with P3, makes the row at `π−θ` the order-reversed row at `θ` up to a
sign, so on the stored path it trades bandwidth for a less regular access
pattern — which is why P4 says it must be measured — while on the generated
path it **halves the recursion**, which is that path's dominant cost. The same
irregularity buys much more there. Hence the ordering: F′ before G, and let
F′'s numbers say what G is still for.

*Risk:* medium–high (this is the only step that touches numerically delicate
code). *Effort:* large.

### Step H — threading

*Gates: nothing. Implements: P7's other half, and the field plan's layered
layer.*

**Possibly the largest single lever, on the T7 measurements.** Reading the
Wigner-sized array scales from 11.0 GB/s on one thread to 41.9 GB/s on eight
(and *down* to 19.7 GB/s on sixteen, so simultaneous multithreading hurts and
the thread count should be cores, not hardware threads). Since the Legendre
stage is bound by per-core load/store throughput rather than by DRAM (P1), it
should scale close to linearly until it meets the DRAM roof — roughly a 4× to
5× ceiling on this machine, against tier-1 batching's 2.3×. This inverts the
document's original ordering, in which threading was the last thing to do; see
§9.

Thread the batched transform over `m`-blocks or colatitudes, on per-thread
buffers from step E. Fix `Wigner::ComputeAll`'s loop form (F8) and settle the
nesting rule: one parallel region at a time, so that the layered layer's
parallel-over-slices and `ComputeAll`'s internal loop cannot nest.

**Done (T9), over colatitudes, on the unbatched transform.** The two directions
differ and the difference decides the implementation: the inverse transform's
colatitudes write disjoint rows of the field and share only read-only input, so
they divide between threads with no reduction at all; the forward transform's
colatitudes all contribute to every coefficient, so each thread accumulates
into a private buffer and the partial sums are added at the end. The
accumulator is kept per thread between calls, for the same reason step E keeps
the work buffers.

Threading is an explicit per-call policy defaulting to sequential ([C9]).
`Execution::Parallel(threads)` names a count; `Execution::Parallel()` takes
OpenMP's, which respects `OMP_NUM_THREADS`.

**The nesting rule is enforced, not documented.** A transform asked to run in
parallel from inside an existing parallel region runs sequentially instead, and
`Wigner::ComputeAll` does the same. Exactly one level threads, whichever level
asks first.

*Measured, `n = 2`, complex, time per call:*

| `lMax` | 1 thread | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| 128, forward | 1.75 ms | 1.89× | 2.98× | **4.02×** | 3.13× |
| 128, inverse | 1.77 ms | 1.83× | 3.07× | 4.15× | **4.57×** |
| 256, forward | 12.16 ms | 1.63× | 2.51× | **2.68×** | 2.23× |
| 256, inverse | 12.54 ms | 1.69× | 2.64× | **2.98×** | 2.97× |

**This closes P1's story.** Single-threaded, the Legendre stage ran at 19% of
the read bandwidth one core can achieve, bound by load/store throughput. At
`lMax = 128` on eight threads it reaches **39–40 GB/s**, against a measured
machine roof near 42: the stage is now genuinely bandwidth-bound, which is what
P1 assumed it was to begin with. There is nothing further to win there without
reducing the traffic itself, which is what step G is for.

At `lMax = 256` the gain is smaller, 2.7–3.0×. Two candidates, not separated:
eight private accumulators of 1 MB each compete for the 16 MB L3 with the
streaming, and the forward transform's reduction is serialised through one
critical section. The forward direction is consistently the slower of the two
at eight threads, which points at the reduction; a tree reduction or a
partitioned one would be the thing to try.

Sixteen threads is not better than eight and is often worse: this machine has
eight cores and sixteen hardware threads, and for memory-bound work the second
thread on a core adds contention rather than throughput. Callers should ask for
cores.

*Against batching:* threading gives 2.7–4.2× where tier-1 batching measured
2.3× (P2), which is why the ordering question in §9 was raised. The two should
compose — threads over colatitudes, batch over fields — but they compete for
the same cache, so composing them is a measurement rather than an assumption.

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

E  plan cache  ─→  H  threading  ─→  F  batching  ─→  F′ on the fly  ─→  G  Wigner storage

           F must land before field-algebra phase 2 ([C9], step F)
           F′ needs F's supplier seam, and may retire most of G ([C10])
```

*The diagram above is the revised order.* The document originally ran
E → F → G → H; the first benchmark run said E → H → F → G, and that is what was
done (§9). F′ was added after F was planned and before it was built, on the
strength of the step-H measurements — see [C10].

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

*Built in T7 as `benchmarks/TransformBenchmark.cpp`, run with the
`TransformBenchmark` target. It is not a test and ctest does not run it.*

There is currently no benchmark harness, and steps F–G are unarguable without
one. What is needed is small:

- forward and inverse transform, `lMax ∈ {32, 64, 128, 256}`, `n ∈ {0, 2}`,
  real and complex, `k ∈ {1, 8, 64}` — the `k > 1` rows need step E's flag
  policy before they can even be planned, so the harness lands with a `k = 1`
  baseline and grows the batched rows as F does;
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

Seven questions were put here; all are answered, an eighth ([C8]) arrived from
the field plan and is answered too, and a ninth ([C10]) was raised by the author
after step H and is answered below. Six accepted the recommendation — of
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

### [C9] The batched transform is public, strided, and explicitly threaded

Raised after T3, on the observation that batching has uses beyond the layered
3D fields it was planned for — ensembles, realisations, time levels, and any
collection of fields sharing a grid and an upper index. The plan had treated
step F as an internal lever for phase 5. It is instead a public facility, which
makes three things decisions rather than implementation details.

**The batch is `(count, stride, dist)`, not "`k` contiguous slices".** This is
FFTW's advanced-interface form and costs nothing at the FFT stage, which is
where the caller's layout is actually read (P6). It covers with one primitive:

| source | `stride` | `dist` |
|---|---|---|
| radial slab, tensor components in `ComponentMajor` | `1` | `FieldSize` |
| tensor components in `PointMajor` | `nComponents` | `1` |

The second row is the consequence worth recording: the field plan states that
transforms *require* `ComponentMajor`, with an explicit repack from
`PointMajor` (§8, phase 2). A stride-aware batch removes that requirement. The
caller's stride appears in only two places — reading samples into the FFT,
which FFTW absorbs, and writing coefficients out — because the Legendre stage
works on our own buffers, whose layout we choose regardless. Whether strided
access is *fast enough* to beat repacking is a benchmark question; what changes
is that the interface no longer forces the answer.

Full FFTW `nembed` generality is rejected: the angular layout is fixed by the
grid, so most of it would be unreachable.

**Threading is an explicit per-call policy, defaulting to sequential.** The
library never creates threads unless asked. The rule, which step H must
enforce and which is now user-visible rather than internal: *exactly one level
threads*. The layered layer threads over slices and calls the transform
sequentially; a caller with one large batch asks the transform to thread; never
both, and never nested with `Wigner::ComputeAll` (F8).

**A batch shares grid, `lMax` and `n`.** The Wigner block is precisely what is
being amortised, so fields at different upper indices cannot batch together.
For a rank-2 tensor that means batching over radii *within* each
`n ∈ {0, ±1, ±2}`, not across components. Worth stating in the header, because
"batch all my fields" is the natural thing to expect and it is wrong.

*Sequencing:* the signature is fixed now; the implementation stays at step F,
after the plan cache (E) and the benchmarks (§5). Nothing in field-algebra
phase 1 calls a transform except one `k = 1` round trip, so waiting costs no
rework. What T4 must do is write the row pack and unpack as an explicit seam
rather than inline, since that is what step F widens.

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

### [C10] Wigner values may be generated on the fly, as a second path

Raised by the author while planning step F: the library computes the Wigner
values up front and reads them as needed, and computing them inside the
transform is the other option. **A second path, not a replacement**, and the
decisions it forces are three.

**It goes in step F′, after F1 and before G.** The alternative was to fold it
into step F. Against that: step F carries the deadline that field-algebra
phase 2 is blocked on it, and widening the step that has the deadline is the
wrong trade. So step F lands the batched primitive alone and unblocks phase 2;
F′ follows as its own step. What F must nevertheless do *now* is write the
inner loop against the supplier seam described in step F, because that is the
only part of F′ that is expensive to retrofit.

**It is selected at grid construction.** Not a per-call argument: the grid
would have to carry the table anyway for the calls that want it, which forfeits
the entire memory saving. Not a template parameter: it would infect every
downstream type for a choice that is about one object's storage. The public
transform signature of [C9] is untouched either way, and `Impl` stays immutable
with an empty table.

*The mechanism was named as a constructor here and is a **policy value**
instead*, settled when T11 was planned. `Chunking` had meanwhile set the
precedent — a value taken alongside the planner flag, because what it describes
is a property of the machine or of the object rather than of the call — and
`WignerValues::Stored()` / `::Generated()` reads the same way. The deciding
argument is composition: `ForBand` is itself a named constructor, so the named
form needs a second one to reach an oversampled generating grid, and a third
for whatever named constructor comes next. The substance of the decision is
unchanged, which is that the choice is made once, at construction, and is not
visible in the type.

**The reason to build it is not primarily speed at high `lMax`.** That was the
expectation; the numbers point elsewhere. The memory saving is unambiguous and
large at every size, and the speed argument is strongest not at extreme `lMax`
but wherever the stored path has hit the DRAM roof — which step H shows it has
already done, at `lMax = 128` on eight threads, inside the range the library is
actually used in. The full argument, and the two preconditions that decide
whether the trade is worth making at all, are in step F′.

---

## 7. Deferred to measurement

These four are not unanswered; they are answered by "the benchmark decides",
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

- **[C11] The forward transform's parallel decomposition above ~16 threads.**
  Raised by P8. Step H's thread-private accumulators over colatitudes are
  correct at eight threads and do not survive 64–128: the accumulators become
  the dominant memory traffic and the colatitude axis is only `lMax + 1` long.
  The candidates are m-block threading (no accumulator and no reduction, but at
  128 threads each thread's contiguous run within a degree is about four
  doubles, so it fetches half-used cache lines), threading over the batch axis
  (free on a streaming grid, ruinous on a stored one — P8), and a
  two-dimensional split over both. **This cannot be settled on the laptop**, and
  choosing on laptop numbers would repeat exactly the error §9 records the first
  benchmark run correcting. Deferred until the harness has been run on a target
  machine; T10 therefore keeps step H's decomposition and fixes only the
  serialisation.

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
rather than the `CanonicalComponentExpansion` classes it used before: those
were superseded and have since been deleted, so an oracle that phase 4 is meant
to reuse could not keep depending on them.

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

**T4 — step C, the transform I/O contract.** *Done.* F1 (zero `out`), F3 (no
new-array execute on caller storage), F5 (size checks throw in all build
modes), F6. After this the field plan's slice contract is honest. F3's fix
creates the row pack/unpack seam that step F batches, written as the explicit
`PackRow`/`UnpackRow` helpers ([C9]).

*Measured, as this step asked.* The copy is invisible, as P1 predicted:
at `lMax = 256, nMax = 2` a forward transform takes 11.9 ms with the copy
against 13.0 ms without, and an inverse 12.9 ms against 12.3 ms — differences
inside run-to-run noise in both directions. There is no case for revisiting it
before step F, where the copy folds into the pack stage and stops existing.

*On reproducing F3.* It could not be reproduced on this machine. The test
`AcceptsUnalignedCallerStorage` asserts, via `fftw_alignment_of`, that its
caller storage really is in a different FFTW alignment class from the
`fftw_malloc`'d planning buffers — so it exercises the case FFTW's contract
calls invalid — and the old code passes it anyway, under `Estimate` and under
`Measure`, at exact equality. This build tolerates the mismatch. The fix
therefore removes documented undefined behaviour rather than an observed
failure, which is what §1 claimed for it; the value of the test is that it
pins the contract and will fail on a build where the classes matter.

**T6 — hand over to phase 1.** *Done.* All seven steps of
`field-algebra-plan.md` §7; see its §11 for what phase 1 has and has not.

**T5 — step B, the grid handle.** *Done.* `shared_ptr<const Impl>`,
`Identity()`, no default constructor (F4). Expected to be "mechanical but
broad"; it turned out to be mechanical and *narrow* — no call site changed at
all. Every existing consumer either constructs a grid locally or holds
`_Grid&`, so widening what a copy costs touched nobody. The one adjustment
outside the grid was making `GridBase::ProjectFunction` `const`, without which
a value-semantic grid could not be used through a `const` handle.

`Impl` is immutable after construction and shared behind a `shared_ptr`, so
concurrent use needs no synchronisation; step E's plan cache will be the first
mutable member and will bring its own lock. `Identity()` returns
`const Impl*`, with `Impl` private — usable through `auto`, and not comparable
across grid instantiations, which is the right restriction.

**With this, steps A–D are complete and field-algebra phase 1 is unblocked.**
E–H remain, all gated on the §5 benchmark harness that does not yet exist.

**T9 — step H, threading.** *Done.* Over colatitudes, with the nesting rule
enforced rather than documented. See step H for the measurements.

**T8 — step E, the plan and buffer cache.** *Done.* See step E above for what
it changed and what it measured.

**T10 — step F, the batched transform primitive.** *In progress.* Tier 1 only:
the public `(count, stride, dist)` descriptor with its explicit threading
policy ([C9]), the batched FFT stage landing data in `[θ][m][κ]` order (P6),
the inner loop rewritten as an axpy of length `k` *against the supplier seam*
([C10]), and internal chunking of the caller's `k` rather than handing it
straight to the inner loop. The single-field signature becomes a `k = 1`
wrapper ([C3]). **Field-algebra phase 2 unblocks here.** Tier 2 (GEMM) is not
in this task and waits on step G.

Four details settled while planning it:

- **The call takes two batch descriptors, not one** — one for the field range
  and one for the coefficient range. Even when the layout *kind* is the same on
  both sides, `dist` differs, being `FieldSize` on one and `CoefficientSize` on
  the other. Equal `count` is a checked throw.
- **Caller strides enter only at the pack and unpack seam.** [C9] already
  states that the Legendre stage works on our own buffers, whose layout we
  choose. So the batched FFT writes `work.out` in `[m][κ]` order (out-stride
  `k`, out-dist `1`, free per P6) and the Legendre stage accumulates in `[j][κ]`
  order. Both are batch-innermost and unit-stride, which is the layout P2's
  2.3× was measured on, and it is obtained regardless of what the caller hands
  in — which is what `PackRow`'s comment already promised.
- **The chunk size is computed per core**, as
  `perCoreL3 / (2 · 16 · nCoefficients)`, the factor of two being headroom for
  the table stream and the FFT output. That reproduces P2's measured optimum of
  8 at `lMax = 256` on the laptop and gives 2–4 per core on the target rather
  than collapsing to 1, which is what a `L3 / threads` formula would do (P8).
  The cache figure is a settable value with a documented default, since P2 says
  the optimum is not a constant to hard-code.
- **The forward reduction becomes partitioned rather than critical.** After the
  colatitude loop each thread sums every thread's partials for its own block of
  coefficients, in parallel, with no critical section. This is step H's own
  suggested fix; it is folded in here because T10 rewrites that loop anyway, and
  because without it every batched forward measurement above a few threads
  measures the critical section instead of the transform. The decomposition
  itself is *not* changed — that is [C11], and it needs a target machine.

*Measured, and it is not measurable.* The partitioned reduction was A/B'd
against the critical section back to back on an idle laptop, three runs each,
from sources identical but for that one block. At the rows where it should
matter most — `lMax = 256`, forward, eight threads — the difference is 0.5%.
**The noise floor is larger than the effect**: the *inverse* rows and the
*single-threaded* rows, neither of which touches the reduction at all, differ by
6–13% between the two builds, which is drift and nothing else. That is what the
arithmetic predicts — at eight threads and `lMax = 256` the critical section is
eight passes over a 1 MB array against a colatitude loop streaming 136 MB, about
6% of the traffic and comfortably hidden. The change is therefore justified by
P8's 128-thread argument and by no laptop number, and it is recorded that way.

Two things follow. The measurement's value is the **noise floor**: about 10% on
this machine between back-to-back runs of identical code, so no laptop figure
below that threshold means anything, and figures from different sessions mean
less still — the same binary measured 1.75 ms and 2.23 ms at `lMax = 128` in
different power states. And it **fails to separate step H's two candidates** for
the `lMax = 256` shortfall, since the reduction's contribution is below noise.
That question stays open and needs the target machine, which is [C11]'s point.

*The harness grew a section filter* (`TransformBenchmark threading`) so that an
A/B costs one section rather than the whole run. Without it the comparison is
expensive enough to be skipped, which is how unmeasured changes get made.

**State at the pause, and what picking it up means.** Two of T10's four parts
are in, the suite passes at 76/76, and nothing is half-written: the tree builds
and every test is green at each of the commits below.

*Done.*

- The partitioned reduction, as described above.
- The supplier seam: both transform loops now take the Wigner row from `d[l]`
  once per degree instead of walking a single iterator across the whole block.
  This is a pure refactor with no behaviour change, and it is the substitution
  point step F′ needs. No new class was required — `ConstGSHView::operator[](l)`
  already *is* the row supplier, and `OffsetForDegree` is closed form, so the
  seam costs a few integer operations per degree. F′ implements the same shape
  over per-thread scratch.
- The benchmark's section filter.

*Items 1 to 3 are now done, and item 4 remains.*

1. **Done.** `Batch{count, stride, dist}` next to `Execution` in `Concepts.h`,
   with `Contiguous`, `Interleaved`, `Strided` and `One`, and `Offset(j, k)`.
   Two additions the plan did not name. `Span(size)` exists because the
   batched entry's size check is an inequality where the unbatched one is an
   equality, and both contracts are now pinned by test. `Disjoint(size)` exists
   because an overlapping descriptor would otherwise produce a wrong answer in
   silence: the transform writes every element of every field it is given.
   Both named layouts satisfy it by construction.
2. **Done, with item 3.** A batched FFT stage that nothing consumes cannot be
   tested — `Workspace` is private — so the two landed together. `Workspace`
   is keyed on `(nPhi, count, flag)` and the plan is built from an
   `FFTWpp::Ranges::Layout` on each side. P6's "for free" is confirmed by
   running code, not just by reading `plan_many_*`.
3. **Done.** The Legendre stage is an axpy of length `c` per `(l, m)` into a
   `[j][κ]` scratch, scattered to the caller's `(stride, dist)`.

   Two things the sketch above got wrong, both in the library's favour.
   Zeroing the output does not become batch-aware — it *goes away*, because
   the scatter assigns rather than accumulates, which is F1's guarantee
   obtained without touching an element the call does not own. And the inverse
   needs a **gather**, not just a scatter: it reads coefficients through the
   caller's stride, and doing that per colatitude would put the stride in the
   hot loop, so a chunk's coefficients are gathered into `[j][κ]` once. That is
   the same trade T4 measured as invisible in the other direction.

   The tests compare a batch against separate unbatched calls and demand
   **exact** equality. Tier 1 widens the inner loop without reordering any
   sum, so a difference would mean it had been restructured rather than
   widened. 89 tests pass in Debug, in Release, and under ASan and UBSan.
4. **Done.** Chunking, the `k = 1` wrapper ([C3]), and the batched rows in the
   harness.

   The chunk is a `Chunking` value set at grid construction, next to the
   planner flag, because the cache figure it needs is a property of the machine
   rather than of the call. `Automatic()` assumes a modest 8 MiB last-level
   cache; `ForCache(bytes)` takes the caller's; `Fixed(count)` defeats the
   heuristic, which is what a benchmark sweeping chunk widths needs. Not a
   preprocessor macro: the library is header-only, so a knob defined
   differently in two translation units would give `GaussLegendreGrid` two
   inline bodies and let the linker choose in silence, and it could not be
   swept without a rebuild.

   **P8's formula needed two corrections, both found by checking it against
   its own anchors.** It divides by *the threads actually running*, not by a
   fixed per-core figure: P2's optimum of eight was measured sequentially, so
   that one thread had the whole 16 MiB, and the 2 MiB per-core share of the
   same machine would predict one. And it rounds to nearest rather than
   truncating — both anchors land just below an integer, 7.94 and 1.98, so
   truncation gives seven and *one*, and the second is precisely the collapse
   to a chunk of one that P8 raised the formula to avoid. With both, the
   formula returns eight and two, which is what §8 always claimed it did.

   *Measured on the laptop, sequential, two runs back to back.* At
   `lMax = 256` both agree: the optimum is `k = 8` at **2.4–2.6×**, falling to
   1.3–1.9× at 16 and 1.05–1.20× at 32. That is P2 reproduced, including the
   reversal beyond the optimum. `ForCache(16 MiB)` — this machine's real L3 —
   returns exactly 8 there, so the heuristic lands on the measured optimum when
   given a true cache figure. `Automatic()` returns 4 and reaches 2.05–2.26×,
   which is the undershoot the conservative default is for.

   At `lMax = 128` the optimum is **not resolvable on this machine**: the peak
   swaps between 8 and 16 between the two runs, at 2.44× against 1.81× one way
   and 2.02× against 2.43× the other. Run one alone would have supported
   capping the automatic chunk near eight; run two refutes it. Recorded because
   the temptation to act on the first run was real, and because it is the
   noise floor of §8 doing exactly what that note warns about.

*Two things not to rediscover.* Zeroing the output must respect the batch
descriptor rather than filling the whole range — a `PointMajor` batch is
interleaved with components that are not part of the call, and
`std::ranges::fill(out, ...)` would destroy them. And the batched entry's size
check is a *span* check, `size >= (n-1)·stride + (count-1)·dist + 1`, not the
equality `CheckSize` makes today; the `k = 1` wrapper keeps the equality, since
its contract is that the range *is* the field.

*And one caution about measuring.* See the noise floor above. Any laptop figure
under about 10% is nothing, cross-session figures are worth less than that, and
benchmarks must run with nothing else on the machine — a first attempt at the
A/B above was thrown away because a build and a test run overlapped it.

**On ThreadSanitizer, so that it is not attempted again from scratch.** It
cannot validate this library as things stand. GCC's `libgomp` carries no TSan
annotations, so the tool cannot see OpenMP's barriers and reports every
parallel region as racing with whatever ran before it: on one run it flagged
`Wigner::ComputeAll`'s loop, which writes provably disjoint slices, and
`shared_ptr`'s *atomic* reference count. Those are false positives and there
are hundreds of them. TSan also needs `vm.mmap_rnd_bits=28` on recent kernels
or it aborts before `main`. Two things follow: the thread-safety of step H's
decomposition and of T10's batching rests on argument rather than on
measurement, and the only way to get a real answer is a Clang build against a
`libomp` compiled with `LIBOMP_TSAN_SUPPORT`. Worth doing once, before the
layered layer adds a second level of threading; not worth doing repeatedly.

**T11 — step F′, Wigner values on the fly.** *First commit done.* Two commits.
First the boundary recursion replacing `lgamma`/`exp`, which stands on its own,
applies to the stored path, removes the `signgam` race described in step F′, and
lands with a tolerance test against the present values. Then the generating
supplier behind step F's seam, selected by a policy value, validated against
the stored path and measured at `lMax ∈ {64, 128, 256}`, `k ∈ {1, 8}`, on one
and eight threads. Those numbers decide what step G is still for.

*The second commit is planned below, and it is two commits rather than one.*

*Where the first commit lands, read off the code so it need not be found
again.* `WignerDetails::WignerMinOrder` (`Wigner.h:56`) is the closed form,
three `lgamma` and an `exp`; `WignerMaxOrder` (`:86`) is it again with a sign.
There are two kinds of call site and only one of them matters:

- **`Wigner.h:378`, `:398`, `:419`, `:476`** — the `m = ±l` terms, evaluated
  once per degree inside `Compute`'s degree loops. This is the cost: at
  `lMax = 256` it is `2 × 255 × 257 ≈ 131k` evaluations per upper index, which
  is where the plan's "comparable to the whole transform" comes from.
- **`Wigner.h:354`, `:358`** — the seed row at `l = |n|`, a loop over *m* at
  fixed `l`. Only `2|n|+1` values per `(n, θ)`, so five of them for a rank-2
  application. Not worth touching *for cost*.

  **It had to be touched anyway, and the note above was wrong to imply
  otherwise.** Those two lines call `WignerMinUpperIndex`/`WignerMaxUpperIndex`,
  which are the same two closed forms under other names, so they reach `lgamma`
  too — from inside `ComputeAll`'s parallel region. Moving only the degree loop
  off the closed form would have left the `signgam` race exactly where it was,
  and this task claims to remove it. The row is generated instead from the
  exact binomial recursion `C(2l, k) = C(2l, k−1)(2l−k+1)/k`, since at `l = |n|`
  the closed form is `sqrt(C(2l, l+m)) · s^{l∓m} · c^{l±m}`. Being short is what
  makes that cheap to do, not a reason not to.

*The recursion, checked against the closed form rather than taken on trust.*
Writing `WignerMinOrder(l, n) = sqrt((2l)! / ((l−n)!(l+n)!)) · s^{l+n} · c^{l−n}`
with `s = sin(θ/2)`, `c = cos(θ/2)`, the ratio at fixed `n` is

```
d(l) / d(l−1) = sqrt( 2l(2l−1) / ((l−n)(l+n)) ) · s · c
```

which is the plan's stated ratio, confirmed. It is valid only for `l > |n|`;
at `l = |n|` one factor of the denominator vanishes and the seed is needed,
where the square root is one and the value is `s^{2|n|}` at `n = +|n|` or
`c^{2|n|}` at `n = −|n|` — one `pow` per `(n, θ)` instead of per `(l, θ)`.
Carrying `WignerMinOrder(l, ±n)` as two running values gives both boundaries,
since `WignerMaxOrder(l, n) = (−1)^{n+l} · WignerMinOrder(l, −n)`.

*Four things not to be caught by.* `Compute` walks degrees ascending
(`Wigner.h:403`), so a running value fits the existing loop with no
restructuring. The boundary terms are guarded by `l <= mMax` and are never
wanted again once `l` passes it, so the recursion stops rather than having to
be advanced unused. The `AtLeft`/`AtRight` special cases at `Wigner.h:63–71`
return exact 0/1 and must survive, since the ratio is `0 · ∞` there.
And the orthonormalisation `sqrt((2l+1)/4π)` is applied to the whole block
*afterwards* (`Wigner.h:490`), so the recursion runs on unnormalised values —
which is also the contract a generating supplier must meet in the second
commit.

*What the first commit measured.* Grid construction at `lMax = 256, nMax = 2`
on eight bound threads, the two binaries run alternately eight times: **0.166 s
before against 0.154 s after**, or about **1.08×**. Small, and *below* the ~10%
noise floor §8 records for this machine — but the ranges are disjoint over six
consecutive alternations, before never under 0.162 and after never over 0.158,
which is the only reason it is quoted at all. The size is what the arithmetic
predicts: the boundary is two values per degree out of `2l+1`, so however
expensive a `lgamma` is, it is being paid on under 1% of the table.

The accuracy is the other half. Against the closed forms, over
`lMax = 128`, `|n| ≤ 3` and six colatitudes, the worst relative discrepancy is
`2.2 × 10⁻¹³` in double and `8.7 × 10⁻¹⁷` in long double — very nearly the same
multiple of `ε`, about 1000, in both. It grows linearly in `lMax` and not with
the number of entries compared, which says the drift belongs to **the closed
form**: it exponentiates a logarithm of size `O(l log 4)` and loses bits in
proportion, where each recursion step multiplies by a factor of order one. So
the test's tolerance is written as `20 · lMax · ε`, and the new values are the
better ones. `CheckAdditionTheorem` is the end-to-end check that this did not
disturb the interior: boundary values feed the two-term recursion at every
higher degree, and it still holds `Σ_m d^l_{Nm} d^l_{N′m} = δ_{NN′}(2l+1)/4π`
to `1000 ε` over all `|N|, |N′| ≤ 40`.

*What the second commit substitutes into, read off the code.*
`ConstGSHView::operator[](l)` is already the row supplier, and both transform
loops take their row from `d[l]` once per degree (T10). The two constraints
that seam carries are stated in step F and are still satisfied: degrees
ascending and contiguous from `|n|`, and each `(n, iTheta)` visited once per
pass.

The substitution is smaller than that makes it sound, and the reason is worth
stating because it decides the shape of everything below. **`ConstGSHView`
carries no storage.** It is `(lMax, mMax, n, const Real*)` over the index
arithmetic it inherits from `GSHIndices`, and `OffsetForDegree` depends on
`mMax` and `n` but *not* on `lMax`. So a view over generated scratch has the
same type as a view over the table, and the consumer loop — both directions,
batched, threaded — does not change at all. Exactly two lines read the table,
`GaussLegendreGrid.h:217` and `:442`, and nothing else in the library, the
tests, the benchmarks or the examples does.

*The depth is a ladder, and the numbers choose the rung.* This document asked
for the deepest and should not have, because the first rung already delivers
what the step is for.

| | scratch per thread at `lMax = 256` | per value | restructure |
|---|---|---|---|
| **A — whole `(n, θ)` block** | 528 KB, L2-resident | ~11 flops, 6 L1 loads | none |
| **B — fused, a row at a time** | 12 KB, L1-resident | the same | the recursion as a three-row state machine |
| **C — fused and θ-blocked** | `B` × 12 KB | ~5 flops, 2 L1 loads | that, plus the colatitude loop and the FFT staging in blocks |

**A is what T11 lands.** `Compute` already writes a whole `(n, iTheta)` block
through a `GSHView`; pointing that view at per-thread scratch rather than at
`_data` is the entire change, so the arithmetic, the order and the rounding are
the same ones the stored path uses. That makes the acceptance test **bit-exact
equality against a stored grid** rather than agreement to round-trip tolerance,
which is a much stronger check and is available only at this rung.

Both of the things F′ exists for are complete at A. **The memory win**: no
table, so 648 MB becomes 4 MB of scratch at eight threads and `lMax = 256`, and
43 GB becomes 1 GB at 128 threads and `lMax = 1024`; `ForBand` oversampling
stops costing anything, where today 2× oversampling is about 4× the table.
**The scaling win**: each thread generates the values for its own colatitudes,
so there is no shared resource — no DRAM stream, no NUMA question, no
first-touch to arrange. That is the whole of P8's argument, and B and C add
nothing to it.

What B and C buy is single-thread speed: B saves A's L2 round-trip, and C
halves the flops by amortising `denom`, `f2` and the four `sqrtInt` loads
across a block of colatitudes. Both are worth having if the stage is close;
neither is worth its restructure if A is already ahead, and **if A on eight
threads does not beat the stored path on eight threads, C's factor of two will
not rescue it and the step should be reconsidered rather than deepened.** That
is the measurement this task exists to take.

*Why B is the risky one, so that it is not attempted casually.* The recursion's
rows are not independently addressable. `Compute` relies on the relative
alignment of rows `l`, `l-1` and `l-2` in the stored layout, writing the two
boundary terms first so that the interior iterators line up at `m = -(l-2)`;
re-deriving that as a rotating three-row state machine is where the numerically
delicate part of this step actually lives.

*The two commits.*

**First, a pure refactor of `Wigner.h`,** verified by every existing value
being bit-identical. `Compute`'s body lifts to
`WignerDetails::ComputeBlock(GSHView<Real, MRange> d, Int n, Real theta,
std::span<const Real> sqrtInt, std::span<const Real> sqrtIntInv)`, a pure
function of its arguments, and `Wigner::Compute` becomes the two-line wrapper
that builds the view over `_data`. `PreCompute` lifts to
`WignerDetails::PreComputeTables(lMax, mMax, nMax)` so that the grid and
`Wigner` share one definition of the `sqrt` tables rather than growing two.

**Then the grid side.** `Impl` gains the policy, the two `sqrt` tables —
`2·lMax + 1` entries each, the tiny tables that replace the 648 MB one — and
holds `wigner` as a `std::optional`, empty when generating. The two seam lines
branch outside the colatitude loop; `AccumulateRow` and `SynthesiseRow` take
the supplier as a parameter, which costs a parameter and no duplication, since
both are already generic over `work`. The scratch follows `Accumulator` and
`CoefficientScratch` exactly: `thread_local`, grow-only, sized
`GSHIndices<MRange>(lMax, lMax, n).Size()` for the **call's** `lMax` — so a
truncated call generates only the degrees it uses, which the stored path cannot
do.

*Tests.* Bit-exact equality of both directions, batched and unbatched, stored
against generated, at several `(lMax, n)`, including a truncated call and
`MRange = NonNegative`; the existing round trips run on a generated grid; and
construction cost and resident size, which is where the memory claim is either
true or not.

*Benchmark.* A `generated` section over `lMax ∈ {64, 128, 256}`, `k ∈ {1, 8}`,
one and eight threads, with the two paths interleaved inside one process, per
§8's noise-floor rule. The harness revision bumps, since the wrapper refuses a
binary older than the script driving it.

### What rung A measured, and it is not what the step hoped for

*Laptop, eight cores, `n = 2`, complex, time per field, two runs agreeing to
better than the noise floor.*

**The memory case is won outright.** At `lMax = 256, nMax = 2` construction
goes from **0.248 s and 648 MB** to **0.015 s and 0 MB**. Nothing here is in
doubt and nothing later can take it away.

**The speed case fails, in every configuration measured.**

| `lMax` | `k` | threads | direction | stored | generated | ratio |
|---|---|---|---|---|---|---|
| 64 | 1 | 1 | forward | 0.145 ms | 0.459 ms | 3.16× |
| 128 | 1 | 8 | forward | 0.428 ms | 1.084 ms | 2.53× |
| 256 | 1 | 1 | forward | 12.50 ms | 27.17 ms | 2.17× |
| 256 | 1 | 8 | forward | 4.26 ms | 8.08 ms | 1.89× |
| 256 | 8 | 1 | forward | 6.29 ms | 10.47 ms | 1.67× |
| 256 | 8 | 8 | inverse | 4.24 ms | 8.59 ms | 2.03× |

Two things the shape of that table says. **Batching narrows the gap** — 3.2× at
`k = 1` down to 1.4× at `k = 8` — because generation amortises over a chunk
exactly as a table stream does. And **threads do not narrow it**: the ratio at
eight threads is what it is at one. The premise that the stored path is pinned
at the DRAM roof while the generated one scales with cores is *not visible on
eight cores*, which is the only machine this has run on.

Against the falsification test §8 set — *if A on eight threads does not beat
the stored path on eight threads, C's factor of two will not rescue it* — the
answer is that it does not, and the arithmetic agrees: at `lMax = 256, k = 1`
on eight threads the generated path spends 3.8 ms more than the stored one, so
halving the recursion's flops cannot close it unless the consumer's own work is
under half a millisecond, which it is not. **B and C should not be built on
these numbers.** What is left of F′ is the memory, which is large, and the
many-core and NUMA argument of P8, which no laptop can test.

### The measurement's other finding, which is worth more than its first

*The chunk width matters more than the path does, and it wants to differ by
direction.* At `lMax = 256, k = 8` on eight threads, taking the whole batch as
one chunk rather than the heuristic's:

| direction | path | heuristic chunk | whole batch |
|---|---|---|---|
| inverse | stored | 4.24 ms | **1.93 ms** |
| inverse | generated | 8.59 ms | **2.47 ms** |
| forward | stored | 4.27 ms | **8.44 ms** |
| forward | generated | 8.08 ms | **8.59 ms** |

The inverse direction is **2.2× faster** with the whole batch, and the forward
is **2× slower** — on *both* paths, which is what makes it a fact about the
chunk rather than about generating anything. The control was run for exactly
that reason.

The cause is step H's decomposition. The forward transform gives each thread a
private accumulator of `chunk × coefficientSize`, which at `lMax = 256` and
`k = 8` is **8.4 MB per thread**, so eight threads ask for 68 MB of a 16 MB L3;
the inverse has no accumulator and spends a wider chunk on nothing but
amortising the Legendre stage. `Chunking::Count` is handed the same
`bytesPerField` in both directions and knows nothing of this.

**This is [C11] arriving early, and on the wrong machine.** §7 says the
forward's thread-private accumulators are the wrong shape above about sixteen
threads and that the laptop cannot show it. It can: the accumulators are
already the dominant term at eight threads once the chunk is wide, and the
directions have visibly different optima. That is a defect in the chunking
policy — one heuristic serving two decompositions — and it is worth its own
task ahead of anything further in F′, because **it is worth more than F′ was**:
a 2.2× on the inverse direction of the stored path, which is the path
everything actually uses.

**T7 — the benchmark harness of §5.** *Done.* Then steps E–H, each gated on the
measurement that justifies it. **The first run of the harness changed what
those steps should be**; see §9.

The two stages are separated without instrumenting the transform. The
coefficient loop is filtered on the call's truncation degree, so a call at the
smallest legal degree does the full FFT work and almost no Legendre work; the
difference from a full-degree call is the Legendre stage, measured on the
production code rather than a replica.

---

## 9. What the first benchmark run changed

Recorded separately because it revises three findings this document argued
from, and because the revision is a decision for the author rather than
something to act on unilaterally.

**Confirmed.** P1's conclusion: the Legendre stage is 44× the FFT stage in wall
time at `lMax = 256`. Every optimisation should target it. Step C's copy is
invisible (measured at T4). Grid construction is cheap — 0.23 s and 648 MB
resident at `lMax = 256, nMax = 2`, matching F9's arithmetic.

**Corrected.** The stage is *not* bandwidth-bound. It runs at 19% of the
single-thread read rate, and is limited by load/store throughput — about four
memory operations per Wigner value, only one of which reaches DRAM. Three
hand-optimised rewrites of the kernel did not move it.

**Consequently changed.** Tier-1 batching delivers 2.3×, not `k`, and reverses
once the coefficient array leaves last-level cache. Threading looks like the
larger lever, 4–5×. Tier-2 GEMM is where the remaining gain is, but for the
load/store reason rather than the intensity reason, which makes step G's layout
change a prerequisite rather than a follow-on.

**Resolved by doing E and H.** The document ordered E → F → G → H; the numbers
said E → H → F → G, and that is what was done. E and H are complete; F and G
remain, with F still carrying the deadline that field-algebra phase 2 is
blocked on it ([C9]). Threading measured 2.7–4.2× against batching's 2.3×, so
the reordering was worth it, and at `lMax = 128` the stage is now
bandwidth-bound, which changes what is left to win: F's value is now mostly
that phase 2 needs the interface, and G's is that it reduces the traffic
itself.

**The original wording, for the record.** On these numbers the ordering that
buys the most, soonest, is closer to **E → H → F → G**: the plan cache and flag policy first because they are
prerequisites and cheap, then threading for its 4–5×, then batching for its
2.3× and because phase 2 needs the interface, then the layout change and GEMM.
Against that, F is the one with a deadline — field-algebra phase 2 is blocked
on it ([C9]) — and threading a single transform interacts with the layered
layer's parallel-over-slices, which is the nesting rule step H has to settle
anyway. Not reordered unilaterally.
