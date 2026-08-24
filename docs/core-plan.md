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

**T10 — step F, the batched transform primitive.** *Done.* Tier 1 only:
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

**T11 — step F′, Wigner values on the fly.** *Done, both commits.* Two commits.
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

---

## 10. Where the efficiency work stands, and what is left worth doing

Written after T11, when the options had multiplied and the path had stopped
being obvious. It supersedes the ordering advice of §9, which was written in a
regime this section shows we have left.

### The decisive fact: there are two regimes, and they bind differently

**Unbatched, `k = 1`, which is what every measurement before T11 used.** At
`lMax = 256` on eight threads the forward transform moves 136 MB of Wigner
values in 4.26 ms — **32 GB/s against a machine roof near 42**, so 76% of it.
At `lMax = 128` step H already reached 39–40 GB/s, which is 95%. In this
regime, better scheduling, better decomposition and better loops are worth at
most about **1.3×** between them. Everything beyond that has to reduce the
traffic.

**Batched, `k = 8` with a chunk that fits, which is what phases 2–5 will always
do.** The inverse at `lMax = 256` runs at **1.93 ms per field**. The table is
read once per chunk, so that is about 17 MB per field — **8.8 GB/s, a fifth of
the roof**. The arithmetic is of order 68 Mflop per field, so roughly
35 Gflop/s across eight cores, **under a tenth of what they can do**.

So in the regime that production will actually run in, **neither bandwidth nor
arithmetic is saturated**. What binds is the loop structure itself: P1's
corrected finding of about four memory operations per Wigner value, only one of
which reaches DRAM.

### What that does to the options

| lever | at `k = 1` | at `k = 8` | cost | confidence |
|---|---|---|---|---|
| direction-aware chunking | — | **2.2× on the inverse** | small | measured, T11 |
| reduced precision [C6] | ~2× | small | medium | high on traffic; needs a tighter oracle first |
| symmetry reduction [P3, P4] | ~2× | small | large | high on traffic; irregular access |
| transform-major layout + GEMM [P5] | ~1.3× | **potentially large** | large | the only lever aimed at what binds |
| generation, F′ rungs B and C | — | — | large | measured to lose (T11) |
| polar truncation | ~1.5–2× | ~1.5–2× | medium | **the estimate is wrong, and the option is dropped**: §11.7 measures ~1.15× at lMax = 256 |

Two entries need saying out loud.

**The traffic-reduction levers are worth much less than this document assumed**,
because its arithmetic was done in the `k = 1` regime. Halving the bytes per
stored value helps a stage at 76% of the roof; it does very little to one at
20% of it. Batching has already bought most of what precision and symmetry
would buy again, and they do not compose — they are two ways of spending the
same headroom.

**Polar truncation was never in the option set and should be.** `d^l_{nm}(θ)`
falls off like `sin^{|m|}(θ)`, so at high `|m|` the values near the poles are
negligible to any chosen tolerance and those `(m, θ)` pairs can be skipped
outright. It is the only option here that makes the problem *smaller* rather
than making the machine work better: it reduces storage and arithmetic
together, in both regimes, and it is independent of the layout question. It is
a standard lever in production spherical-harmonic libraries and this plan
missed it.

### The order to work in

1. **Direction-aware chunking.** *Done.* `Chunking::Count` divided the cache by
   the thread count, which is the *forward* transform's rule: that direction
   gives every thread a private accumulator of `chunk × coefficientSize`. The
   inverse has no accumulator — its one gathered block is shared and read-only —
   so the same rule starved it, and T11 measured the cost at **2.2×** on the
   batched inverse.

   The parameter is now the number of **copies** of the block that will be live
   at once: the threads for the forward direction, one for the inverse. Both of
   P8's anchors survive unchanged, since both were forward-shaped.

   *Measured, same session, same machine, eight bound threads.* At
   `lMax = 256, k = 8` the inverse goes from **4.23 ms to 2.12 ms per field**,
   a **2.0×**, landing within 10% of what the whole batch as one chunk
   achieves; at `lMax = 128` it goes from 0.31 to 0.25 ms. The forward is
   unchanged at 4.29 ms, which is the point — the same measurement says the
   whole batch there is 2× *slower*, because eight private accumulators of
   8.4 MB ask for 68 MB of a 16 MB cache.

   One thing the rule still cannot see, and the target-machine run should: on a
   multi-CCD or multi-socket machine a *shared* block is pulled into each cache
   domain that touches it, so the inverse's single copy is really one per
   domain — eight CCDs on the deployment target.
2. **The target-machine run.** It settles [C11], and it is the only thing that
   can say whether F′'s memory-and-NUMA case survives. One correction to P8
   while waiting: an EPYC 9334 has *more* memory bandwidth per core than the
   development laptop, not less, so "the stored path is at the DRAM roof while
   a generated one scales with cores" may not hold there either. What survives
   for F′ is NUMA on a 648 MB table, and the 5.4 GB one at `lMax = 512`.
3. **Then the transform-major restructure** — step G's layout together with
   step F's tier 2. It is the only lever pointed at what actually binds in the
   batched regime, and **it subsumes [C11]**: a per-`m` decomposition has
   disjoint outputs, so it needs no accumulator and no reduction, which is
   exactly the forward transform's weak point. This document has carried the
   layout, the GEMM and the decomposition as three separate items; they are one
   piece of work.

   *Written up in full in `gshtrans-reference.tex` §12* — the matrix form, the
   arithmetic-per-memory-operation argument, the layout change, the free
   transpose at the FFT, the real-GEMM-on-complex-data trick, what it does to
   threading and NUMA, what it costs in test strength, and an honest estimate
   (two or three, not an order of magnitude, because the shape is skinny).
4. **Evaluate polar truncation** before [C4]'s symmetry-versus-precision
   question, since it dominates both and is independent of them.

   *Done, and it is dropped.* §11.7 measures it: about **1.15×** at
   `lMax = 256` for the version a GEMM can take, against the 1.5–2× this
   document assumed, and it does not reach that figure at any degree measured.
   Two further facts settle it — the time saved would be less than the
   arithmetic saved, since truncation shrinks the GEMM's inner dimension and
   M5 measured efficiency falling with it; and the tolerance would enter
   [C12]'s cross-kernel oracle, weakening the check that made §11 safe to
   build exactly as the saving grew. The item is closed rather than deferred:
   it was evaluated, which is what this line asked for, and the evaluation
   went against it.

F′'s rungs B and C are not to be built, and [C4] drops below all of the above.

**With item 4 answered, this document has no scheduled work left.** §11 is
built through M6; the target-machine run of item 2 is wanted and gates
nothing. What remains is `thoughts.md`'s list, which is a different document's
business.

### The strategic caveat, which matters more than the ordering

The core is at 76% of the memory roof unbatched and has one clear structural
lever batched. **Phases 2 to 5 of the field algebra have not been started**,
and phase 2 has been unblocked since T10. A further 2× on the transform is
worth less to what this library is *for* than having tensor fields at all, and
the efficiency work has reached the point where each further step costs more
and returns less than the one before it.

So: take item 1 because it is nearly free, take the server numbers because
they are being waited on anyway, and then build phase 2 — leaving the
transform-major restructure as a well-specified piece of work to pick up when
the field layer's real problems ask for it.

### What happened next, 2026-08-23

*The advice above was taken, and the paragraph it rests on is now out of date
in the way one hopes for.* Phases 2 to 5 are built, and so are the
contravariant derivative, the layered fields, the tangential bundle and the
radial seam — `field-algebra-plan.md` §§12–21. The field layer is no longer
the thing that is missing.

Three consequences for this document, none of which change its conclusions.

- **Step F has consumers now**, several of them. The batched transform is used
  by the tensor layer's per-upper-index groups and by the layered layer's
  radial stacks, and §17.5 measured what it is worth there — including the
  case where it collapses to nothing, which is the forward direction at
  `lMax = 256` on eight threads, because the thread-private accumulator no
  longer fits. That is the sharpest argument yet for item 3, and it is a
  measurement rather than the abstract factor of two this section had.
- **Item 2, the target-machine run, is still outstanding**, and is still the
  cheapest thing on the list. The one attempt returned a machine description
  and a STREAM figure and no transform tables at all.
- **The strategic caveat still holds, pointed the other way.** It said the
  field layer's problems were worth more than another 2× on the transform.
  They were, and they have been done; what is left in the field layer is
  small and mostly editorial. So the balance has genuinely shifted towards
  item 3 — which should be started from the server numbers rather than from
  the laptop ones.

### And what happened after that, the same day

*Item 2 is no longer a gate, and item 3 is planned.* The clause above —
"started from the server numbers rather than the laptop ones" — is
**superseded**: `earth-tunya` is waiting on an IT update with no date on it,
and blocking the largest remaining piece of work on an unschedulable
dependency costs more than starting without it.

§11 is the plan. Its central decision is that **both kernels are kept**, so
the loop path and the matrix path coexist as a construction-time policy
rather than one replacing the other — which turns the single server
comparison this document has been waiting for into a comparison available on
any machine, and turns §12's loss of the bit-exact tests into a cross-kernel
oracle. Item 4, polar truncation, stays where it is and stays separate
([C16]).

---

## 11. The transform-major restructure, in detail

Written 2026-08-23, when §10's item 3 was picked up. `gshtrans-reference.tex`
§12 is the mathematics and the arithmetic argument — the matrix form, the
free transpose at the FFT, the real-GEMM-on-complex-data trick, and the
honest estimate — and is not repeated here. This is the work order and the
decisions.

**Two things changed before anything was written**, and both come from the
author rather than from the analysis.

### 11.1 The gate is dropped, and the shape changes with it

§10 item 2 made the target-machine run a prerequisite: start the restructure
from the server numbers rather than the laptop ones. **That gate is removed.**
`earth-tunya` is waiting on an IT update with no date on it, the numbers have
been outstanding across several sessions, and blocking the largest remaining
piece of work on an unschedulable dependency is worse than starting without
it. The measurement is still wanted and §11.6 still asks for it; it is no
longer a precondition.

What replaces it is better than what it replaced. **Both kernels are retained
permanently** ([C12]), so the comparison this document has wanted from a
single server run becomes a comparison anyone can make on any machine at any
time, including this one, today. A gate on one measurement becomes a
mechanism for many.

### 11.2 Decisions taken

**[C12] Both paths are retained, selected at construction by a policy value.**
Not a migration and not a flag day: the loop kernel of steps F and H stays
exactly as it is, and the matrix kernel is a second path beside it. The
precedent is exact — [C10] made stored-versus-generated Wigner values a
construction-time policy for the same reason, that a per-call choice would
forfeit the storage saving and a template parameter would infect every
downstream type with a decision about one object's storage.

Three consequences, and the second is the one worth the most.

- **The layout stops being a change and becomes a choice.** The table is
  built once, in whichever layout the grid was asked for, and the two are
  mutually exclusive per grid — 648 MB at `lMax = 256, nMax = 2` is not a
  thing to hold twice.
- **The loop path becomes the matrix path's oracle**, which inverts §12's
  "it costs the exact-equality tests". That paragraph is right that a GEMM
  sums in an order its kernel chooses, so the batched-against-unbatched tests
  cannot demand bit-exactness *on the matrix path*. But with both kernels
  present there is a comparison that does not exist today: the same inputs
  through two independent implementations of the whole Legendre stage, to a
  tolerance. That covers the layout, the indexing, the FFT ordering, the
  accumulation and the threading decomposition — which is very nearly
  everything this work can get wrong.

  *Stated precisely, because overclaiming it would be easy.* The two paths are
  not independent in the `d` values: both read the same recursion, so a wrong
  value is wrong in both. That half is already pinned elsewhere, by
  `CheckWignerConvention.h` against the `l = 1` table and `CheckLegendre.h`
  against `std::sph_legendre`. What the cross-kernel check adds is everything
  built *on* those values, and that is the part being written here.
- **The cost is two Legendre kernels to maintain**, and it should be named
  rather than absorbed. The inner loop is the thing being duplicated, in both
  directions, with the real and complex cases in each. That is the price of
  the option, it is paid every time either is touched, and it is accepted
  because §11.1's argument needs both to exist at once.

**[C13] The policy is `TransformKernel`, and it names the kernel rather than
the layout.** `TransformKernel::Loop()` and `TransformKernel::Matrix()`,
appended to the grid's constructor and to `ForBand` after `WignerValues`,
defaulted to `Loop()` so that nothing existing moves.

It governs the layout too, since neither kernel works with the other's, and
naming it for the kernel is what makes that legible: the caller chooses an
algorithm and inherits a storage order, not the reverse. It also makes the
one incompatibility statable in one sentence instead of reading as two
storage policies disagreeing —

```cpp
auto grid = Grid(lMax, nMax, flag, chunking,
                 WignerValues::Stored(), TransformKernel::Matrix());

// Refused at construction, and the message says why:
auto bad = Grid(lMax, nMax, flag, chunking,
                WignerValues::Generated(), TransformKernel::Matrix());
```

**`Matrix` and `Generated` are incompatible**, which is a fact about the
recursion rather than an unimplemented case. [C17] states why and says what
happens when a caller asks for both.

**[C14] The GEMM comes from an optional BLAS, found and not fetched.**
`GSHTRANS_WITH_BLAS`, in the pattern `field-algebra-plan.md` §19.5 [R1]
established and §21 exercised: a build without it compiles, links nothing
extra, and simply does not offer `TransformKernel::Matrix()` — the way a build
without `Interpolation` does not offer `SplineDerivative`.

*The first draft of this paragraph ended "the refusal is at construction and
names the option", which contradicts the sentence before it.* If the factory
does not exist there is nothing left to refuse at construction. The two are
different mechanisms for different failures and [C17] separates them.

*This corrects `gshtrans-reference.tex` §12, which is stale.* It says "Eigen
is already a dependency and has a competent GEMM; starting with Eigen and
measuring against a linked BLAS is the cheap order." Eigen left the tree when
GaussQuad dropped it, and nothing in the library names it now. So there is no
free GEMM to start from and the choice had to be made rather than defaulted
into.

Found rather than fetched, unlike the four siblings, because BLAS is not
header-only and not ours: a system OpenBLAS, MKL or Accelerate is what a
caller on a real machine already has, and fetching and building one would be
this library taking on a responsibility it has declined for every other
dependency.

**It is the first non-header-only dependency**, which is a real deviation from
`thoughts.md` §5's stated preference, and it is why the option exists rather
than the dependency. Everything the library does today stays available without
it.

*And the honest caveat, which §12 already states and this decision does not
escape.* Skinny-`N` is the shape general BLAS kernels handle worst, `N = 2k`
is between 2 and 16, and the 513 products per upper index are individually
small enough that call overhead is not negligible. A batched BLAS interface is
the right shape and is not universally available. So a linked BLAS is expected
to be *better* than a hand-written kernel on this shape, not good at it.

**[C17] Two refusals, by two mechanisms, because they are two failures.**
[C13] and [C14] each produce a way for a caller to ask for something they
cannot have, and the first draft of this section treated them as one. They are
not alike: one is a build that cannot do what was asked, the other is two
supported choices that do not compose.

***No BLAS in the build: the factory does not exist.*** `TransformKernel`
carries no `Matrix()` in that configuration, so the failure is a compile error
at the call site rather than a throw at grid construction.

The deciding argument is that the choice is reversible in only one direction.
Adding the factory back later with a runtime throw is strictly widening and
breaks no existing caller; going the other way — removing a factory people
have written against — is not. So the tighter option is the one to start from,
and if the looser one is ever wanted it is available without cost.

Two supporting reasons. The precedent is exact and is already visible in the
test count: a build without `Interpolation` runs 270 tests against 282 because
the spline tests are *not built*, rather than built and skipped. And the
library reaches for this instinct repeatedly — `field-algebra-plan.md` [D8]
chose `Component<0, 1>()` on a tangential tensor being **absent** over being
zero, for exactly the same reason.

*The cost is real and the usual escape hatch does not work here.* A caller
taking their kernel from a configuration file needs their own `#ifdef`. The
obvious remedy — a `constexpr bool` to branch on with `if constexpr` — **does
not rescue it**, because in a non-template function the discarded branch is
still parsed, so a `TransformKernel::Matrix()` sitting in the dead arm still
fails to compile. That caller does not exist yet, and the library code that
does need the branch — the tests, the benchmark, and any later wisdom
mechanism — is one place each. Recorded so that the day it bites, the fix is
known to be a widening one.

***`Matrix` with `Generated`: the grid throws.*** Both are supported by the
build and neither is wrong on its own, so there is nothing to remove; the
combination has to be refused where the two values meet.

**The reason is stronger than §12's, which undersells it.** That section says
generation produces values "in exactly the wrong order" and a supplier feeding
a GEMM would have to transpose. The sharper statement: **the recursion's
output for one `(n, θ)` spans every order at once.** So isolating the single
`(n, m)` block a per-order product wants means either keeping all of it — which
*is* the table, and not having one is the whole purpose of `Generated` — or
re-running the recursion once per order, `2·lMax + 1` times over. It is not an
ordering inconvenience with a transpose as its price; there is no way to
generate what the GEMM needs without paying one of those two.

**Neither silent substitution is available, and it is worth saying which
breaks what.**

- *Honouring `Matrix` and ignoring `Generated`* — building the table anyway —
  can exhaust memory on a caller who chose `Generated` precisely because they
  cannot afford 648 MB at `lMax = 256`, or 5.4 GB at 512.
- *Honouring `Generated` and ignoring `Matrix`* — quietly using the loop —
  lies to the benchmark. That is fatal **specifically under [C12]**: the entire
  justification for carrying two kernels is being able to compare them, and a
  policy that reports the other kernel's numbers under this one's name destroys
  the mechanism this section is built on. A warning does not repair it, because
  warnings go to stderr and no benchmark harness reads stderr.

The second argument also disposes of "fall back to `Loop` with a warning" as a
general answer, and it is worth noting that it is not fastidiousness about
error handling. It is a consequence of [C12] specifically: this library
tolerates a silent substitution less than most, because measurement is what
the design is for.

**[C15] The reflection symmetry is part of this work but not part of its first
measurement.** D&T (C.118) relates the matrices at `±m` by reversing the
colatitude and a sign alternating with `l`, so splitting the input into even
and odd parts about `θ = π/2` halves both the stored table and the arithmetic.
Inside a GEMM that is a clean halving of `K`, which is why P4's objection —
that the same reduction trades bandwidth for irregular access — does not carry
over from the loop.

It is nevertheless M6 and not M3. Landing it in the same change as the kernel
would confound the one measurement the whole restructure exists to produce:
the GEMM against the loop, on the same table size, same arithmetic, same
everything but the kernel. Halve the arithmetic in the same commit and there
is no longer a number that says whether the restructure was worth doing.

**[C16] Polar truncation is not part of this.** §10 item 4 remains its own
piece of work. It is independent of the layout, it applies to both kernels,
and folding it in would make the A/B unreadable for the same reason [C15]
defers the reflection. *§11.7 measures what it would be worth, now that there
is a kernel it would fit, and the answer is much less than §10 supposed.*

### 11.3 What has to change, and one cost §12 understates

The reference note's §12 lists three changes: the Wigner layout, the FFT
output landing in `m`-major order, and the real-GEMM trick. All three stand.
One of them carries a cost that section does not price.

**The intermediate is larger than §12's figure suggests.** It says the
intermediate is `nφ × nθ × k` complex, "which is just the field, 2.1 MB at
`lMax = 256` and `k = 1`". True, and `k = 1` is the case that does not matter:
the matrix kernel exists for the batched regime, and at `k = 8` the same
intermediate is 17 MB — the whole of this laptop's L3, and it is live at the
same time as the table stream and the output block.

That is a genuine tension with the present design rather than a detail. The
loop kernel FFTs one colatitude at a time *precisely* to keep the intermediate
small; the matrix kernel must run all the FFTs first, because a per-`m`
product needs a contiguous `(θ, κ)` block. So the restructure trades a small
working set for a large one, and the chunk rule that `Chunking` applies may
well want a different constant on this path — possibly a smaller optimum `k`,
which cuts against the GEMM's preference for a wider `N`.

**This is a measurement, not an objection**, and it is the first thing M5
should look at. It is recorded here because it is the most likely way for the
restructure to underperform its estimate, and because finding it in the
numbers without having predicted it would waste a session.

### 11.4 The steps

Ordered so that each is separately testable and the one measurement that
matters is not confounded.

**M1 — the layout.** `Wigner` gains `[n][m][l][θ]` as a construction option,
alongside the present `[n][θ][(l,m)]`. Total size is unchanged — summing
`n_L(m)·nθ` over `m` gives the same `66,045 × 257` doubles per upper index at
`lMax = 256`. The existing `Storage` tag orders the `(n, θ)` axes only and is
a different and coarser thing; this does not replace it.

*The test is a walk.* Build both layouts on the same grid and require every
`(n, m, l, θ)` value to be present in each and bit-identical. That is a
complete check of the layout in isolation, it needs no transform, and it is
what makes a later disagreement between the kernels attributable to the kernel.

*Done*, as `GSHTrans/src/WignerMatrices.h`. A class of its own rather than a
layout flag on `Wigner`, because the two have unrelated interfaces — `Wigner`
hands out a `(l, m)` block for one `(n, θ)` through `ConstGSHView`, this hands
out an `(l, θ)` matrix for one `(n, m)` as a span a BLAS call can take — and a
single class serving both would have to branch in the accessor the inner loop
calls. The grid already holds its table in an `optional`, empty on a
generating grid, so a second one beside it is the shape that was already
there.

**The values are generated, never transposed.** Building a `Wigner` table and
transposing it would need both live at once — 1.3 GB at `lMax = 256, nMax = 2`
to end with 648 MB. So the recursion runs into per-thread scratch one
`(n, θ)` block at a time, through `WignerDetails::ComputeBlock` with the same
arguments the stored path passes, and is scattered into place. Same recursion,
same seeds, same evaluation order, same orthonormalisation — so bit-identity
is a property of the construction rather than a hope about it.

*The five tests are not vacuous, and that was checked rather than assumed.*
Perturbing the scatter from row-major to column-major — the single most likely
way to get this wrong — fails four of the five. The fifth is the one asserting
that a matrix's height falls linearly in `|m|`, which is a shape property and
correctly does not notice. Restored, all five pass. The suite goes from 296 to
301, green in Debug and Release.

**The measurement M1 owes, and it is the one cost this step has.** A block
computed for one `(n, θ)` scatters across every matrix at stride `nθ` in the
degree, so nearly every value written touches its own cache line. Both layouts
built on the same angles, eight threads, this laptop:

| lMax | nMax | block layout | matrix layout | ratio | size |
|-----:|-----:|-------------:|--------------:|------:|-----:|
|   64 |    2 |       3.9 ms |        1.8 ms | 0.45× |   10.5 MiB |
|  128 |    2 |      20.4 ms |       24.0 ms | 1.17× |   81.9 MiB |
|  256 |    2 |     156.0 ms |      297.2 ms | 1.90× |  647.5 MiB |

So the scatter costs up to **1.9×** on construction at operator sizes, and the
penalty grows with `lMax` as the write set leaves cache — which is the shape
the argument predicts. It is paid once per grid, it is 140 ms in absolute
terms, and **it is not being fixed**: blocking the scatter over colatitudes so
that a cache line is filled by consecutive `θ` is the known lever, and taking
it now would mean tuning an access pattern before M2 and M3 have had their say
about what it should be. Recorded so that if grid construction ever becomes
the complaint, the answer is already written down.

*The size column also settles a claim this section makes twice*: 647.5 MiB at
`lMax = 256, nMax = 2` in the matrix layout, against the 648 MB the document
quotes for the block one. The transposed triangle really is the same triangle.

*One piece of documentation drift found on the way.* The suite stood at **296**
before this step, not the 282 that `field-algebra-plan.md` §21.3 records. That
number has been overtaken — the three examples of the last commit are part of
it and not all of it — and is left for whoever next touches that section to
correct at its source.

**M2 — the FFT restructure.** All FFTs before the Legendre stage, landing
`m`-major. One `plan_many` over `nθ k` rows with output stride `nθ k` and
output distance 1 writes `[m][θ][κ]` directly, so the transpose is free for
the same reason it was free at tier 1 (P6).

*Separately testable, which is why it is its own step:* it produces the same
numbers as the interleaved path in a different order, so it can be checked
against the current FFT stage before any matrix exists. §11.3's working-set
question is measurable here too, ahead of the kernel.

*Done*, as `ForwardFourierStage` on the grid, beside the transforms rather
than inside them: [C12] keeps the loop kernel exactly as it is, so this is a
second entry point that nothing calls yet and M3 will.

**The oracle is a direct sum, not the loop kernel's FFT stage.** This section
asked for the latter. The former is strictly stronger and cost a dozen lines:
two paths through the same FFTW plan agree even when the plan is the wrong
transform, and say nothing about the sign of the exponent, the normalisation,
or where negative orders sit. A naive DFT written out in the test says all
three, and it passed first time, so those three conventions are now pinned
rather than inherited.

*Seven tests, and four of them bite.* Perturbing the pack from θ-major to
k-major fails three. The two that survive — the interleaved-batch and
sub-range tests — compare two calls perturbed identically, which is worth
recording as a limit rather than a gap: a relative test cannot catch an error
common to both of its sides, and it is the direct sum that covers them.

#### What the measurement changed, which is most of what this step produced

**§12's "the transpose itself is free" is true only away from one hazard, and
badly false at it.** The output write has stride `howMany` complex doubles, so
when `howMany · 16` is a power of two the writes for successive orders collide
in the same cache sets. Measured directly at `nPhi = 520`, per transform:

| howMany | 62 | **64** | 66 | 254 | **256** | 258 | 2046 | **2048** | 2050 |
|---|---|---|---|---|---|---|---|---|---|
| cost | 0.80 | **2.96** | 0.83 | 0.90 | **7.32** | 0.88 | 1.76 | **8.53** | 1.65 |

Four to eight times, at every power of two and nowhere else. **This is the
hypothesis `field-algebra-plan.md` §17.7 raised for `RadialMajor` and
rejected**, having found `nR = 255, 256, 257` identical — correctly, because
its tiling already handled it. Nothing tiles here: FFTW writes straight
through. So the same hypothesis is right in one place and wrong in the other,
and the difference is whether anything sits between the write and the cache.

**A small block beats one big call, which is the opposite of §12's design.**
That section specifies a single `plan_many` over `nθ k` rows. Measured, one
thread, `k = 8`:

| lMax | stage, block 3 | stage, all θ | FFT workspace, block 3 | all θ | full forward |
|-----:|---------------:|-------------:|-----------------------:|------:|-------------:|
|   64 |        0.11 ms |      0.23 ms |                 0.1 MiB | 2.1 MiB |     0.67 ms |
|  128 |        0.59 ms |      1.15 ms |                 0.2 MiB | 8.2 MiB |     6.00 ms |
|  256 |        3.25 ms |      5.46 ms |                 0.4 MiB | 32.6 MiB |    52.01 ms |

So the small block is **1.7× to 2× faster and wants eighty times less
workspace**. A strided write over a few hundred bytes stays inside a couple of
cache lines; one over tens of kilobytes does not. The default is therefore a
block of three, guarded to shrink further if the product would alias, and
`thetaBlock` is a hint in the sense `Chunking` is rather than an instruction.

**§11.3 priced the working set at half of what it is.** It put the
`m`-major intermediate at 17 MB at `lMax = 256, k = 8` and treated that as the
cost. There are *two* buffers of that size in a full-height call — the packed
input as well as the transposed output — so it was 33 MB, per thread. Blocking
removes that objection rather than answering it: the intermediate `out` stays
16.3 MiB because the problem fixes it, and everything else drops to 0.4 MiB.

**And the headline for M3.** The Fourier stage is **3.25 ms against a 52 ms
forward transform** at `lMax = 256, k = 8` on one thread — about six per cent.
So the restructure's cost is not in this half, and whatever M3 measures will
be the Legendre stage's doing. That is the number M5 should hold on to.

*The copy blocking costs is not the problem either*, which was the other thing
worth ruling out: split three ways, the copy is 0.7–1.9 ms against the FFT's
1.65–18 ms at `lMax = 256`, and it falls as the block grows. Both sides are
contiguous in `(θ, k)`, so it is a `memcpy` per order and not a gather.

**M3 — the matrix kernel.** Per-order `dgemm` with `N = 2k`, one stored matrix
serving both directions through a transpose flag. `TransformKernel::Matrix()`,
`GSHTRANS_WITH_BLAS`, and the construction-time refusals of [C13] and [C14].

*The test is [C12]'s oracle:* the same inputs through both kernels, forward
and inverse, real and complex, batched and unbatched, agreeing to a tolerance.
Plus the existing round-trip and analytic tests run on a matrix grid, which
costs nothing to add and covers what the cross-kernel check cannot.

**Split in two, as T2 was.** M3a is the forward direction and all the
infrastructure — the build option, the policy, the grid wiring, the BLAS
binding — and M3b is the inverse, which reuses every part of it.

#### M3a, done

`GSHTRANS_WITH_BLAS`, `Blas.h`, `TransformKernel` in `Policies.h`, and
`ForwardMatrixKernel` on the grid. The forward transform agrees with the loop
kernel to `1e-13` relative across complex and real fields, batched and
unbatched, on a full grid and truncated below it, and on the reduced `m ≥ 0`
scalar grid. Perturbing the negative-order mapping fails the three complex
cases and correctly not the two real ones, which have no negative orders.

**The option is three-state, not two.** `AUTO` — the default — uses a BLAS if
one is there and says so if not; `ON` demands one and fails the configure;
`OFF` never looks. `ON` has to mean something different from the default
because a caller who typed it deserves a failed configure rather than a silent
downgrade to the loop kernel, which is [C17]'s argument applied to the build.

**A third refusal, which [C17] did not anticipate: BLAS has no `long
double`.** The interface is `s`, `d`, `c`, `z` and nothing wider, so a grid at
that precision cannot have the matrix kernel however it is configured. It is
refused at construction like the other two, and the dispatch carries an
`if constexpr` so that the body is never instantiated there. This is a real
narrowing — the loop kernel serves `long double` and always will — and it
belongs beside [C14]'s other costs rather than in a footnote.

**Off is a supported configuration and was checked, not claimed** — and the
check earned its place immediately: the first version of the precision
refusal was outside the `#ifdef` and broke the no-BLAS build, because
`BlasDetails` does not exist there. With the option off the suite runs 308
tests against 314, the six being the cross-kernel ones.

**The weights are applied to the intermediate, not folded into the matrix.**
They cannot be folded: the same matrix serves the inverse, which carries no
weight, and one stored matrix serving both directions is the whole reason the
layout is worth having. So it is one extra pass over the intermediate, which
is 16 MiB at `lMax = 256, k = 8`.

**The scatter cannot be avoided either.** The coefficient block is triangular
— the stride between consecutive degrees at one order is `2l + 1` and so is
not constant — and a GEMM writes with one leading dimension. So the product
goes to a small per-order buffer and is scattered from there, which is the
same volume of writing the loop kernel's own scatter does.

#### What it is worth, measured

Forward, complex, `n = 2`, one thread, BLAS pinned to one thread:

| lMax | k | loop | matrix | speedup |
|-----:|--:|-----:|-------:|--------:|
|   64 | 1 | 0.15 ms | 0.06 ms | 2.5× |
|   64 | 8 | 0.67 ms | 0.25 ms | 2.7× |
|  128 | 1 | 1.90 ms | 0.58 ms | 3.3× |
|  128 | 8 | 5.88 ms | 2.07 ms | 2.8× |
|  256 | 1 | 10.27 ms | 5.17 ms | 2.0× |
|  256 | 8 | 49.15 ms | 17.86 ms | 2.8× |

**Two to three and a bit, which is exactly what §12 predicted** — and worth
saying, because this document's measurements have more often contradicted its
predictions than confirmed them.

**But the gain is there unbatched too, which §12 predicted it would not be.**
That section says the restructure "buys nothing at all in the regime where the
stage is bandwidth-bound — which it is, unbatched, on eight threads at
`lMax = 128`". The qualifier is the whole of it: *on eight threads*.
Sequentially the stage is bound by per-core load/store throughput and not by
DRAM — P1's finding — so it sits far below the roof, and the GEMM's better
arithmetic per memory operation helps at `k = 1` as much as at `k = 8`. §12's
prediction is a claim about the threaded regime and **remains untested until
M4**, which is where the matrix kernel learns to thread. Nothing here
contradicts it.

**A threaded BLAS loses, and by more than nesting explains.** §11.4's M4 note
says the BLAS must be single-threaded when GSHTrans threads over orders,
because two thread pools multiply. Measured, it is worse than that: with
GSHTrans **sequential**, so with no nesting at all to blame, giving OpenBLAS
eight threads costs 23.72 ms against 16.97 at `lMax = 256, k = 8`, and 0.15
against 0.10 at `lMax = 64, k = 1`. The 513 products per upper index are
individually skinny — `N = 2k` is at most 16 — so thread launch dominates what
the threads could win. **So the obligation is simpler and stronger than M4
stated it: the BLAS should be single-threaded here always, not merely when
GSHTrans threads.**

*One caveat on all of the above.* Run to run these move by 10–15%, which is
this laptop's noise floor and is smaller than every difference quoted.

#### M3b, done, and with it M3

`InverseFourierStage` and `InverseMatrixKernel`. **One stored matrix serves
both directions**, which was the layout's central claim and is now code: the
forward multiplies by `D` and the inverse by `Dᵀ`, and the difference is the
transpose flag in the BLAS call. Nothing is copied and there is no second
table.

Twelve cross-kernel tests now, six each way, plus a round trip on a matrix
grid — which the cross-kernel tests cannot supply, since two kernels wrong in
the same way would agree with each other.

**The inverse has a hazard the forward does not, and it is the band of orders
no coefficient reaches.** The transform is over `nPhi` orders whatever the
degree, so the orders between `lMax` and `nPhi − lMax` are read by the FFT and
written by nothing, and the intermediate is kept between calls — so they hold
whatever the last call left. Zeroing that band is what the loop kernel is
doing when it clears its whole row buffer. Removing it fails two tests, both
of them the ones at a degree below the grid's, where the band is widest; the
full-degree cases pass without it by the accident of the buffer being zero
already, which is exactly why the truncated cases are in the suite.

#### The inverse gains less than the forward, and the shape says why

Inverse, complex, `n = 2`, one thread:

| lMax | k | loop | matrix | speedup |
|-----:|--:|-----:|-------:|--------:|
|   64 | 1 | 0.19 ms | 0.05 ms | 3.8× |
|   64 | 8 | 0.66 ms | 0.24 ms | 2.8× |
|  128 | 1 | 1.96 ms | 1.34 ms | 1.5× |
|  128 | 8 | 5.36 ms | 2.98 ms | 1.8× |
|  256 | 1 | 13.91 ms | 6.32 ms | 2.2× |
|  256 | 8 | 44.14 ms | 30.06 ms | **1.5×** |

Against the forward's 2.0× to 3.3×, and worst exactly where the forward was
best. **The likely cause is structural rather than incidental, and it is not
in §12 at all.** The two directions have the same matrices and opposite
shapes:

    forward:   (n_L × n_θ) · (n_θ × 2k)     M = n_L,  K = n_θ,  N = 2k
    inverse:   (n_θ × n_L) · (n_L × 2k)     M = n_θ,  K = n_L,  N = 2k

A GEMM amortises its operand loads over the inner dimension, so `K` is what
decides how well it does. Forward, `K = n_θ` is a few hundred at every order.
Inverse, `K = n_L` **falls linearly in `|m|` and reaches one**, so half the
products have almost no inner dimension to amortise over and the call is
nearly all overhead. §12 records the two shapes and does not notice that this
follows from them.

*Stated as the likely cause and not the established one.* The other candidates
were priced and are too small: the zero band is 7 orders of 230 KB at
`lMax = 256`, and the coefficient gather is the mirror of the forward's
scatter and so cannot explain a difference between them. Confirming it wants
the per-order products timed against their own `K`, which is M5's work and not
a reason to hold M3.

*If it is confirmed, the answer is a batched BLAS call* — one invocation for
the whole set of orders rather than 513 — which §12 already names as the right
shape and not universally available. That is M6's neighbourhood, not this
step's.

**Suite 314 to 320 with BLAS, 308 without.** Green in Debug, in Release, and
under ASan and UBSan.

**M4 — threading over orders.** Products at different `m` write disjoint
outputs, so there is no accumulator and no reduction in either direction —
which is how this subsumes [C11], by deleting the thing [C11] was a question
about rather than by answering it.

Two cautions from §12, both to be built in rather than discovered. The work
per order is not constant: `n_L(m)` falls linearly in `|m|`, so a static split
over orders is about twice as unbalanced as it looks and the schedule must
divide work rather than count orders. And a thread's share of the table is a
set of whole `(n, m)` blocks, so first touch during the table build can be
made to match the split exactly — which the current colatitude split cannot,
and which is the better NUMA story this path has to offer.

*A third, which M1 created and which is not in §12 at all.* `WignerMatrices`
builds in parallel over `(n, θ)`, because that is the axis the recursion runs
on, and **scatters** into the `(n, m)` blocks. So first touch is by colatitude
and does not match the split M4 wants, exactly as it does not today. The fix
is cheap and known — a first-touch pass parallel over `(n, m)` before the fill
— and it belongs here rather than in M1, since M4 is where the split is
decided and a first-touch pass that matches no split is worth nothing.

*And M4's answer to it is that the two wants are in conflict, which is the
part worth carrying to the target machine.* The schedule chosen below is
`dynamic`, because no static work model survives contact with what M3b found.
But a dynamic schedule has **no fixed thread-to-order mapping**, so there is
nothing for a first-touch pass to match, and the NUMA story §12 offers is not
available while it is in use. The alternative — a static schedule over orders
dealt round-robin after sorting by height, which is both work-balanced and
fixed — is buildable and is not built, because its whole benefit is on a
machine with more than one memory domain and this laptop has one. **The
schedule question and the NUMA question are one question, and it is a
target-machine question.**

**A fourth obligation, on the BLAS rather than on this code, and it has to be
written on the seam.** M4 threads over orders with OpenMP and calls the GEMM
from inside that region. A BLAS with a thread pool of its own then multiplies
the two: eight OpenMP threads each entering a call that wants eight threads
asks for sixty-four on eight cores. That is precisely what `Execution`'s "both
at once is worse than either" rule exists to prevent, and the library's
`omp_in_parallel()` guard cannot see it, because the second pool is not
OpenMP's.

So: **GSHTrans threads over orders and requires the BLAS to be single-threaded
when it does.** A BLAS built on the same OpenMP runtime satisfies this for
free — `libgomp` defaults to one active level, so the inner region runs serial
precisely because the outer one is open — and one built on its own pthreads
does not, and needs `OPENBLAS_NUM_THREADS=1` or its equivalent from the
caller. This is the same kind of obligation `field-algebra-plan.md` §19.2
wrote onto `RadialOperator`: something no caller can honour without being
told.

*Checked on the development machine rather than assumed, and it was not what
installing the package suggested.* Ubuntu ships OpenBLAS in pthread, OpenMP
and serial builds as alternatives of one library. Installing
`libopenblas-openmp-dev` **does not select it**: pthread carries priority 100
against OpenMP's 95, both alternatives sit in auto mode, and the link resolves
to pthread as before. Both the runtime and the link-time symlink have to be
set explicitly. Worth recording because the failure is silent — the build
succeeds, the answers are right, and only the timings are wrong.

**The road not taken, recorded so the choice reads as one.** The alternative
is to let a threaded BLAS parallelise each product and not thread over orders
at all. It is rejected because the 513 products per upper index are
individually skinny — `N = 2k` is between 2 and 16 — which is the shape BLAS
threading handles worst, and because threading over orders is what gives this
path the NUMA story above. Neither reason would survive a much larger `k`, and
neither is likely to meet one.

#### M4, done

`OverOrders` runs the per-order body threaded or not, and both kernels are
written against it. The orders read a shared, read-only table and a shared,
read-only intermediate and write output no other order writes, so there is
**no accumulator and no reduction in either direction** — [C11]'s question
deleted rather than answered, which is what §10 item 3 meant by subsuming it.

**`schedule(dynamic)`, because the obvious static split is wrong twice over.**
§12 warns once: `n_L(m)` falls linearly in `|m|`, so counting orders is about
twice as unbalanced as it looks. M3b found the second and sharper reason. A
static split weighted by `n_L` assumes time proportional to `n_L` — but the
inverse's inner dimension *is* `n_L`, and a GEMM's efficiency falls with its
inner dimension, so time is superlinear in `n_L` and a linear model mis-splits
in the direction it was correcting. Dynamic holds no model and so cannot hold
a wrong one. Determinism is unaffected, since the orders write disjoint
output, and a test asserts the threaded answer is **bit-identical** to the
sequential one — which is the right standard within one kernel even though it
is not available across two.

**The scratch is per thread and taken inside the region; anything shared is
captured from outside it, and must be.** `OrderScratch` and `MatrixScratch`
are both `thread_local`, so a buffer filled before the region belongs to the
master thread and reaching for it again inside would find an empty one. That
is a trap the shape of this code sets and the comment says so.

**The Fourier stage had to thread too, which M2 did not anticipate.** M2 left
it sequential and measured it at six per cent of the transform — true, and
true only while the Legendre stage was sequential as well. Once the orders
threaded it became a **29 per cent Amdahl term**, and the matrix kernel
stopped scaling past two threads. Threading it over colatitude blocks is five
lines, since the blocks are disjoint on both sides and the workspace and its
plan are already `thread_local`. What that changed, `lMax = 256, k = 8`,
eight threads: forward 3.23× to **4.59×**, inverse 1.40× to **1.92×**.

#### What M4 measured, and it settles §12's open prediction

Forward and inverse, complex, `n = 2`, speedup of matrix over loop:

| lMax | k | 1 thread | 2 | 4 | 8 |
|-----:|--:|---------:|--:|--:|--:|
| 128 | 1 | 2.79× / 2.25× | 2.71× / 2.23× | 4.11× / 2.72× | 3.93× / 2.65× |
| 128 | 8 | 3.01× / 1.83× | 2.45× / 2.20× | 3.50× / 2.66× | **4.27× / 3.07×** |
| 256 | 1 | 3.16× / 1.90× | 2.19× / 1.55× | 1.40× / 1.17× | **1.20× / 1.10×** |
| 256 | 8 | 3.02× / 1.64× | 4.10× / 1.78× | 5.18× / 2.08× | **4.59× / 1.92×** |

**§12's prediction is confirmed, and precisely.** It says the restructure
"buys nothing at all in the regime where the stage is bandwidth-bound — which
it is, unbatched, on eight threads". At `lMax = 256, k = 1` on eight threads
the answer is **1.20× and 1.10×**, which is nothing, and it is the only cell
in the table where that is true. M3a's contrary-looking sequential result was
never a contradiction and this is why.

**The forward exceeds §12's estimate, and the reason is [C11] rather than the
GEMM.** That section's honest expectation was "a factor of two or three in the
batched regime". Sequentially the forward gives 3.02×, inside it. Threaded it
gives **4.59×**, outside it — and the extra is not the matrix kernel getting
better but the loop kernel getting worse: at `lMax = 256, k = 8` the loop
forward scales only 1.48× from one thread to eight, which is exactly
`field-algebra-plan.md` §17.5's collapse, eight private accumulators of 8.4 MB
in a 16 MB cache. The matrix kernel has no accumulator to collapse and scales
2.24×. So the threaded comparison measures the GEMM *and* the deletion of
[C11]'s problem together, and the plan said it would.

**The inverse still lags, consistently with M3b.** 1.92× against the forward's
4.59× at the operator size, and the gap is the same one M3b attributed to the
inner dimension falling to one. Nothing here contradicts that and nothing here
confirms it either; it is still M5's to settle.

**M5 — measure.** Both kernels, both directions, batched and unbatched, over
the `lMax` range the `transforms` and `batching` sections already walk, as a
new named benchmark section. §11.3's working set first.

*What the numbers have to be read against.* §12's estimate is a factor of two
or three in the **batched** regime and **nothing at all unbatched**, where the
stage is already at 39–40 GB/s against a roof near 42. A run at `k = 1`
showing no gain is the prediction coming true, not the restructure failing.

**So the reporting is part of the step, not a presentational afterthought.**
Three requirements, each guarding against a specific way this measurement can
be misread or can mislead.

- **Two tables, batched and unbatched, not one table with `k` as a row.**
  A single table puts a column of `1.0×` at `k = 1` next to the gains at
  `k = 8`, and the null result reads as failure to anyone scanning it — which
  will include whoever reads it a year from now.
- **Carry achieved GB/s against the roof the `stream` section already
  measures.** This is the requirement that does the real work, and it is a
  column the existing tables in this document do not have. A caption asserting
  that no gain is expected unbatched is something a reader must take on trust;
  a row showing 40 GB/s against a 42 GB/s roof **demonstrates** that the null
  result is the ceiling. That is the same move §17.5 of
  `field-algebra-plan.md` made when it explained the batched forward's
  collapse to 1.12× by the chunk arithmetic rather than filing it under
  bandwidth — the number that explains the number.
- **Run each kernel in its own process invocation.** This is a hazard in the
  A/B rather than in its presentation. The policy is construction-time, so
  comparing kernels means two grids, and at `lMax = 256` that is 648 MB of
  table each: both live costs 1.3 GB, while building them in sequence leaves
  the second starting on a cold cache with different first-touch placement.
  §17.7 has already been caught by exactly this, where allocating a buffer
  inside the timed loop made a measurement 3–5× worse and never win — and it
  took the numbers looking wrong to notice. The section-name mechanism already
  supports one kernel per invocation, so this costs nothing. Failing that,
  run both orders and say so if they disagree.

*Done*, as the `kernels`, `kernels-loop` and `kernels-matrix` sections, with
the per-order `Gflop/s` table beside them. All three requirements above are
built in.

#### Three defects in the harness, found by using it

**The default build directory measures nothing, and says so now.**
`cmake -S . -B build` leaves `CMAKE_BUILD_TYPE` **empty**, so the benchmark
built there runs about ten times slow — and every figure is internally
consistent, so nothing looks wrong. The first run of this section reported
speedups of *forty*. The harness now prints a warning when built without
`NDEBUG`. This is not new and is not confined to the new section: the
pre-existing `transforms` section in the same build reports `lMax = 256` at
151 ms against 15 ms in a Release build.

**The read-scan roof was measuring floating-point latency, not memory.**
`ScanBandwidthGBs` accumulated into one variable, which is a dependency chain
the compiler may not reassociate, so at one thread it reported **12.7 GB/s
against a triad of 36.2 on the same machine** — a threefold understatement. It
came right at four threads and above, where several chains overlap. Fixed with
four accumulators; the scan and the triad now agree at one thread, 36.7 against
35.6.

*This matters beyond the harness, and the damage is bounded.* §10's central
claim — 39–40 GB/s "against a machine roof near 42", so 95% of it — was taken
at eight threads, where the old helper was already right. **That claim
stands.** Any single-thread comparison made against that column did not, and
this document made none.

**And one the library owned rather than the harness.** The sequential rows
reported table traffic at twice single-core bandwidth, because with an OpenMP
BLAS linked, "GSHTrans sequential" did not mean "BLAS sequential": with no
outer region open, the BLAS took the whole machine. `OPENBLAS_NUM_THREADS` is
inert in that build, and there is no portable call to set a BLAS's thread
count. `OverOrders` now opens a team **even for the sequential case**, so
every GEMM is issued from inside an OpenMP region and a same-runtime BLAS is
nested and serial without anyone setting anything. Verified: sequential
`lMax = 256, k = 8` gives 19.6 ms at `OMP_NUM_THREADS=16` and 20.8 at 1, where
before the fix the two differed.

#### The measurement

Speedup of matrix over loop, with the matrix kernel's table traffic and the
machine's roof beside it:

*Unbatched, `k = 1`, `lMax = 256`:*

| threads | fwd | inv | GB/s | roof |
|--------:|----:|----:|-----:|-----:|
| 1 | 2.95× | 1.94× | 25.8 | 39.9 |
| 4 | 2.03× | 1.50× | 46.3 | 44.4 |
| 8 | **1.53×** | **1.14×** | **41.0** | **41.0** |

*Batched, `k = 8`, `lMax = 256`:*

| threads | fwd | inv | GB/s | roof |
|--------:|----:|----:|-----:|-----:|
| 1 | 2.95× | 1.69× | 6.8 | 39.6 |
| 4 | 5.35× | 2.18× | 14.6 | 44.0 |
| 8 | **4.10×** | **2.22×** | **13.8** | **42.2** |

**The roof column does exactly what it was put there to do.** Unbatched on
eight threads the matrix kernel runs at **41.0 GB/s against a 41.0 GB/s
roof** — it is *at* the ceiling, so the 1.53× is everything there was to win
and no further arrangement can add to it. Batched it runs at a third of the
roof, so the gain there is not bandwidth at all. Two regimes, separated by a
number rather than by an assertion, which is what §11.4 asked for.

*A refinement of §12, which was right in substance and wrong in detail.* It
predicted no gain unbatched **because the loop kernel is already at the
roof**. Measured, the loop kernel is at 65% of it and the *matrix* kernel is
at 100%, which is why there is still a 1.53× rather than nothing. The
conclusion — that nothing further is available in that regime — is
unaffected and is now demonstrated rather than predicted.

#### M3b's hypothesis is refuted

M3b attributed the inverse's smaller gain to its inner dimension falling to
one. The per-order table settles it, and against that hypothesis. Equal flops
both ways, so `Gflop/s` compares efficiency directly, at `lMax = 256, k = 8`:

| n_L | 255 | 193 | 129 | 65 | 17 | 3 |
|---|---|---|---|---|---|---|
| forward | 94.9 | 74.6 | 74.5 | 74.4 | 71.2 | 54.9 |
| inverse | 92.3 | 74.8 | 76.3 | 74.0 | 71.0 | 56.5 |

**Identical, at every height.** Efficiency does fall with `n_L` — from 95 to
55 — but it falls *the same way in both directions*, because `n_L` is `M` for
the forward and `K` for the inverse and a small dimension costs the same
either way. There is no asymmetry in the products to explain an asymmetry in
the transforms.

**What actually explains it is the loop kernel, not the matrix kernel.** At
`lMax = 256, k = 8` on eight threads the matrix kernel's two directions are
9.86 ms and 9.82 ms — the same to within noise. The **loop** kernel's are
40.4 ms and 21.8 ms: its inverse is nearly twice as fast as its forward,
because the inverse carries no thread-private accumulator and §10's
direction-aware chunking already gave it what it needed. So "the inverse gains
less" was never a statement about the matrix kernel. It is the loop kernel's
inverse having less room to improve, and the right reading is that the matrix
kernel makes the two directions cost the same where they did not before.

*One residual, named and not established.* Sequentially the matrix kernel is
still asymmetric — 19.9 ms forward against 32.7 inverse at `lMax = 256,
k = 8` — and the per-order table says it is not the products. The leading
candidate is that the inverse's products write into the 16 MiB intermediate,
520 separate blocks each touched once, while the forward's write a small
buffer reused at every order and scatter separately. It disappears at four
threads and above, so it is not worth chasing.

**M6 — the reflection**, per [C15], if M5 says the path is worth deepening.

#### M6 assessed, and the assessment changes what it is for

[C15] made this step conditional on M5, so M5's numbers get a say before any
of it is built. Two things were established first.

**The relation holds, exactly, in this library's convention.** D&T (C.118) in
the stored values reads

    d^l_{nm}(π − θ) = (−1)^{l+n} d^l_{n,−m}(θ)

and it is verified rather than assumed: worst absolute difference **3.8e-15**
on values of order one, over every `(n, m, l, θ)` at `lMax = 12, nMax = 2`.

**But it relates different orders, not one order at mirrored colatitudes, and
that changes the arithmetic argument.** For a scalar field it is the familiar
same-`m` symmetry that halves the colatitude sum. At `n ≠ 0` it maps `(n, m)`
to `(n, −m)`, so what it gives is that **the matrix at −m is the matrix at +m
with its columns reversed and an `l`-alternating sign**. Working the forward
sum through,

    f^n_{l,−m} = (−1)^{l+n} Σ_j w_j D^(n,m)_{lj} F_{−m}(θ_{j̄})

— the same matrix, applied to the `−m` Fourier data in reversed colatitude
order, with a sign flip on the output rows.

*So §12's summary of the saving is half right.* It says the reduction halves
"both the stored table and the arithmetic". **The table halves; the arithmetic
does not.** The same `2·lMax + 1` products are still done, of the same shapes.
What is saved is the table, in size and in traffic.

**What is gained instead is better than what was claimed, and M5 is why.** The
pair at `±m` can go in **one** GEMM rather than two, with the two right-hand
sides side by side — which doubles `N` from `2k` to `4k`. M5's per-order table
is what makes that interesting: efficiency there runs at 55–95 Gflop/s, far
below peak, and `N = 2k` between 2 and 16 is the skinniest dimension in the
problem and the one §12 already names as what general BLAS kernels handle
worst. Doubling it attacks exactly the limit M5 measured.

**The regime argument, from M5's roof column.**

- *Unbatched, many threads:* the matrix kernel is **at the roof** (41.0 GB/s
  against 41.0). Halving table traffic is the only lever that can help there,
  and it could give up to another 2×.
- *Batched, many threads:* a third of the roof, so traffic is not the
  constraint and halving it buys little. The `N`-doubling is what would help
  here, and it is unmeasured.
- *Memory:* 648 MB to 324 MB at `lMax = 256, nMax = 2`, and 5.4 GB to 2.7 GB
  at `lMax = 512`. This is regime-independent and is the clearest benefit of
  the three.

**The recommendation, and it is a decision rather than work.** M6 is worth
doing, but for **memory and for `N`**, not for the arithmetic halving [C15]
recorded — and the case is weakest in exactly the regime phases 2–5 run in.
Against it: this is the numerically delicate step, [C15] and step G both say
so, and it touches every part of M1 to M4 at once — a second table layout,
both kernels pairing orders, a colatitude-reversed right-hand side, and a
sign pass on half the output.

*What makes it safe to attempt whenever it is wanted* is that the oracle
already exists. [C12]'s cross-kernel comparison would catch any sign or
index error immediately, which is the argument for having kept both kernels
and is worth noting as it paying off a second time.

#### M6, built

`WignerMatrices::Reflected` stores non-negative orders only and checks the
angles support it; `Sign(l, n)` is the reflection, in one place, so no caller
writes it out. Both kernels take the `−m` product from the `+m` matrix, with
the colatitudes reversed and the sign applied — on the way out in the forward
direction, on the way in in the inverse.

*The oracle paid for itself.* All fourteen cross-kernel tests passed **first
try**, and flipping the sign or dropping the reversal fails six of them each.
That is the whole argument for [C12] arriving: this is the numerically
delicate step, and it was checked against an independent arrangement of the
same sum rather than against a tolerance on a property it asserts about
itself.

**The table halves, as predicted**: 647.5 → 325.0 MiB at `lMax = 256,
nMax = 2`, and 5150 → 2580 MiB at `lMax = 512`. A shade over half because
order zero is its own reflection and is stored once either way, which is what
the count test asserts.

#### The first version was slower sequentially, and the measurement said why

Pairing `±m` into one GEMM — the `N`-doubling this section argued for — needs
both right-hand sides adjacent, so **both** have to be copied into a buffer.
Interleaved A/B on an idle machine, matrix-kernel forward at `lMax = 256,
k = 8`: **1.5× better on eight threads and 1.5× worse on one**.

That asymmetry is the answer. If the `N`-doubling were buying anything it
would show sequentially, where arithmetic efficiency is what is left; it
showed the opposite. So the threaded gain is the **halved traffic**, and the
paired copy is pure cost — consistent with M5, which found the batched regime
at a third of the roof and therefore not traffic-bound, and the threaded one
closer to it.

**So the pairing was taken out and only the halving kept.** Each order does
two GEMMs against one matrix: `+m` multiplied where it already lies, `−m`
through a reversed copy, which is the one copy the reflection genuinely
requires since no BLAS takes a negative stride. `N` stays at `2k`.

#### What M6 is worth, interleaved, three passes

Matrix-kernel time at `lMax = 256, k = 8`, before against after:

| | before | after | gain |
|---|---|---|---|
| forward, 1 thread | 19.7, 19.5, 19.2 ms | 16.8, 16.9, 16.9 | 1.16× |
| forward, 8 threads | 8.63, 8.34, 8.65 ms | 6.15, 6.18, 6.10 | **1.39×** |
| inverse, 1 thread | 32.7, 33.6, 33.2 ms | 21.6, 21.4, 21.3 | **1.55×** |
| inverse, 8 threads | 9.37, 8.87, 8.89 ms | 6.58, 6.86, 6.51 | **1.36×** |

*The sequential inverse gains most, and not from the reflection.* Splitting
the pair let the `+m` product write **straight into its block of the
intermediate** instead of through a result buffer — which is how the kernel
worked before M6 for the `+m` half and is a copy M3b had introduced without
noticing. The reflection is what made the inefficiency visible.

Against the loop kernel, `lMax = 256, k = 8` on eight threads, the restructure
now stands at **5.8–6.0× forward and 2.9–3.2× inverse**, from 4.0–4.2× and
2.3× before this step.

**[C15]'s claim was wrong in a way worth keeping.** It said the reduction
halves the table *and the arithmetic*. The table halves. The arithmetic does
not, and the attempt to convert the symmetry into arithmetic — the paired
GEMM — measured worse. What the reflection buys is traffic and memory, which
is what §12's own DRAM argument should have predicted.

---

### 11.7 Polar truncation, priced before it is planned

§10 put polar truncation in the option table at **~1.5–2×, in both regimes**,
and called it "the only option here that makes the problem *smaller*". The
matrix kernel is where it would fit: the intermediate is `[m][θ][k]`, so a
band of colatitudes at fixed order is **contiguous**, and truncating it is a
pointer offset and a smaller `K` — nothing else changes. It is awkward in the
loop kernel and natural in this one.

So it is worth pricing before it is planned. Measured on the built table,
`nMax = 2`, as the fraction of values retained when the colatitudes are
trimmed to where the values exceed a tolerance times the largest in their
matrix:

| lMax | rectangular, 1e-15 | 1e-8 | per-degree, 1e-15 | 1e-8 |
|-----:|---:|---:|---:|---:|
|  64 | 94.7% | 91.0% | 92.6% | 86.8% |
| 128 | 90.6% | 87.0% | 86.2% | 80.5% |
| 256 | 86.7% | 83.9% | 80.1% | 75.4% |

and the per-degree ceiling out to larger degrees, at `nMax = 0`: **75.1%** at
`lMax = 512` and **71.3%** at 1024.

**Three things follow, and together they say this is not the lever §10 thought.**

- **The buildable version saves about 13% at `lMax = 256`, not a factor of
  1.5.** A GEMM takes a rectangle, so the band has to be one band per order —
  and a rectangle is set by its *widest* row. The per-degree trim is
  triangular, needs a product per row-block rather than one per order, and
  still only reaches 1.25× there.
- **It does not improve much with degree.** The ceiling is 1.40× at
  `lMax = 1024` and flattening. §10's 1.5–2× is not reached at any size
  measured. The reasoning was right — `d^l_{nm}` does fall off like
  `sin^{|m|}θ` — but the decay to a tolerance worth having is slower than the
  argument suggests: at `1e-15`, `sin^{256}θ` is still above threshold over a
  third of the range.
- **And the time saved would be less than the arithmetic saved**, which is the
  decisive point. Truncation shrinks `K`, and M5 measured GEMM efficiency
  falling with the inner dimension — 95 Gflop/s at `n_L = 255` against 55 at
  3. So a 13% cut in flops buys less than 13% in time, and possibly nothing.

*What it would cost to build*, for the record, since the shape is now clear: a
band per `(n, m)` in `WignerMatrices` and the rule that sets it, an offset and
a smaller `K` in the forward kernel, the same plus zeroing the out-of-band
rows in the inverse, a tolerance policy, and tests. Comparable to M1 and M3a
together — the smaller half of what §11 took, since the infrastructure exists.

*And one cost that is not effort.* The tolerance would enter [C12]'s oracle:
the two kernels would agree only to the truncation tolerance rather than to
`1e-13`, so the check that made §11 safe to build gets weaker exactly as the
saving grows. That argues for a default tolerance tight enough to leave the
oracle intact, which is also the setting that saves least.

**So it is dropped, and the rationale is here rather than in a decision
line.** §10's entry stands corrected. The option is real and the mechanism is
sound; it is worth about 1.15× where this library is used rather than 1.5–2×,
the realised figure would be lower still, and it would be paid for in the one
safety net that made §11 tractable. Anyone who revisits it should start from
the table above rather than from §10's estimate, and should know that the
arithmetic is available at the price of a triangular trim and a product per
row-block — which is a different design from the one costed here.

### 11.5 What this does to the wisdom question

`thoughts.md` §10 proposed a wisdom mechanism and gave it one customer,
`Chunking::Tuned`, with the note that the case would be better made after a
target-machine run. **[C12] gives it a second customer and a better argument
than the server run would have.**

The mechanism's premise is that these choices cannot be settled by reasoning
and vary by machine. `TransformKernel` is now exactly such a choice, it is
made once at construction where a measurement is cheap against building the
table beside it, and — unlike every other knob — it has *two complete
implementations that produce the same answer*, so timing both on the actual
problem is a well-posed thing to do rather than a heuristic.

Nothing here builds it. It is recorded because §10's ordering advice —
`Chunking::Tuned` alone, and wait for the server — was written before this
existed, and the second half of that advice is now spent.

### 11.6 The target-machine run, reframed

Still wanted, no longer gating. What it settles is unchanged and is worth
restating so that dropping the gate does not read as dropping the question:

- whether the loop kernel's forward accumulator collapse at `lMax = 256` on
  eight threads (§17.5 of `field-algebra-plan.md`) gets worse at 64, which is
  the strongest single argument for this restructure;
- whether `Chunking`'s single-shared-L3 assumption survives eight CCDs, where
  the inverse's one shared copy is really one per cache domain;
- and whether the generated path, which lost on this laptop in every
  configuration, wins on a machine with sixteen times the aggregate L3.

When the machine arrives, the benchmark to run is M5's section and the
existing `server` one. Until then the laptop's numbers are real numbers about
a real machine, and two kernels that both exist can be compared on any third.

---

## 12. The wisdom mechanism, in detail

Written 2026-08-24, when `thoughts.md` §10 was picked up. That section is the
assessment — why these choices cannot be settled by reasoning, and why [C12]
gave the idea a second and better customer — and is not repeated. This is the
work order and the decisions, and the first of them is that the knob set is
smaller than §10 supposed.

### 12.1 What is actually tunable, which is three things and not six

`thoughts.md` §10 lists "five choices whose right answer is machine-dependent"
and adds `TransformKernel` as a sixth. Read against the code, most of that
list is not a tuning question, and saying which is the useful part of this
section.

| knob | where set | tunable? |
|---|---|---|
| `Chunking` | construction | **yes** — a scalar, measured at 2.0× (§10 item 1) |
| `TransformKernel` | construction | **yes** — two implementations of one answer ([C12]) |
| `thetaBlock` | per call, with a heuristic | **yes**, and smallest of the three |
| `WignerValues` | construction | **no** — see below |
| FFTW planner `Flag` | construction | **no** — FFTW already has wisdom |
| `Execution` threads | per call | **no** — the caller's, by the one-level rule |
| `RadialMajor` vs `ApplyRadially` | caller's code | **no** — not a library choice at all |

**`WignerValues` is a constraint the caller states, not a knob to optimise.**
Timing it would pick `Stored` on any machine with memory to spare, because
T11 measured `Generated` losing in every configuration — and that would be
the right answer *for time* and the wrong thing to do, because a caller who
chose `Generated` did so to avoid 648 MB at `lMax = 256` or 5.4 GB at 512. A
tuner that overrode a memory constraint on a timing argument would be
substituting its own objective for the caller's. It stays where it is.

**The planner flag is FFTW's business.** Tuning `Estimate` against `Measure`
against `Patient` is what FFTW's own wisdom does, better and with persistence
we would be duplicating. What this library could usefully add is a route to
*import and export* FFTW's wisdom, which `thoughts.md` §5 already records as
an ask against FFTWpp; that is a different feature and it is not this one.

**Threads are the caller's.** `Execution` is per call and the rule that
exactly one level threads means the caller is the only one who knows which
level that is. A tuner choosing a thread count would be choosing for code it
cannot see.

So the mechanism has **two real customers and a small third**, which is fewer
than §10 hoped and still enough: the two are the ones measured to matter most,
at 2.0× and 4.6× respectively.

### 12.2 The obstacle nobody had noticed, and what it forces

**`Chunking` and the planner flag live in `Impl`, beside the table.** So do
`WignerValues` and `TransformKernel`, and for those two it is right — they
decide what the table *is*. For the other two it is not: a chunk is a scalar
the transform reads per call, and the table does not depend on it. Sweeping
four candidate chunks today therefore means constructing four grids and
building four tables, which at `lMax = 256, nMax = 2` is 0.16 s and 648 MB
apiece — to choose an integer.

That is the trap §17.7 of `field-algebra-plan.md` fell into and reported: a
first version that allocated inside the timed loop measured 3–5× worse and
never won, and it took the numbers looking wrong to notice. Here it would not
even be a measurement error, just waste — but it is enough waste to stop the
tuner being the cheap thing §10 wants it to be.

**[C18] The table and the cheap policies are separated, and the handle carries
the cheap ones.** `Impl` keeps what decides the table — `lMax`, `nMax`, the
quadrature, `WignerValues`, `TransformKernel` and the table itself. The
handle gains `Chunking` and the planner flag as members beside the
`shared_ptr`. Then

```cpp
auto tuned = grid.With(Chunking::Fixed(8));
```

is a pointer copy and two scalars, sharing one table, and a sweep costs one
table build rather than four.

Three things follow, and the second is the one to check rather than assume.

- **`With` is offered for `Chunking` and the flag and for nothing else.** The
  two table-deciding policies have no `With`, because there is no table to
  share: changing either means a different table, which is a different grid,
  and the constructor is where you say so.
- **Two grids differing only in chunking share an `Identity()`, and that is
  correct rather than a leak.** Identity is the field layer's test that two
  operands index the same buffers, and they do: same points, same degrees,
  same table, fields interchangeable. A chunk is how the inner loop schedules
  itself and is not observable in any result — the batched tests demand *exact*
  equality against unbatched calls, which is the standing check that it is
  not. So sharing identity is the honest answer and it is also what makes
  `With` useful, since a tuned grid must stay compatible with fields already
  built on the untuned one.
- **`Impl` stays immutable**, which is what keeps a shared grid safe to use
  concurrently without a lock (step B). Nothing here adds a mutable member;
  the tuner produces values and hands them back.

*This is worth doing whether or not the tuner is ever built*, which is the
argument for taking it first. It removes a real cost from any caller sweeping
chunk widths, including the benchmark harness, which today rebuilds a grid per
chunk in its `batching` section.

### 12.3 Decisions taken here

**[C19] `Tune` returns policy values; it does not return a grid and it does
not configure one behind the caller's back.** The shape is

```cpp
auto choice = Tune<Real>(shape, wisdom);        // measures what it must
auto grid   = Grid(lMax, nMax, flag, choice.chunking,
                   WignerValues::Stored(), choice.kernel);
```

rather than `Grid::Tuned(...)`. **This is the same argument [C17] made and it
is load-bearing for the same reason.** That decision refused a silent
substitution — "honouring `Generated` and ignoring `Matrix` lies to the
benchmark, and that is fatal *specifically under [C12]*, because the entire
justification for carrying two kernels is being able to compare them". A tuner
is the piece of machinery most likely to substitute: wisdom naming a kernel
the build cannot offer, a candidate that failed to construct, a tie resolved
in favour of the incumbent. Handing the caller the values makes every one of
those visible in a variable they can print, and makes the wrong ones
impossible to hide. A `Grid::Tuned` would have had to decide each case in
silence.

It also keeps the grid's constructor where it is. That constructor already
takes six arguments and adding a seventh whose meaning is "ignore three of the
others" is the kind of interface that reads as an accident.

**[C20] The kernel comparison is opt-in, sequential, and priced.** §11.4's M5
is explicit that comparing kernels means two grids, that at `lMax = 256` both
live costs 1.3 GB, and that building them in sequence "leaves the second
starting on a cold cache with different first-touch placement". A tuner has to
choose one of those and neither is free.

It builds them **in sequence, destroying each before the next**, so peak
memory is one table rather than two. The cost is three table builds — two to
measure and one for the winner — which at `lMax = 256, nMax = 2` is about
0.6 s against 0.16 s for the untuned grid, and the fairness objection is
answered by repetition rather than by interleaving: each candidate is timed
several times and the best window taken, which is what §8's noise-floor rule
already requires and what the harness already does.

**So kernel tuning is not the default.** `Tune` measures the chunk unless
asked for more, because the chunk needs one table and the kernel needs three.
A caller who wants the kernel chosen says so, and pays a construction cost
they can see in the argument they passed.

*The alternative, recorded because it is the better answer if this ever
matters enough:* both kernels can be compared at a **smaller degree** and the
answer extrapolated. It is rejected for now because M4's own table shows the
answer changing sign with size — 4.59× at `lMax = 256, k = 8` on eight
threads against 1.20× at `k = 1` on the same row — so a proxy measurement is
exactly the kind of reasoning this mechanism exists because we cannot do.

**[C21] The fingerprint is a guard against the obvious mistake, and is not
claimed to be more.** §10 requires that "the fingerprint is part of the key",
because wisdom carried to another machine is worse than none and the failure
is silent. It is right, and the library cannot portably identify a machine:
there is no standard way to read a CPU model or a cache size, and a hostname
is not the property that matters.

So the fingerprint is what is portably available —
`std::thread::hardware_concurrency()`, the size of `Real`, the GSHTrans
version, and whether the build has a BLAS — plus, **where it can be read, a
free-text machine description recorded for a human rather than compared by the
code.** On Linux that is the model name from `/proc/cpuinfo`; where it cannot
be read the field is empty and nothing changes. An entry whose portable part
disagrees is treated as absent, not as an error, so a wisdom file moved
between machines degrades to no wisdom rather than to wrong wisdom.

**And a caller may name their own machine**, which is the escape hatch for the
case the portable part cannot see: the same binary on two nodes of a cluster
with different cache sizes. That tag joins the key. It is offered because the
alternative is a mechanism that is silently wrong in exactly the environment
this library is pointed at.

**[C22] A tie goes to the incumbent, and is recorded as a tie.** §8's noise
floor on the development machine is about ten per cent, and several of the
differences at stake are smaller. A tuner that picks the nominal winner of a
7% difference is picking noise, and worse, it will pick differently on the
next run and the caller will see the choice flapping. So a candidate must beat
the current best by more than a stated margin to displace it, the default
being the incumbent — `Chunking::Automatic()`, `TransformKernel::Loop()` — and
the entry records that the comparison was inconclusive.

That matters beyond tidiness: an entry marked inconclusive is one a later,
quieter run may usefully revisit, while an entry recording a 3% win looks like
knowledge.

**[C23] The store is a text file, versioned, and a corrupt or unreadable one
is not an error.** Text because FFTW's wisdom is text and because a file a
human can read, diff and hand-edit is the difference between a mechanism
people trust and one they work around. Not an error because losing tuning is
not losing correctness: a missing, truncated or unparseable file yields an
empty wisdom and the defaults, which is exactly what a caller who never called
`Load` gets. A throw there would turn a performance convenience into a
deployment failure.

*The one thing that is an error* is asking to `Save` somewhere unwritable,
because that is a request the caller made and can act on.

### 12.4 The shape of the key

The key is a **problem shape**, per §10, and the useful discipline is keeping
it small enough that a caller's second run hits it.

```
(lMax, nMax, precision, MRange, threads, tag)  ->  { chunking, kernel, notes }
```

Two departures from §10's list, both to make hits likelier.

**Direction is not in the key; it is in the value.** The two directions are
known to want different chunks — §10 item 1 measured 2.0× on the batched
inverse from telling them apart — so the entry carries a chunk for each rather
than the key carrying a direction. That halves the number of measurements a
caller needs to make before their entry is complete.

**Batch count is not in the key either.** `Chunking::Count` already takes the
per-field size and the number of live copies and computes a chunk, so what is
being tuned is the *cache figure that formula uses*, not a chunk for one batch
size. Tuning the figure rather than the answer means one entry serves every
`k`, which is what makes the mechanism worth having for an application whose
batch size varies by component.

*That is a small but real change to what `Chunking::Tuned` means*, and it is
better than the alternative: `Chunking::ForCache(bytes)` already exists, the
formula around it has two measured anchors (P8, and §8's corrections to it),
and a tuner that fits the one number the formula does not know is a tuner
working with the model rather than against it.

### 12.5 The steps

**W1 — [C18]'s separation.** `Chunking` and the flag move from `Impl` to the
handle; `With` is added for those two. No behaviour changes and no measurement
is expected: this is the prerequisite that makes everything below cheap.

*The tests are that a `With` grid answers identically to one constructed with
the same chunk* — exact equality, both directions, batched, since a chunk is
not observable in a result — *and that it shares `Identity()` with its parent,
so a field built on one is usable on the other.* The second is the one that
would be missed and the one §12.2 rests on.

*Done*, at 348. **One consequence this section did not price, caught by an
existing test rather than by reading:** the handle grows. It was exactly a
`shared_ptr`, and `TestGaussLegendreGrid` asserted so — `sizeof(Grid) ==
sizeof(shared_ptr<void>)`, pinning F9's "152 bytes to 16". Two more members
break that assertion while changing nothing about what it is for, since F9 is
about a copy not carrying a 648 MB table and not about a handle being one
word. The assertion is now a bound rather than an equality, and says why.

*And a third test the plan did not name*, which is that `With` really moves
the policy rather than returning a copy of the default. Without it the other
two would pass on a `With` that did nothing at all — a chunk is not
observable in a result, which is exactly what makes the no-op invisible.

**W2 — the timing core.** A function that times a candidate on the caller's
actual problem: several windows, the best taken, the spread reported, and
[C22]'s margin applied. Nothing persistent, nothing keyed.

It is separate from the benchmark harness deliberately. The harness reports to
a human and may take minutes; this runs inside a caller's start-up and must
cost a fraction of building the table beside it. What they share is the
discipline, and W2 is where the ten-per-cent rule stops being a paragraph in a
plan and becomes a constant in the library.

**W3 — `Chunking::Tuned`.** The first customer, per §10: sweep candidate cache
figures on one grid, per direction, and return the pair. One table build, and
the knob already measured to matter most.

*The acceptance test is not that it finds the optimum* — there may not be one
resolvable, and §8 records the `lMax = 128` case where the peak swapped
between runs — *but that it never returns something worse than the default by
more than the margin*, checked over several sizes. That is the property a
caller actually needs, and it is testable where "finds the best" is not.

**W4 — `Wisdom`: the key, the fingerprint, load and save.** [C21] and [C23].
Still no kernel tuning; the store's first content is W3's answers.

*A test worth writing before the code:* a wisdom file written on one
fingerprint and loaded under another yields no entry, silently. That is the
failure §10 calls "worse than none" and it should be pinned rather than
trusted to the comparison being written correctly.

**W5 — kernel tuning, opt-in.** [C20]. Sequential construction, three table
builds, and an argument the caller has to pass.

**W6 — does the mechanism pay?** The honest closing step, and the one that
decides whether W4 stays. Tuned against default, on this machine, over the
sizes the harness already walks: if a tuned grid is not measurably better than
`Chunking::Automatic()`, then the conservative default is doing its job and
the persistence layer is machinery without a customer. §10 says as much —
"if it does not, nothing has been built that has to be maintained" — and W1 to
W3 are worth having either way.

### 12.6 What this is not

Not an autotuner, per §10, and the distinction is worth keeping in the code as
well as the prose: FFTW searches a space of plans it generates, this times a
handful of named alternatives. The name is borrowed for the *persistence* —
that a machine's answer is worth writing down — and not for the search.

And the caveat §10 states, which this section does not escape: a tuner
measures which of the options we have is best on a given machine. It does not
say whether the option set is the right one.
