# Lessons

Things that cost real time once and would cost it again. Kept short on
purpose: this is not a design record, and it is not a changelog. The design
lives in the code and its comments; the mathematics lives in
`canonical-components.tex`; what the library does lives in
`gshtrans-reference.tex`.

## Measurement

**Benchmark from a Release tree, and check.** `cmake -S . -B build` leaves
`CMAKE_BUILD_TYPE` empty, and the harness there runs about ten times slow with
every figure self-consistent — it once reported speedups of forty. The harness
warns now, but the habit is the real guard.

**Interleave the two versions in one session.** A laptop's clock scales under
load: repeating one measurement seconds later moves it by tens of per cent,
which is enough to invent a difference that is not there. Measured, the first
pass over a fixed grid runs about nine per cent slow against the sixth,
systematically against whichever candidate goes first. Timing A to completion
and then B is what this forbids.

**A handle copied in an inner loop is a contended atomic.** The grid is a
`shared_ptr` handle so that fields and views can carry it by value, and copying
one is an atomic increment and decrement on a count *every* holder shares.
`TensorExpansion::Coefficient` built a view — and so copied the handle — on
every `(l, m)` read, and every spectral operator reads through it. Alone that
was 8 ns a call. On eight threads, each working on its *own* expansion over
the one grid, it was 380 ns, and eight surface gradients took four times
longer on eight threads than on one: no data was shared, only the count.
Reading the buffer directly gave 1.1 ns and 3 ns, and the eight gradients went
from 0.48 s to 0.0096 s. The rule is the one the radial seam already states
for allocation: a handle is copied when an object is made, never when an
element is read.

**A sequential kernel's time moves by ten per cent with where the code
lands.** Checking that a change to the OpenMP regions cost nothing, the
sequential inverse came out 11 per cent slower, reproducibly, in code the
change had not touched. Reverting half the change moved the loss to the
forward transform instead; forcing the row lambdas inline made both worse.
Then the *old* headers were rebuilt with a no-op edit to the benchmark's
`main`, and they moved by the same amount — 61 to 69 ms one way, 63 to 70 the
other. It is the placement of a hot loop in the binary, so it is a property of
the build and not of the run: interleaving does not average it away, and two
builds that differ anywhere cannot be told apart below about a tenth on those
rows. The threaded rows and the matrix kernel do not show it. Before believing
a small sequential difference, perturb the *baseline* and see how far it moves
by itself.

**Require a margin.** Ten per cent is the noise floor. A tuner that picks the
nominal winner of a seven per cent difference is picking noise, and will pick
differently next run.

**One OpenMP runtime per process.** The OpenMP build of OpenBLAS links
`libgomp`, so a clang build of this library — which uses `libomp` — puts two
runtimes in one process. Results stay correct; timings do not, because the
library holds a same-runtime BLAS to one thread through that runtime, and a
BLAS on another runtime is out of its reach.

**A parallel region with one thread is not a parallel region.** The BLAS was
first kept serial by issuing every GEMM from inside an OpenMP region, on the
argument that a region opened inside a region is nested and the default of one
active level makes a nested region serial. True of a team of two or more. A
team of *one* is inactive: it is not a level, `omp_in_parallel()` is false
inside it, and a region opened from it gets the whole machine — so the
sequential case, the commonest, was the one case in which the BLAS was not
held back. Measured against the OpenMP build of OpenBLAS at `lMax = 512` with
a chunk of eight, a transform under `Execution::Sequential()` used 4.7 cores
and ran 14 per cent *slower* than it does on one. It went unseen twice over:
the development machine's default BLAS is the pthread build, which none of
this touches, and OpenBLAS threads a GEMM only above a size that the default
chunking at `lMax = 256` stays under — so it waits for a large cache, a large
degree, or a tuned chunk, which is to say for the production machine. What
works is `omp_set_num_threads(1)` inside the
region: that sets the thread count of the tasks that go on to open nested
regions, and of nothing else. `Details::InSerialisingRegion` is that, and
`Threading.ARegionOpenedInsideTheSerialisingRegionGetsOneThread` is the test.

## Numerics

**Build tangential test fields as `grad f + r̂ × grad g`.** A field written
`a(θ,φ) θ̂` is almost never smooth, because `θ̂` is not, at the poles, so its
spin-weighted expansion does not converge. A divergence check built that way
was wrong by 0.18 at *every* truncation from `lMax` 8 to 128, which reads as a
bug and is not.

**Under Schulten–Gordon the completeness relation tests nothing**, because the
algorithm normalises by it: it is tautological per row and blind to a bad
match. The runtime check is the recurrence residual.

**Track the 3-j phase; do not read it back.** A near-stretched row at high
degree spans `1e201`, the rescaling flushes the last element to zero, and
whole rows come back negated with every magnitude correct to `1e-16`.
Completeness and the residual are both blind to it — only
cyclic-permutation invariance catches it. SLATEC has the same fragility.

**The chunk heuristic's `copies` argument is not the thread count**, in either
direction. The forward transform gives every thread a private accumulator, so
its copies are its threads; the inverse gathers one block, shared and
read-only, so it has one copy however many threads read it. Serving both with
the thread count starves the inverse, which measured 2.2× slow where the whole
batch would have fit.

**Polar truncation was measured and dropped.** It came out at about 1.15×
against the 1.5–2× assumed. Do not re-propose it without measuring again: it
shrinks the GEMM's inner dimension and efficiency falls with it, so the
realised gain is lower still.

**The Wigner recursion has a ceiling, and it is a formula.** The seed of the
column at order `m` is about `(sin θ)^m`, and the column recovers only once
`l sin θ` reaches `m`, so a seed that underflows *and* matters needs
`lMax · s|ln s| > |ln min|`, and `s|ln s|` peaks at `1/e`: the recursion is
sound iff `lMax < e·|ln min|` — 1926 in double, 237 in single, and measured to
be so within a per cent. Above it columns come back zero where the truth is
order one, with no sign of trouble: the benchmark's `huge` section ran at 2048
for months. `MaxSafeDegree` enforces it, a little low so the seed stays out of
the denormals. The cure, an exponent carried per column, is known and unbuilt.

**Watch the quantity, not a proxy for it.** Schulten–Gordon recurses inward
from both ends and must stop each half where its values stop growing. SLATEC
stops where the recurrence coefficient `|c1|` first rises, and so did this
library. That coincides with the values' behaviour when a row has an allowed
region and not when it is near-stretched, and there one half ran downhill:
four digits lost where it cost least, and 44 of 77 `(l, 2l, l)` tables to
`l = 1000` that could not be built at all. Stopping at the first fall of `|g|`
itself fixed both at no cost, and was better on fat triangles too.

**Single precision is for storage.** A recursion's range is the range of its
arithmetic, and a three-term recurrence loses about `n²ε` along a row however
good the algorithm. So `float` Wigner tables and `float` 3-j tables are
computed in double and rounded. It costs construction time and buys the whole
of double's range and accuracy; single precision then pays where it should, in
what is stored and moved.

## C++

**Sizes and indices are signed.** `std::ptrdiff_t` throughout, and
`std::size_t` only at the point where a standard container is indexed or
sized. A size here is something subtracted from, compared with an order that
may be negative, and multiplied into an offset, and unsigned arithmetic gets
each of those wrong in silence: the difference that should be negative is
enormous instead, and no sanitiser objects, because wrapping is defined
behaviour. Three accessors — `NumberOfCoLatitudes`, `NumberOfLongitudes`,
`FieldSize` — had inherited `size_t` from a `.size()`, and the price was
seventy-odd `static_cast<Int>` at their call sites, each a place where a
negative value would have become a large one. Found by turning
`-Wsign-compare` on in the tests.

**A rule written four times drifts three ways.** "Component α from its orbit
representative" was written out in the flat and layered tensor types and in
both their expansions, and each copy came to be wrong differently: a sign on
real orbits, the conjugate of `i` on imaginary ones, a dropped `i`, a
dangling view. None showed at rank 2, where no pinned member is reached by
conjugating, so every test passed while a real elastic tensor read back with
the wrong sign. It is one pair of functions now, beside the table they
interpret. If a rule has to exist in two places, the test belongs on the
*relation it implements*, over every case, and not on either copy.

**A letter outside a `MultiIndex` alphabet is a hard error**, thrown inside a
constant expression, not a SFINAE-friendly constraint failure. Over
`TangentialSlots`, `requires { Flat<0, 1>; }` is **true** and the use then
fails to compile. Every accessor taking a component's letters as template
arguments has to check them first, in an `if constexpr`.

**`and` short-circuits evaluation, not well-formedness.** A `static constexpr
bool` written as a chain whose later terms name members that the earlier terms
guarantee compiles under GCC and is rejected by clang. Write a concept
instead. This is the class of divergence the clang leg of CI exists for.

**Doxygen attaches a `///<` at the top of a function body to the enclosing
function.** Four transform entry points were "documented" for a while by a
stray comment on a function-local type alias, and the warning log was clean
throughout.

**`#pragma once`, and no guards.** A guard is a name that has to be unique and
is chosen by hand; this tree has two files called `BundleMaps.hpp`, and only
the care taken over forty macro names kept them apart. The pragma fails only
when one file is reached by two paths, which a header-only library included
through one root does not do.

**To make OpenMP optional, ask the compiler and leave the pragmas alone.**
`_OPENMP` is defined exactly when the flag is given, so a macro of the
library's own could only ever disagree with it. Two things were learned on the
way to `OpenMP.hpp`. Serial stand-ins *named* `omp_get_max_threads` and so on
collide, because GCC ships `<omp.h>` whether or not OpenMP is on and a caller
may include it: the wrappers need names of their own. And hiding the pragmas
behind a macro works and cannot be kept, because no formatter can lay out a
macro whose argument says `for` and `if`; an ignored pragma with
`-Wunknown-pragmas` silenced for the length of the header, and only when
`_OPENMP` is undefined, costs the caller nothing.

## Build and CI

**`cmake --build --parallel` with no number is `make -j`, which is no limit.**
Under the Makefile generator every translation unit starts at once, and these
are heavy ones. The development laptop reached its critical temperature and
powered itself off in the middle of a gate run at sixteen jobs, straight after
two other full builds; six jobs still touched 88 °C. Give the number, and on
the laptop give a small one.

**Build the no-dependency configuration before pushing anything that touches
`SphericalGrid.hpp`.** A private member added inside its
`#ifdef GSHTRANS_HAVE_BLAS` block compiles in the default tree and fails only
where BLAS is absent.

**`Blas.hpp` assumes a 32-bit-integer BLAS.** An ILP64 build — MKL's, for
instance — uses the same symbol names, so linking one would not fail to link:
it would pass the wrong thing. There is no portable way to detect it from
inside the header.

**A test that cannot fail is worse than none**, because it is counted. The
2026-09 review found five: a nesting check evaluated after the calls it was
about had returned; a "relative" error divided by a `bool`; an elastic-tensor
test comparing an expression with itself; a residual tolerance a row wrong by
a third could pass; and a CI matrix in which the BLAS kernel's tests compiled
to nothing in every job, all green. The habit that catches them is to watch a
new test fail before trusting it to pass — two written during that work passed
on the broken code at first, and had to be sharpened.

**Before "fixing" an API oddity, look for the test that pins it.** Three
low-severity findings of the same review were decisions and not oversights,
and each had a `static_assert` saying so: no floating-point scalar of another
precision, `Materialise` defaulting to complex, and `conclusive` meaning that
the incumbent was displaced.

## Open questions

Deliberately unresolved, each wanting a decision rather than an edit:

- **`Scheme`** is shaped unlike the other policy values — a class with no
  instances whose named constructors return nested tag types, because the
  three schemes have different interpolant types and the dispatch has to be
  overload resolution. Either rename it or accept the exception.
- **`Chunking::Count`** takes `int copies` where everything around it uses
  `std::ptrdiff_t`.
- **`SphericalGrid.hpp`** is the largest file by a wide margin.
- **The two layered field types** duplicate a good deal between them.

