# FFTWpp: where things stand, and what's next

Written 2026-08-21, at the end of the work driven by `fftwpp-plan.md`.

## Done

Six commits on `develop`, all pushed, CI green on every job. `git log
main..develop` has the detail; each commit message says why, not just what.

Every item in `fftwpp-plan.md` is addressed:

1. **Install and export rules.** `find_package(FFTWpp)` works and yields
   `FFTWpp::FFTWpp`. `FindFFTW.cmake` is rewritten around imported targets
   with components, and ships inside the install so consumers need nothing
   extra. Dependencies use `FIND_PACKAGE_ARGS`. `tests/package` is a
   standalone project that consumes the installed package, and CI runs it.
2. **The planner mutex is ours.** `PlannerMutex()` / `PlannerLock` guard every
   planner entry point in `Core.h`. Taken at the leaf only, so the higher
   layers cannot deadlock against it. Execution stays unlocked.
3. **Alignment.** `AlignmentOf`, `SameAlignment`, and on `Plan`:
   `CanExecuteOn`, `ExecuteChecked`, `InputAlignment`, `OutputAlignment`.
4. **FFTW-internal threading is opt-in**, behind `FFTWPP_USE_FFTW_THREADS`,
   with `ThreadSession` as the RAII form and `ThreadsEnabled` as the
   compile-time query.
5. **The Ranges interface is unchanged in signature.** Additions only.

Plus: CI on GitHub Actions (GCC and Clang, Debug and Release, `-Werror`,
macOS, a threads build, install-and-consume, ASan, TSan, clang-format,
Doxygen); the test suite from 31 cases to 99; the examples now self-check and
run as tests; README and `docs/testing.md` rewritten.

Since then, `CleanUp()` has been demoted. It was shown as the way to finish a
program in every example, the README and the mainpage, and it is not: FFTW's
persistent state is reachable for the life of the process, so leaving it is not
a leak, while calling it discards accumulated wisdom and leaves every live plan
undefined -- including plans owned by unrelated code in the same process.
`Ranges::Plan` now keeps an atomic count, `LivePlanCount()` reports it, and
`CleanUp()` and `CleanUpThreads()` throw rather than do it silently.

Two things found along the way that were not in the plan:

- **The library could not compile against libc++ at all**, because of
  `ranges::fold_left_first` and `views::zip_transform` — libc++ does not ship
  those until LLVM 22 and 23. macOS was excluded outright. Six call sites
  replaced; nothing else needed C++23, so the requirement dropped to
  `cxx_std_20`. The macOS CI job is what keeps it there.
- **`Normalisation()` was wrong for an R2C plan** — it used the output
  dimensions, which are the halfcomplex `n / 2 + 1`. Only ever correct because
  callers happened to ask the backward plan. A test now requires both
  directions of one logical transform to agree.

## Examples

Renamed to the number-name form on 2026-08-21 and extended from five to
eight; `06-guru_layouts`, `07-plan_reuse` and `08-threads` are new. All are
registered as tests, so they are run rather than merely compiled.

Writing `07-plan_reuse` turned up a false negative: `CanExecuteOn` compared
alignment classes even for a plan created with `Unaligned`, which assumes
nothing about alignment, and so refused buffers the plan could run on. Fixed
in both `Plan` and `GuruPlan`.

## Loose ends

- GSHTrans still carries its own process-wide `PlannerMutex`. It is now
  redundant and can be deleted there. That is a change in the other
  repository, not this one.
- `develop` has not been merged to `main`.
- `NumericConcepts` is still fetched at `GIT_TAG main`. Pinning it to a tag
  would make builds reproducible, but that needs a tag to exist upstream.

## The guru interface: built

Implemented on 2026-08-21, along the lines assessed below. `FFTWpp::Dim`,
`Ranges::GuruLayout`, `Ranges::GuruPlan` and the `TransformAlong` builders
live in `FFTWpp/src/Guru.h`; the raw precision-dispatching wrappers are in
`Core.h`. Split-complex is left out, as recommended.

Two things worth remembering from doing it:

- `GuruPlan` is constrained to *views*, not ranges, and that is load-bearing.
  When a container could be deduced, the by-value constructor parameter took a
  copy, the plan was built on the copy's storage, and every transform wrote
  where the caller could not see it. Every round trip still passed, because a
  copy round-trips perfectly well. Only the cross-check against the same
  transform driven by hand caught it.
- Introducing a `Ranges::Internal` namespace shadowed `FFTWpp::Internal`, so
  an unqualified `Internal::` inside `Ranges` resolved to whichever happened
  to be visible. There is now one namespace of that name; keep it that way.

The assessment below is kept because it records why the scope is what it is.

## The guru interface: the original assessment

Asked whether it is worth exposing. Short answer: doable, worthwhile, but
worth scoping much more narrowly than "expose the guru interface" — and the
valuable part is not the intimidating part.

### What it actually adds

Guru generalises the advanced interface in two ways: each dimension gets its
own input *and* output stride, and the `howmany` loop becomes
multi-dimensional rather than a single `(count, dist)` pair.

The second is the one that bites. Take a 3D array `(n0, n1, n2)`, row-major,
and transform along the middle axis. The transforms start at offsets
`i0 * n1 * n2 + i2` — two independent strides, and `Layout` has one `dist`.
Not expressible. Along axis 0 or axis 2 it is fine; along any interior axis it
is not. Guru does it in one plan with `howmany_rank = 2`:

```
dims        = {{n1, n2, n2}}
howmanyDims = {{n0, n1 * n2, n1 * n2}, {n2, 1, 1}}
```

"FFT along axis k of an N-dimensional array" is the commonest thing people
want, and it is exactly where the current wrapper runs out. That is the case
for doing it.

### But the workaround just got safer

That middle-axis transform can already be done as `n0` new-array executions of
a rank-1, `howMany = n2` plan at offsets. As of this work that path is
*checkable* rather than silently undefined: `ExecuteChecked` reports when an
offset lands in a different alignment class. So guru is a performance and
elegance win here, not a capability gap. FFTW can plan the whole nest better
than a hand-driven loop can, but the gap is smaller than it first looks.

### What to build, in order

1. A `Ranges::Dim{n, inStride, outStride}` and a guru plan constructor taking
   `dims` and `howManyDims`. Roughly 150 lines of `Core.h` wrappers plus a
   constructor and validation in `Plan.h`. Purely additive; no API break.
2. **The part that actually matters:** a helper that *computes* the dims from
   a shape and a list of axes, so nobody hand-derives strides — something like
   `Ranges::TransformAlong(shape, {1})`. The guru interface is intimidating
   because you do stride arithmetic in your head; a library that does it for
   you removes the entire reason it is intimidating. If only one thing gets
   built, build this.

### Mirroring the format type-safely

Asked separately whether guru's descriptor could be mirrored in a type-safe
way. It can, and the type safety comes from somewhere slightly unexpected.

`fftw_iodim` is three bare ints, and its failure mode is transposing the input
and output strides, or swapping the two arrays. Strong typedefs per field would
be heavy; C++20 designated initialisers do the job for almost nothing:

```cpp
struct Dim {                      // mirrors fftw_iodim
  std::ptrdiff_t n;
  std::ptrdiff_t inStride;
  std::ptrdiff_t outStride;
};
```

```c
/* C */    fftw_iodim dims[]    = {{n1, n2, n2}};
           fftw_iodim howmany[] = {{n0, n1*n2, n1*n2}, {n2, 1, 1}};
```
```cpp
/* C++ */  auto layout = Ranges::GuruLayout{
               .transform = {{.n = n1, .inStride = n2, .outStride = n2}},
               .batch = {{.n = n0, .inStride = n1*n2, .outStride = n1*n2},
                         {.n = n2, .inStride = 1,     .outStride = 1}}};
```

Three things fall out of that shape:

- **The two dim lists become unswappable**, by living as named members of one
  `GuruLayout` rather than as two same-typed arguments. That is the mistake
  types can actually catch here; the stride values they cannot.
- **`ptrdiff_t` hides the guru/guru64 split.** Dispatch to
  `fftw_plan_guru64_*` when a value exceeds `INT_MAX`, else `guru`. The C API
  makes that the caller's problem; a wrapper should not.
- **The descriptor belongs to the plan, not the view.** `fftw_iodim` couples
  both sides in one struct, so it cannot hang off a one-sided `Layout`.
  `GuruPlan(inRange, outRange, guruLayout, flag, direction)` is the honest
  mirror, and it leaves `Layout`/`View` untouched as the ergonomic default.

The real type safety is not in the struct at all, though. It is in never
writing a stride -- see `TransformAlong` above. Once that exists `Dim` is an
implementation detail most users never see.

### What to leave out

Split-complex arrays (separate real and imaginary storage). They do not fit a
`std::complex` range at all — different value-type story — and there is no
consumer. Half-doing it would be worse than not doing it.

### Decide up front

Whether the dim type carries `int` or `ptrdiff_t`. `guru64` is nearly free if
the choice is made at the start and painful to retrofit.

### The real risk

Guru lets you write layouts where two indices alias the same address, which is
undefined and silently wrong rather than a null plan. A full overlap check is
a lattice problem, but a cheap sufficient one catches the realistic mistakes:
sort dims by `|stride|` and require each stride to be at least the span of
everything inner to it. That is `O(rank log rank)`, and should be documented
as a necessary condition rather than a proof.

This is an argument *for* wrapping guru rather than against it. It is the
place where the validation culture added in this round pays off most.

### Recommendation as given

Build it when something pulls on it, not before. GSHTrans is well served by
the advanced interface, per `fftwpp-plan.md`, so there is no consumer waiting.
An untested generalisation with no user is a liability, and guru is exactly
the surface where untested means wrong.

The exception is item 2 above: `TransformAlong` is useful on its own merits,
makes a common thing easy, and forces the descriptor design to be settled.
That is the piece worth doing first, and it can be done for the advanced
interface's expressible cases before any guru code exists at all.

In the event the whole thing was built, minimally, on the grounds that it can
be refined if use cases arise.
