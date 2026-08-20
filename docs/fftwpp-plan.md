# FFTWpp: notes and suggested work

A hand-over document, on the same footing as `gaussquad-plan.md`: read from
`da380/FFTWpp@main` as fetched on 2026-08-20, with the GSHTrans side measured
rather than assumed.

## What it is, and what depends on it

About 2100 lines, header-only, across `Core.h` (an `fftw_malloc` allocator and
`FFTWpp::vector`, plan execution, cleanup), `Options.h` (flags and directions),
`Plan.h`, `Views.h`, `Wisdom.h` and `Utility.h`.

GSHTrans uses roughly a dozen names, and the load-bearing ones are
`Ranges::Layout`, `Ranges::View` and `Ranges::Plan`. Those wrap FFTW's advanced
interface — `plan_many_dft` and friends — and they are the reason the library's
batched transform is described by a `(count, stride, dist)` descriptor that
costs a descriptor rather than a repack. That was verified by running code, not
assumed: the batched FFT lands its output directly in the order the Legendre
stage wants, and both tensor layouts (component-major and point-major) are
transformable in place because of it.

**This is the dependency to keep.** Displacing it would mean rewriting the
advanced-interface handling by hand, and the design decision it enabled —
`core-plan.md` [C9], dropping the repack requirement from the tensor layer —
rests on it working.

## What was found

### The planner is not serialised, and it needs to be

FFTW's planner is documented as not re-entrant: `fftw_plan_*` may not be called
concurrently from several threads. There is **no mutex, lock or thread-related
machinery anywhere in FFTWpp** — a case-insensitive search over the whole
header set returns nothing but a comment.

GSHTrans discovered this and compensates: it carries a process-wide
`PlannerMutex` and takes it around every plan construction, while leaving
execution unlocked. That works, but it is the wrong place for it. A consumer
who plans from several threads without knowing to do this has a latent race
that will be rare, machine-dependent and very hard to attribute.

**Suggested:** the wrapper should own that lock. A `std::mutex` local to the
planning entry points, taken in `Plan`'s constructors and in `GenerateWisdom`,
would make the wrapper safe by construction and let every consumer delete their
own. It costs nothing at execution time, which is where the work is.

### There is no support for FFTW's own threading

No `fftw_init_threads`, no `fftw_plan_with_nthreads`, no
`fftw_cleanup_threads`. A consumer linking `libfftw3_omp` or `libfftw3_threads`
gets no way to use it through this wrapper.

For GSHTrans this is not a problem and arguably a virtue: the library threads
over colatitudes itself and each thread executes its own plan on its own
buffers, so FFTW-internal threading would be nested parallelism and unwanted.
But that is a policy this library happens to share, not a general one, and a
consumer with one very large transform has no route.

**Suggested:** either add the three calls behind an opt-in, or state in the
documentation that FFTW-internal threading is out of scope and the caller
should parallelise outside. Silence is the only bad option, because it looks
like an oversight.

### Wisdom is better supported than GSHTrans's own notes claimed

`Wisdom.h` provides `ImportWisdom`, `ExportWisdom`, `ForgetWisdom` and two
`GenerateWisdom` overloads taking layouts. An earlier GSHTrans note said there
was no route to persistent wisdom; that was wrong, and the note has been
corrected. GSHTrans does not use them because step E of its own plan removed
the `WisdomOnly` policy — the grid used to pre-generate wisdom for exactly two
transform shapes and then refuse to plan anything else, which made every
batched shape unplannable.

**Nothing to do here**, beyond perhaps a worked example: import at start-up,
export at exit, and what happens when a plan shape is absent from the wisdom.

### Alignment is exposed but not queryable

`Options.h` exposes `Unaligned` (`FFTW_UNALIGNED`), so a consumer can plan for
storage whose alignment it does not control. What is missing is a query —
FFTW's `fftw_alignment_of` — and any checked form of the new-array execute.

`Plan.h`'s own comment states the constraint: new-array execution is valid only
for buffers with the same layout *and alignment characteristics* as the
planning buffers, and neither FFTW nor FFTWpp checks it. GSHTrans treats this
as a hazard rather than a feature: it copies every row through the plan's own
buffers in both directions rather than executing on caller storage, which was
measured to be invisible against the Legendre stage.

**Suggested, in order of usefulness:** expose `fftw_alignment_of`; then offer a
checked new-array `Execute` that throws or asserts when the alignment classes
differ. That turns silent undefined behaviour into a diagnosable error, and it
would let a consumer skip a copy when it is safe to do so and know when it is
not.

### There are no install rules

As with GaussQuad: no `install()` and no `export()`, so `find_package(FFTWpp)`
cannot work, and a downstream `INTERFACE` target that links FFTWpp cannot be
exported. **This and the same gap in GaussQuad are jointly what stop GSHTrans
shipping a CMake package**; it currently installs its headers and no config,
and consumption is by `add_subdirectory` or `FetchContent`.

## Suggested work

1. **Install and export rules.** Same twenty lines as GaussQuad, same payoff.
2. **Own the planner mutex.** Correctness, and it lets consumers delete theirs.
3. **Expose `fftw_alignment_of`, and a checked new-array execute.** Turns a
   documented piece of undefined behaviour into something diagnosable.
4. **Decide about FFTW-internal threading** — support it or say it is out of
   scope.
5. **Nothing about the Ranges interface.** It is the part GSHTrans leans on
   hardest and it has not needed changing.

## A note on the wider build

Both of these, and GSHTrans itself, currently lean on "fetch everything".
`FetchContent` is the right default for getting started and the wrong one for
being depended on: it makes every consumer re-download and rebuild, it prevents
a distribution or an environment module from supplying the library, and — as
above — it makes proper packaging impossible.

The cheap improvement, applied consistently, is to declare dependencies with
`FIND_PACKAGE_ARGS` so that `FetchContent_MakeAvailable` tries `find_package`
first and only fetches when nothing is installed. GSHTrans now does this. It
costs one line per dependency and means an installed copy is used when one
exists, which is what makes the install rules above worth having.
