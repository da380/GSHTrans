# GSHTrans: design

How the library is put together, and why it is put together that way. This
document is about the software. The mathematics, the conventions and the
numerical methods are in `gshtrans-reference.tex` (and its PDF); how to use the
library is in `README.md` and the numbered `examples/`; what is still open is
in `open-issues.md`.

The intended use sets the priorities throughout: real rank-2 and rank-4 tensor
fields, layered Earth models, Release builds, and a many-core shared-memory
target (dual-socket, 64 cores, large L3). Where a choice trades convenience for
"never return a wrong number silently", the library takes the second.

---

## 1. Layers

The library is header-only. `GSHTrans/GSHTrans.hpp` includes everything;
`GSHTrans/Core.hpp` includes the core alone. The extensionless `GSHTrans/All`
and `GSHTrans/Core` forward to those two for code written against an older
layout.

| Layer | Headers | What it provides |
| :--- | :--- | :--- |
| Core | `SphericalGrid.hpp`, `GaussLegendreGrid.hpp`, `Wigner.hpp`, `WignerMatrices.hpp`, `Indexing.hpp`, `Views.hpp`, `Policies.hpp`, `Tuning.hpp`, `Blas.hpp`, `OpenMP.hpp`, `Utility.hpp`, `Concepts.hpp` | quadrature grid, Wigner tables, the batched transform and its two Legendre kernels, policies, threading, measurement-based tuning |
| 3-j | `3j.hpp` | Wigner 3-j symbols and coupling matrices; independent of grids and fields |
| Spin fields | `SpinField/` | the lazy, index-checked algebra of spin-weighted fields |
| Tensors | `Tensor/` | tensor fields in the canonical basis, symmetry orbits, tensor algebra |
| Expansions | `Expansion/` | the spectral side: spin and tensor expansions, `ð`/`ð̄`, raising and lowering, surface gradient, intrinsic derivative, interpolation |
| Layered | `Layered/` | fields on a stack of spheres; radial grids with interfaces; radial operators; `RadialMajor` |

The 3-j header stands alone: it includes nothing from the rest of the library.

---

## 2. Grids and transforms

### 2.1 The grid is a handle

`GaussLegendreGrid` (a `SphericalGrid` with Gauss–Legendre colatitudes) is a
value-semantic handle to a shared, immutable implementation that owns the
quadrature, the Wigner table and a per-thread cache of FFTW plans and work
buffers. Fields, views and expansions carry the handle by value, so they can be
returned and stored without lifetime rules, and two grids are the same grid
when `Identity()` agrees.

`grid.With(...)` changes the chunking policy or planner flag without
rebuilding the table and keeps the identity, so fields are interchangeable
between a grid and its retuned copy.

**A handle is copied when an object is made, never when an element is read.**
Copying a `shared_ptr` is an atomic increment and decrement on a count every
holder shares. A per-element copy costs a few nanoseconds on one thread and
hundreds on many, because every thread contends on that one count even when no
data is shared. Accessors such as `TensorExpansion::Coefficient` therefore
index the underlying buffer directly rather than building a view per read.

`ReleaseThreadCaches()` frees the calling thread's plans and buffers. The cache
otherwise lives as long as the thread, which matters to a caller that wants to
call `FFTWpp::CleanUp()`.

### 2.2 The primitive is a batch

Every transform is of **a batch of fields that share grid, degree and upper
index**; the single field is a batch of one. The Wigner block for a given upper
index is what batching amortises, so fields of different spin are separate
calls — `TensorField` is that loop. There is no "mixed upper index" batch
because the kernels would have nothing to share across it.

A `Batch` describes where the fields are, independently for input and output:

- `Contiguous`, `Interleaved` and `Strided` cover every layout that is affine
  in the field index (element `j` of field `k` at `j * stride + k * dist`).
- `Batch::At(offsets, stride, size)` places each field at an arbitrary offset,
  and `Subset(which, size)` selects fields of an existing batch (for instance
  the solid radii of a layered field).

The kernels touch a layout only through `Count`, `Stride`, `Offset`, `Span` and
`Disjoint`, and always gather into chunk-local scratch before FFTW or the BLAS
sees anything. That is why an arbitrary layout costs nothing extra: no kernel
knows about it. Decisions behind `At`:

- **The batch owns its offset table**, shared and immutable. A `Batch` stays a
  value that can be copied and returned, and a dangling table is impossible.
  The alternative (a non-owning span, as in FFTW's guru interface) would make
  a stale table a transform that *writes* to the wrong place in silence.
- **Disjointness is checked once**, at construction, by sorting; this is why
  `At` takes the field size. A call then only checks that the fields it is
  given fit.
- **Equality is by value**, so two batches that describe the same layout are
  equal whichever table they hold.

An FFTW-style plan object (validate once, allocate outside the parallel
region, hold tuning results) was considered and not built: it cuts across
keeping policies on the grid handle with `With()`, and `At`'s construct-time
validation captures most of its benefit.

### 2.3 Two Legendre kernels

The transform is an FFT in longitude and a Legendre sum in colatitude. The
Legendre stage has two implementations, chosen when the grid is made:

- `TransformKernel::Loop()` — the default — interleaves the two stages per
  colatitude block.
- `TransformKernel::Matrix()` does every FFT first and then one pair of GEMMs
  per order against a table laid out `[n][m][l][θ]`. It needs a BLAS, stores
  only half the Wigner values (the other half by the reflection
  `d^l_{nm}(π − θ) = (−1)^{l+n} d^l_{n,−m}(θ)`, so it requires a symmetric
  grid), and has no accumulator or reduction in either direction because
  different orders write disjoint output.

Both are kept permanently. The loop kernel is the matrix kernel's oracle — the
same sum through an independent arrangement — and which one is faster depends
on the machine, the degree and the batch size, which only a measurement on
that machine can settle (`Tuning.hpp`, `TransformBenchmark kernels`). The
derivation of the matrix form and its measured gains are in the reference.

### 2.4 Chunking

A batch is processed in chunks sized to fit cache. `Chunking::Count` takes the
bytes per field and the number of **copies** of the working set that will be
live at once. That is not the thread count: the forward transform gives every
thread a private accumulator (copies = threads), while the inverse shares one
read-only block (copies = 1). Using the thread count for both starves the
inverse.

`Chunking::MaximumCount` is 63, not 64: a chunk of 64 lands on one of the
power-of-two FFT strides that cause cache-set collisions, and measures 8–20 %
slower at small degree.

### 2.5 Wigner values: stored or generated

`WignerValues` chooses at construction whether the grid holds the table or
generates each block inside the transform. A stored table and a generated
one agree **bit for bit**, which the tests check; this is why both routes
compute their values the same way, including in single precision (§5.2).

---

## 3. Threading and errors

### 3.1 One level, and only when asked

Threading is a per-call `Execution` policy that defaults to sequential: the
library creates no threads unless asked. Exactly one level threads. A call
asked to run in parallel from inside an existing parallel region runs
sequentially; the rule lives in one place, `Execution::TeamSize`.

OpenMP is optional. `OpenMP.hpp` is the only file that includes `<omp.h>`; it
wraps the few runtime calls the library makes under names of its own (GCC
ships `<omp.h>` whether or not OpenMP is enabled, so stand-ins named
`omp_get_max_threads` could collide with a caller's include). The pragmas stay
as pragmas: the headers that use them silence `-Wunknown-pragmas` for their own
length when `_OPENMP` is undefined. The library asks the compiler through
`_OPENMP` rather than defining a macro of its own that could disagree with it.
`OpenMP::Available` reports which was built.

### 3.2 Keeping the BLAS serial

The GEMMs here are skinny, and threading them inside a threaded transform
loses. A BLAS on the same OpenMP runtime is held to one thread by issuing every
GEMM inside `Details::InSerialisingRegion`, which sets the nested thread count
to one. Merely being inside a parallel region is not enough: a team of one is
*inactive*, so a region opened from it is not nested and gets the whole
machine — and the sequential case is the commonest. A BLAS with its own thread
pool, or on another OpenMP runtime (e.g. a clang build against a `libgomp`
OpenBLAS), is out of the library's reach and must be set to one thread by the
caller.

`Blas.hpp` assumes 32-bit BLAS integers. An ILP64 BLAS uses the same symbol
names, so linking one would not fail; it would pass the wrong arguments. There
is no portable way to detect this from a header.

### 3.3 Preconditions throw, in every build

Size, degree, upper-index, layout and range checks throw
`std::invalid_argument` in Release as well as Debug.
The library is meant to be run in Release, and an `assert`-only check there is
no check. `assert` is kept only for internal invariants that no caller input
can reach, and where an asserting member initialiser would run before a
throwing check, the throwing check is placed first so that both builds throw
the same thing.

### 3.4 Exceptions and OpenMP regions

An exception that escapes an OpenMP region terminates the process. Every
region that can throw — both Fourier stages, both kernels, the loop over
orders, the Wigner table constructions, the radial-line loops — runs its body
through `Details::ExceptionCapture` (`Utility.hpp`), which keeps the first
exception and rethrows it after the region closes. A thread whose set-up
failed still passes through the worksharing loop and does nothing, as OpenMP
requires. The result is that a call that throws sequentially throws the same
exception threaded. Regions whose bodies only copy or allocate nothing are not
wrapped.

---

## 4. The field algebra

### 4.1 The upper index is in the type

A spin-weighted field carries its upper index `N` as a template parameter, so
every rule of the index algebra (the table in `README.md`) is checked when an
expression is written. An unlawful combination fails overload resolution at
the call site, which also makes it testable: `static_assert(!requires { u + v; })`.

Rules that are facts rather than table rows:

- **Conjugation reverses the upper index.**
- **`real`, `imag` and `RealValued` exist only at `N = 0`**, the only index at
  which real-valuedness is invariant under rotation of the local frame. The
  rule is closed under every node, so it forbids nothing but a real-valued
  terminal at nonzero spin.
- **One precision per expression tree.** A floating-point scalar of a different
  precision is refused, because `2.0 * floatField` would narrow silently. An
  integer scalar is accepted: it is exact in every precision.
- `pow`, `exp` and the like are `Map(f, F)` at `N = 0`, with `F`'s result type
  constrained so that a wrong callable is a missing overload, not an error deep
  inside a node.

### 4.2 Lazy, pointwise, alias-safe

Operators build expression nodes; nothing is evaluated until assignment,
`Materialise` or `Integrate`. Every node is pointwise and index-preserving —
the value at `(iTheta, iPhi)` reads its operands only there — so `u = expr`
with `u` inside `expr` is safe without a temporary. Operations that are not
pointwise (derivatives, raising and lowering, radial operators) live in the
spectral and layered layers, where they produce new storage.

Lifetimes inside expressions:

- A node holds an owning operand by moving it in, and a **view by value**, just
  as it holds another expression. A view held by reference would dangle when
  an expression of named views is returned from a function.
- `Component()` is `const&`-qualified and deleted on rvalues — always on
  `TensorField`, which owns, and on a node only when it has taken ownership of
  a field (`TensorDetails::HoldsStorage`). So `Transpose(t).Component<0,1>()`
  over a named `t` works, and the same call on a temporary field does not
  compile instead of dangling.
- `Materialise` defaults to `ComplexTensor`: an expression node carries no
  reality to follow, and complex is never wrong.

### 4.3 Tensors: canonical components and orbits

A tensor field is stored by its canonical components in the basis
`e_-, e_0, e_+`; a component with lower indices summing to `N` is a
spin-weight-`N` field. Symmetries (permutation symmetry, reality) partition the
components into **orbits**; only one representative per orbit is stored — the
member of smallest flat index, whose first non-zero slot is `-` when reality
is in the group. Every other component is derived from its representative by a
sign and possibly a conjugation; an orbit that reality pins to real or
imaginary values is stored in a separate real buffer.

That derivation is written **once**, as `DerivedComponent` and
`DerivedCoefficient` beside the orbit table in `Tensor/Orbits.hpp`, and used by
the flat and layered types in both domains. The relation, not any one copy of
it, is what the tests check: every component of every tensor type obeys every
generator and the reality condition, spatially and spectrally. A rule written
in several places drifts apart, and does so first in the cases (real tensors of
rank three and above) that low-rank tests never reach.

Consequences:

- A real tensor cannot live on an `NRange = NonNegative` grid: some of its
  representatives sit at negative `N` (at rank ≤ 2, all of them have
  `N ≤ 0`).
- A contraction or symmetrisation sums the terms its operand actually
  represents, and is represented if any of them is. A partly-unrepresented sum
  is therefore not silently zero.
- `Permute<Image>(T)^{a0 a1 …} = T^{a_Image[0] a_Image[1] …}`; an image that is
  not a permutation is not an overload. `TensorSymmetry` requires its
  generators to be signed permutations.
- A letter outside a `MultiIndex` alphabet is a hard error (thrown in a
  constant expression), not a SFINAE-friendly failure, so every accessor that
  takes component letters as template arguments checks them first with
  `if constexpr`.

### 4.4 The spectral side

Spectral operations are written against concepts (`SpinCoefficients`) rather
than the owning `SpinExpansion`, so a tensor's component can be raised,
lowered, evaluated or interpolated as it stands. The derivative operators —
`ð`/`ð̄`, the contravariant derivative (the surface gradient of Dahlen & Tromp),
the intrinsic derivative and the bundle maps — and how they differ from one
another are covered in the reference. An expansion too short to carry a
gradient's output returns zero at the lowest degree that can hold it rather
than refusing.

### 4.5 Interpolation

`Interpolate(field, scheme)` returns a callable of `(θ, φ)` that is itself a
valid field constructor argument. The spectral scheme is exact for a
band-limited field; the local schemes (bilinear, bicubic) are built on a
padded copy of the grid: a wrap column at `φ = 2π`, three periodic ghost
columns at each end so the spline's end conditions are away from anything
evaluated, and two polar rows computed from the expansion. This is why making
a local interpolant costs a forward transform. A local interpolant takes the
samples it already has, not a second evaluation of the expression.

---

## 5. Numerics

### 5.1 The Wigner recursion and its ceiling

`d^l_{nm}(θ)` is computed by recursion in `l` from a seed of order
`(sin θ)^m`. The seed underflows, and the column recovers only once
`l sin θ ≳ m`, so the recursion is sound only for
`lMax < e·|ln min|`. `MaxSafeDegree<Real>()` enforces
`⌊e(|ln min| − |ln ε|)⌋` — 1827 in double, 30 747 in long double — slightly
below that so the seed never goes denormal. `Wigner`, `WignerMatrices` and the
grid all refuse a degree above it (the grid because a generating grid has no
table to do the refusing). Above the ceiling the table would be wrong by order
one with no sign of trouble. Extending the range (an exponent carried per
column) is known and not built; see `open-issues.md`.

The seed row's binomial is formed directly while it fits (`|n| ≤ 508` in
double) and from logarithms above, without `lgamma`, whose global state is a
data race under threads.

### 5.2 Single precision is for storage

A recursion's range and accuracy are those of its arithmetic, and a three-term
recurrence loses about `n²ε` along a row however good the algorithm. So
`float` Wigner values (stored, transform-major, generated and in the spectral
interpolant alike) and `float` 3-j rows are computed in double and rounded on
store
(`WignerRecursionReal`). `float` then shares double's ceiling and accuracy and
pays off where it should, in memory and bandwidth. The generated path narrows
from the same double values as the stored one, which is what keeps them
bit-identical.

### 5.3 3-j symbols

The Schulten–Gordon recursion builds each row inward from both ends — the
stable direction — and matches the halves. Three points matter:

- **Each half stops where its values first fall**, not where the recurrence
  coefficient first rises (the SLATEC rule). The two coincide when a row has a
  classically allowed region and not when it is near-stretched, which is the
  top of every coupling sum; there the coefficient rule lets one half run
  downhill and loses digits or overflows.
- **The halves are joined through `hypot`**, so their scale ratio is never
  squared and cannot overflow where the entries themselves do not.
- **The phase is tracked, not read back.** A near-stretched row can span 200
  orders of magnitude; rescaling can flush its last element to zero, and a row
  whose sign was taken from it comes back negated with every magnitude
  correct. Neither completeness nor the residual detects that; only
  cyclic-permutation invariance does, which is why the tests use it at high
  degree.

Every row is checked against the recurrence that defines it before it is
returned; the tolerance, `ResidualTolerance`, is 50 (n + 1) ε, set from the
worst residual measured (about 2–3 units). The completeness relation is not
used as a check, because the algorithm normalises by it.

### 5.4 Measured and rejected

Polar truncation of the Legendre sum was measured at about 1.15× (against an
expected 1.5–2×) and lowers GEMM efficiency by shrinking its inner dimension,
so it is not implemented.

---

## 6. Layered fields

A layered field is a stack of spherical fields on a `RadialGrid`, stored
radius-major so that each radius is one field of a contiguous batch. A radial
grid may carry an element partition with **interfaces** (a repeated radius
holding the two one-sided values).

The radial operators (finite differences, Lagrange, element and spline
derivatives, `Resample`) are **conveniences**: user codes with serious radial
discretisations do their own. They are therefore held to one rule — never a
wrong number in silence — and no more:

- Finite-difference stencils are built within each element, one-sided at its
  ends; an element narrower than the stencil is refused.
- Lagrange differentiation refuses any interface (pointing at
  `ElementDerivative`); both refuse a repeated radius on a grid that cannot say
  what it means.
- `Resample` answers the top node of a target element from below, so a model
  resampled onto its own mesh is the identity.
- Barycentric weights are formed on nodes rescaled to an interval of length
  four, so physically scaled radii (metres) do not overflow.
- An operator that exposes `Radial()` is checked against the stack's grid by
  identity; a bare callable cannot be. An argument and a result on different
  radial grids are refused.
- `Resample` and `SplineDerivative` check every piece's node count against
  what the scheme needs before any work starts, so the error names the
  element.
- `ApplyToLines` never hands an operator aliasing input and output spans; it
  supplies a scratch line when the caller's two buffers are one.

`RadialMajor` holds coefficients `[(l, m)][r]` so each radial line is
contiguous, and keeps the radial grid its lines run along. Two routes lead
there: transform then transpose (`RadialMajor(Expand(f))`), or transform
straight into the lines (`ExpandToLines`, which is
`Batch::Interleaved(nR, nR)` on the output). They are bit-identical. The
direct route saves a whole radius-major copy of the coefficients — gigabytes
at lMax = 512 with a few hundred radii — but is not reliably faster: its
scatter has stride `nR`, and at power-of-two `nR` under threads it loses to the
tiled transpose. There is therefore no automatic choice; it is offered as the
memory-saving route.

---

## 7. Policies and tuning

Policies (`Execution`, `Chunking`, `FFTWpp` planner flag, `WignerValues`,
`TransformKernel`, interpolation `Scheme`) are values. Those fixed by the table
are chosen at construction; the rest can be changed with `With()` or passed per
call.

`Tuning.hpp` measures the alternatives on the caller's own problem and returns
**values**, never a configured grid, so nothing is substituted behind the
caller's back. A challenger must beat the incumbent by 10 % to displace it,
because that is the measured noise floor; `conclusive` means exactly "the
challenger displaced the incumbent". Tuning is cheap enough to run at start-up
up to lMax of a few hundred, which is why there is no persistence layer.

---

## 8. Conventions in the code

- **Sizes and indices are signed** (`std::ptrdiff_t`); `std::size_t` appears
  only where a standard container is indexed or sized. Sizes here are
  subtracted, compared with orders that may be negative and multiplied into
  offsets, and unsigned arithmetic gets each of those wrong in silence.
- **Optional dependencies are absent, not disabled.** Without a BLAS there is
  no `TransformKernel::Matrix()`; without `Interpolation` there are no local
  schemes, `SplineDerivative` or `Resample`. Asking for one is a compile error
  at the call site, not a run-time throw. The macros are
  `GSHTRANS_HAVE_BLAS` and `GSHTRANS_HAVE_INTERPOLATION`.
- **`#pragma once`**, not include guards: guards are hand-chosen names that
  must be unique (there are two `BundleMaps.hpp`), and a header-only library
  reached through one root never meets the case the pragma gets wrong.
- **Prefer a concept to a `static constexpr bool` chain** whose later terms
  name members the earlier terms guarantee: `and` short-circuits evaluation,
  not well-formedness, and clang rejects what GCC accepts.
- **Doxygen attaches a `///<` at the top of a function body to the enclosing
  function**, so such comments silently document the wrong entity.

---

## 9. Testing and measurement

These are rules for working on the library. Each one is here because ignoring
it has produced a wrong conclusion.

**Testing**

- A new test should be seen to fail on the code it is meant to catch. A test
  that cannot fail is worse than none, because it is counted.
- Prefer property tests where the property is cheap to state ("every component
  obeys every generator") to example tests of one member.
- Before "fixing" an API oddity, look for the `static_assert` that pins it; the
  oddity may be a decision.
- Build tangential test fields as `grad f + r̂ × grad g`. A field written
  `a(θ,φ) θ̂` is generally not smooth at the poles, and its expansion does not
  converge.
- Under Schulten–Gordon the completeness relation tests nothing (it is the
  normalisation); use the recurrence residual and cyclic permutations.
- Run the suite at CI's width, `OMP_NUM_THREADS=4`, before pushing. Assert
  bit-equality only between two routes under the *same* policy; different
  thread counts sum in different orders. A threaded transform reproduces
  itself exactly between calls, so a difference that moves between calls is a
  race; one that moves between processes may be `FFTWpp::Measure` choosing a
  different plan (use `Estimate` to rule that out).
- Build the no-dependency configuration
  (`-DGSHTRANS_WITH_BLAS=OFF -DGSHTRANS_WITH_INTERPOLATION=OFF`) before
  pushing changes to `SphericalGrid.hpp`; a member added inside an
  `#ifdef GSHTRANS_HAVE_BLAS` block breaks only there.
- ThreadSanitizer is not used: `libgomp` carries no TSan annotations, so every
  barrier is reported as a race.

**Benchmarking**

- Benchmark a Release build. An empty `CMAKE_BUILD_TYPE` runs about ten times
  slower with internally consistent figures; the benchmark warns when built
  without `NDEBUG`.
- Interleave the candidates in one session. A laptop's clock scaling moves a
  repeated measurement by tens of per cent, and the first pass runs
  systematically slow.
- Treat ten per cent as the noise floor for sequential kernel timings: moving
  a hot loop elsewhere in the binary shifts it by that much, so two builds that
  differ anywhere cannot be told apart below it. Perturb the baseline to see
  how far it moves by itself before believing a small difference.
- One OpenMP runtime per process. Mixing `libomp` and `libgomp` (a clang build
  against the OpenMP OpenBLAS) gives correct results and meaningless timings.
- Give `cmake --build --parallel` a number. With none, the Makefile generator
  starts every translation unit at once, and these are heavy.
