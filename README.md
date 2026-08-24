# GSHTrans

A C++23 header-only library for fields on the sphere: generalised spherical
harmonic transforms, and a lazy, index-checked algebra over spin-weighted
fields.

A spin-weighted field carries an **upper index** `N` — a spin weight — which
says how it transforms under a rotation of the local frame
`e_± ↦ e^{∓iψ} e_±`. The library tracks `N` in the type system, so the index
arithmetic of an expression is checked when it is written rather than assumed
when it is run.

## Status

The library is complete through the layers below.

* **The numerical core** — `GaussLegendreGrid`, `Wigner`, `Indexing`, `Views`:
  the batched transform, the plan cache, threading, and Wigner values either
  held as a table or generated on the fly, as a construction-time policy. Two
  Legendre kernels are carried permanently, `TransformKernel::Loop()` and,
  where a BLAS is present, `TransformKernel::Matrix()` — the second worth 3–6×
  where it is worth anything and storing half the table, with the first kept
  as its oracle. `Tuning.h` chooses between them by measuring the caller's own
  problem.
* **The field layer** — spin fields, tensor storage and algebra, the reality
  reduction, and the spectral side; the contravariant derivative, which is
  D&T's surface gradient; layered (three-dimensional) fields, the radial seam
  and `RadialMajor`; tangential tensors, the intrinsic derivative and the
  bundle maps; ready-made radial derivatives and resampling; the element
  partition on `RadialGrid`; and interpolation of a field as a callable of the
  two angles.
* **The Wigner 3-j symbols**, which touch no grid, no transform and no field:
  the Schulten–Gordon recursion, checked against the recurrence that defines
  it.

`docs/canonical-components.tex` is the authority on the mathematics and the
conventions, and the code defers to it by name. `docs/gshtrans-reference.tex`
describes the library that exists. `docs/lessons.md` is a short record of the
things worth not learning twice.

## The spin-field algebra

```cpp
using Real = double;
using Grid = GSHTrans::GaussLegendreGrid<Real, GSHTrans::All, GSHTrans::All>;

// Quadrature headroom for a product of two band-limited fields.
auto grid = Grid::ForBand(8, 2, 2.0);

auto u = SpinField<2, Grid>(grid, [](auto theta, auto phi) { ... });
auto v = SpinField<2, Grid>(grid, [](auto theta, auto phi) { ... });

// conj(u) carries upper index -2, so the product lands at zero and can be
// integrated. Nothing is evaluated until Integrate walks the expression.
const auto pairing = Integrate(conj(u) * v);
```

Every node — terminal, view or expression — satisfies the `SpinWeighted`
concept, which fixes `UpperIndex`, `Value`, the grid handle, `operator[](iTheta,
iPhi)` returning by value, and `EvaluateInto(span)`. Operators return lazy
nodes; evaluation happens on assignment, on `Materialise`, or on `Integrate`.

`examples/FieldExample.cpp` is the worked version of the above.

### The index rules

| expression | upper index | requires |
| :--- | :--- | :--- |
| `f + g`, `f - g` | `N` | `N_f == N_g` |
| `f * g` | `N_f + N_g` | — |
| `f / g` | `N_f` | `N_g == 0` |
| `-f` | `N` | — |
| `conj(f)` | `-N` | — |
| `abs(f)`, `abs2(f)` | `0` | — (result is `RealValued`) |
| `real(f)`, `imag(f)` | `0` | `N == 0` |
| `s * f`, `f * s`, `f / s` | `N` | — (a complex `s` promotes the value kind) |
| `s / f` | `0` | `N == 0` |
| `Map(f, F)` | `0` | `N == 0` |
| `Integrate(f)` | — | `N == 0` |

Two of these are worth stating as facts rather than as table rows, because the
superseded layer had them wrong. **Conjugation reverses the upper index**: `f`
at `N` conjugates to a field at `−N`. And **`real` and `imag` exist only at
`N = 0`**, because that is the only upper index at which they are covariant.

A third constraint lives on the concept rather than on any operator: a node may
be `RealValued` **only** at `N = 0`. Real-valuedness is not preserved by the
frame rotation, so it is not a property any component of any tensor can have at
nonzero spin weight. The rule is closed under every node in the algebra, so it
costs no expressiveness — it forbids only a real-valued *terminal* at `N ≠ 0`.

Unlawful combinations fail overload resolution at the call site, so they can be
written as negative tests: `static_assert(!requires { u + v; })`.

`pow`, `exp`, `log` and anything else of that kind are `Map(f, F)` at `N = 0`
rather than named nodes.

### Aliasing

Every node is *pointwise and index-preserving*: the value at `(iTheta, iPhi)`
reads its operands only at `(iTheta, iPhi)`. So `u = expr` is safe when `u`
appears in `expr`, and needs no temporary. Operations that are not pointwise —
gradients, raising and lowering, radial derivatives — deliberately live in the
spectral layer instead.

## Grids and transforms

`GaussLegendreGrid<Real, MRange, NRange>` is a **value-semantic handle** over a
shared immutable implementation: copying one is a pointer copy, and two grids
are the same grid when `Identity()` matches. It owns the quadrature, the Wigner
table, and a per-thread cache of FFTW plans and their work buffers.

* `GaussLegendreGrid(lMax, nMax, flag, chunking, values, kernel)` — a grid of
  that resolution.
* `Grid::ForBand(lBand, nMax, oversampling)` — a grid for fields of band
  `lBand` with quadrature headroom to degree `oversampling * lBand`. The 3/2
  dealiasing rule is `oversampling = 1.5`; `2.0` is exact for a single product.

`nTheta = lMax + 1`, and `nPhi` is the least fast FFT length `≥ 2·lMax + 1`, so
the orders `m = ±lMax` are resolved separately. Samples are stored
colatitude-major with longitude fastest: flat index `iTheta * nPhi + iPhi`.

The transform primitive is **a batch of `k` same-spin fields**, with the single
field as `k = 1`:

```cpp
grid.ForwardTransformation(lMax, n, in, inBatch, out, outBatch, policy);
grid.ForwardTransformation(lMax, n, in, out);   // the k = 1 wrapper
```

* A batch is a `Batch{count, stride, dist}` — element `j` of field `k` lives at
  `j * stride + k * dist`. `Batch::Contiguous`, `Batch::Interleaved` and
  `Batch::Strided` name the usual cases. Fields and coefficients take separate
  descriptors, since `dist` differs between them.
* A batch shares grid, degree **and upper index**: the Wigner block for that
  index is exactly what batching amortises. "Batch all my fields" is the
  natural expectation and the wrong one.
* Threading is an explicit per-call `Execution` policy, **defaulting to
  sequential**. The library creates threads only when asked, and exactly one
  level threads: an operation asked to run in parallel from inside an existing
  parallel region runs sequentially instead.
* Size, degree and upper-index checks throw `std::invalid_argument` in **all**
  build modes, release included.

Real-valued fields, with their reduced `m ≥ 0` coefficient storage, exist only
at `n = 0`; a real transform at `n ≠ 0` throws.

**Two Legendre kernels, chosen at construction and both kept.**
`TransformKernel::Loop()` is the default and is what the library has always
done. `TransformKernel::Matrix()` is the transform-major arrangement of the same
sum: every FFT first, then one `dgemm` per order against a table laid out
`[n][m][l][θ]`. It needs a BLAS — `GSHTRANS_WITH_BLAS`, which
is `AUTO` by default and found rather than fetched — and a build without one
does not offer it at all.

It is worth **5.8–6.0× forward and 2.9–3.2× inverse** batched at `lMax = 256`
on eight threads, and little unbatched there, where the loop kernel is already
at the memory roof and nothing is available. Products at different orders
write disjoint output, so it carries no accumulator and no reduction in either
direction, which is where the forward transform's threading used to lose. It
also stores **half** the Wigner values — 325 MB against 648 at `lMax = 256`,
2.6 GB against 5.2 at 512 — the negative orders coming from
`d^l_{nm}(π − θ) = (−1)^{l+n} d^l_{n,−m}(θ)`.

The pair is kept rather than one replacing the other, so that the loop kernel
is the matrix kernel's oracle — the same inputs through two independent
arrangements of the same sum — and so that which kernel suits a machine is a
question that machine can answer. `benchmarks/TransformBenchmark kernels` is
that question.

A BLAS with a thread pool of its own must be told to use one thread: these
products are skinny and threading them loses. A BLAS on the same OpenMP
runtime needs nothing, since every GEMM here is issued from inside a region.

## Interpolating a field

`Interpolate(field, scheme)` turns a field — or an expansion, or any lazy
expression — into a callable of `(θ, φ)`. It models `ScalarFunctionS2`, which
is what a field constructor takes, so remeshing onto another grid is one line:

```cpp
auto at    = Interpolate(field, Scheme::Bicubic());
auto value = at(theta, phi);
auto moved = SpinField<N, Grid>(otherGrid, Interpolate(field));
```

| scheme | exact? | per point |
| :--- | :--- | :--- |
| `Scheme::Spectral()` | yes, for a band-limited field | `O(lMax²)` |
| `Scheme::Bilinear()` | no, second order | `O(1)` |
| `Scheme::Bicubic()` | no, fourth order | `O(1)` |

The local two come from `Interpolation` and are absent — not refused — in a
build without it. They are about **4400×** cheaper per point than the spectral
sum at `lMax = 128`, and they want an oversampled grid: at the band limit
bicubic carries 13% error, at 4× oversampling `7e-4`, at 8× `4e-5`. `ForBand`
is how you ask for the room.

Two things about the sphere that a rectilinear scheme does not know are fixed
by handing it a **padded** grid rather than the field's own. The longitudes
stop one step short of `2π`, so a wrap column is added — exactly, since
`φ = 2π` is `φ = 0`. And neither pole is a grid point, so two polar rows are
added, computed from the expansion rather than guessed; that is why building a
local interpolant costs a forward transform, and why its cheapness is per
evaluation rather than per interpolant.

At a pole only one order survives, so the value there is `c·e^{±iNφ}` — a
*row*, not a constant. That is not a defect: a spin-weighted field at a
coordinate pole is genuinely not single-valued, because `e_±` depends on the
azimuth of approach. Measured, the polar and wrap cells come out two to five
times **better** than the interior, which is the padding doing its job.

A colatitude outside `[0, π]` throws; a longitude is reduced modulo `2π`,
which is exact. The two axes are not alike and their boundaries fail
differently.

## Choosing a policy by measuring it

Some of the library's choices cannot be settled by reasoning and vary by
machine. `Tuning.h` times the alternatives on the caller's own problem and
hands back **values**, never a configured grid — so nothing is substituted
behind your back, which matters most for the thing most likely to substitute
silently.

```cpp
auto chunk  = TuneChunking(grid, lMax, n, count, policy);
auto kernel = TuneKernel<Grid>(lMax, nMax, n, count, policy);
auto tuned  = grid.With(chunk.chunking);   // a pointer copy, one table
```

`grid.With(...)` changes the chunking policy or the planner flag without
rebuilding the table — neither decides it — and the result shares `Identity()`
with its parent, so fields are interchangeable between them.

Measured on the development laptop: the kernel choice is worth **1.9–5.8×**
and picks the matrix kernel in every configuration tried; the chunk is worth
at most **1.28×** and nothing at all in thirteen of eighteen, the conservative
default being good. Both are cheap enough to run at start-up, which is why
there is no persistence layer; that would change at around `lMax` 512, where
building the tables to measure them stops being cheap. A candidate must beat the incumbent by 10% to displace it,
since that is the measured noise floor and picking the winner of a 7%
difference is picking noise.

## The Wigner functions

`Wigner`'s stored value at upper index `N`, degree `l`, order `m` is

> `sqrt((2l+1)/(4π)) · d^l_{Nm}(θ)`,

where `d^l_{Nm} = P^N_{lm}(cos θ)` is the generalised Legendre function of
**Dahlen & Tromp (1998) eq. (C.115)** — **upper index first**. This is pinned
by `tests/CheckWignerConvention.h` against the `l = 1` table, and corroborated
by `tests/CheckLegendre.h`, which fixes the `N = 0` row against
`std::sph_legendre`.

Values come from stable recurrence relations, computed in parallel over
`(n, θ)`. Orthonormalisation is the only normalisation offered: there is no
`Normalisation` template axis.

## Wigner 3-j symbols

`3j.h` gives the coupling coefficients, which is what Gaunt integrals and
mode coupling need. The table is the primitive — a single symbol costs a whole
table, so ask for the table:

```cpp
auto table  = Wigner3jMatrix<double>(l1, l2, l3);   // the whole (m1, m3) plane
auto value  = table(m1, m3);                        // m2 = -(m1 + m3)
auto stack  = Wigner3jStack<double>(l1, l3);        // one table per middle degree
```

Values come from the **Schulten–Gordon** recursion: each row is built inward
from both ends of its range — the stable direction, since recursing outward
follows the decaying solution and loses the answer exponentially — matched
where the halves overlap, and normalised from the unitary property. No
closed-form seed, no factorials, and the turning point is found by watching
the recurrence coefficient rather than locating it analytically.

That matters because a one-directional recursion loses these values near
*stretched* triangles, where one degree approaches the sum of the other two —
which is the top of every coupling sum rather than an exotic corner. Measured
against an independent route (cyclic-permutation invariance, which runs the
recursion along different lines) it agrees to `1e-16` at `(200,200,200)` and
`1e-15` at `(1000,1000,1999)`.

Every row is checked against the recurrence that defines it before it is
handed back, so a table that has gone wrong is refused rather than returned.
That check replaced the completeness relation, which stopped testing anything
once the algorithm began normalising by it.

`CouplingElement(m, mp)` and `FillCouplingMatrix` give the same symbols in the
layout normal-mode codes expect — the first order negated with an alternating
phase, which is what Woodhouse's original routine returned. It is a convention, not a second
calculation, and the conversion is its own inverse.

6-j is not implemented. If it is ever wanted, the same paper covers it.

## Building

```
cmake -S . -B build
cmake --build build -j
cd build && ctest
```

Requires a C++23 compiler (GCC 13+), CMake 3.20+, and a local FFTW. OpenMP is
used for the parallel paths. `GaussQuad`, `FFTWpp`, `NumericConcepts` and
`Interpolation` are looked for on the system and fetched by `FetchContent`
only if they are not there. All four are header-only, and none of them brings
Eigen: GaussQuad used to, and no longer does.

| option | default | effect |
| :--- | :--- | :--- |
| `GSHTRANS_BUILD_EXAMPLES` | `ON` | build and register `examples/` |
| `GSHTRANS_BUILD_TESTS` | `ON` | build `tests/` |
| `GSHTRANS_BUILD_BENCHMARKS` | `ON` | build `benchmarks/TransformBenchmark` |
| `GSHTRANS_INSTALL` | `ON` when top level | generate the install and export rules |
| `GSHTRANS_WITH_INTERPOLATION` | `ON` | radial resampling, spline derivatives, and the local interpolation schemes |
| `GSHTRANS_WITH_BLAS` | `AUTO` | the matrix transform kernel. `ON` fails the configure without a BLAS; `OFF` never looks |

**Both optional dependencies are absent rather than disabled.** Without
`Interpolation` there is no `Scheme::Bicubic()` to call and no
`RadialSplineDerivative.h` to include; without a BLAS there is no
`TransformKernel::Matrix()`. Asking for one is a compile error at the call
site rather than a throw at run time, and CI builds with both off so that the
claim is run rather than asserted — 344 tests there against 380.

### Using it from another project

```cmake
find_package(GSHTrans REQUIRED)
target_link_libraries(your_target PRIVATE GSHTrans::GSHTrans)
```

`add_subdirectory` and `FetchContent` work too, and give the same target name.
`GSHTrans/src/Version.h` defines `GSHTRANS_VERSION` for feature tests against
a particular release.

`tests/package` is a standalone project that consumes the installed package;
it is what CI uses to check that the export rules still produce something
usable, which the main build cannot see.

### Testing

The suite runs in Debug, in Release, and under AddressSanitizer and
UndefinedBehaviorSanitizer -- `scripts/test_sanitized.sh address` does the
latter. The examples are registered as tests, so they are run rather than
merely compiled. Every public header is additionally compiled on its own, and
twice, so that one which stops standing alone fails the build rather than
waiting for the first consumer that does not already pull the missing include.

The sanitiser jobs run `Debug` with **leak detection on**, which is not what a
casual local run does; `scripts/test_sanitized.sh address` is the exact
configuration and is worth using rather than a hand-rolled one. The clang leg
is advisory and allowed to fail, since no clang OpenMP runtime is installed on
the development machine — but it earns its place: it is what found a
`static constexpr bool` constraint whose later terms named members a
non-spin-weighted operand does not have, which GCC accepted and clang did not.

ThreadSanitizer is deliberately not offered. The parallelism here is OpenMP,
GCC's `libgomp` carries no TSan annotations, and every barrier and reduction
is therefore reported as a race; that is a limitation of the tooling, not a
finding.

`TransformBenchmark` is not a test and `ctest` does not run it. It takes
section names so that an A/B costs one section rather than the whole run:

```
stream grid transforms threading batching generated
kernels kernels-loop kernels-matrix interpolation tuning server huge
```

`benchmarks/run-server-benchmark.sh` drives it on a target machine. **Build it
Release.** `cmake -S . -B build` leaves `CMAKE_BUILD_TYPE` empty, and the
harness there runs about ten times slow with every figure internally
consistent — it once reported speedups of forty. It now warns when built
without `NDEBUG`.

## Layout

```
GSHTrans/Core          umbrella: grid, Wigner, indexing, policies, tuning, 3j
GSHTrans/Field         umbrella: the spin-field algebra
GSHTrans/Tensor        umbrella: tensor fields and their algebra
GSHTrans/Expansion     umbrella: the spectral side
GSHTrans/Layered       umbrella: three-dimensional fields
GSHTrans/All           all of them
GSHTrans/src/          the headers themselves
docs/                  the theory note, the reference, the lessons
tests/  examples/  benchmarks/  scripts/
```

`examples/` is a numbered series meant to be read in order, each introducing
one thing and assuming the ones before it; `examples/README.md` lists them.
They are registered as tests, so they run rather than merely compile.

## License

BSD 3-Clause. Copyright (c) 2025, David Al-Attar.
