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

The rebuild is complete through the layers below, against two planning
documents in `docs/`. Both record their decisions and the measurements behind
them, including the ones that turned out to be wrong.

* **`docs/core-plan.md`** — the numerical core: `GaussLegendreGrid`, `Wigner`,
  `Indexing`, `Views`. Steps A–H are landed, with the batched transform, the
  plan cache, threading, and Wigner values generated on the fly as a
  construction-time policy. What remains is the transform-major restructure
  (step G with step F tier 2) and polar truncation, both gated on measurements
  from a target machine rather than on design.
* **`docs/field-algebra-plan.md`** — the field layer, and everything built on
  it:

  | | |
  |---|---|
  | phases 1–5 | spin fields, tensor storage, tensor algebra, reality reduction, the spectral side |
  | §16 | the contravariant derivative — D&T's surface gradient |
  | §17 | layered (three-dimensional) fields, the radial seam, `RadialMajor` |
  | §18 | tangential tensors, the intrinsic derivative, the bundle maps |
  | §19 | ready-made radial derivatives, and resampling |
  | §20 | the element partition on `RadialGrid` |
  | §21 | what the `Interpolation` update changed |

`docs/canonical-components.tex` is the authority on the mathematics and the
conventions; both plans defer to it. `docs/gshtrans-reference.tex` describes
the library that exists.

Read `core-plan.md` §8 before touching core code, and
`field-algebra-plan.md` §18.2 and §19.5 before touching the tensor bundles or
the radial seam — those are where the decisions live.

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

* `GaussLegendreGrid(lMax, nMax, flag, chunking)` — a grid of that resolution.
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

## Building

```
cmake -S . -B build
cmake --build build -j
cd build && ctest
```

Requires a C++23 compiler (GCC 13+), CMake 3.20+, and a local FFTW. OpenMP is
used for the parallel paths. `GaussQuad`, `FFTWpp` and `NumericConcepts` are
looked for on the system and fetched by `FetchContent` only if they are not
there. All three are header-only, and none of them brings Eigen: GaussQuad
used to, and no longer does.

| option | default | effect |
| :--- | :--- | :--- |
| `GSHTRANS_BUILD_EXAMPLES` | `ON` | build and register `examples/` |
| `GSHTRANS_BUILD_TESTS` | `ON` | build `tests/` |
| `GSHTRANS_BUILD_BENCHMARKS` | `ON` | build `benchmarks/TransformBenchmark` |
| `GSHTRANS_INSTALL` | `ON` when top level | generate the install and export rules |

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

ThreadSanitizer is deliberately not offered. The parallelism here is OpenMP,
GCC's `libgomp` carries no TSan annotations, and every barrier and reduction
is therefore reported as a race; that is a limitation of the tooling, not a
finding.

`TransformBenchmark` is not a test and `ctest` does not run it. It takes
section names (`stream grid transforms threading batching server huge`) so that
an A/B costs one section rather than the whole run;
`benchmarks/run-server-benchmark.sh` drives it on a target machine.

## Layout

```
GSHTrans/Core          umbrella: grid, Wigner, indexing, policies, 3j
GSHTrans/Field         umbrella: the spin-field algebra
GSHTrans/Tensor        umbrella: tensor fields and their algebra
GSHTrans/Expansion     umbrella: the spectral side
GSHTrans/Layered       umbrella: three-dimensional fields
GSHTrans/All           all of them
GSHTrans/src/          the headers themselves
docs/                  the plans, and the theory note
tests/  examples/  benchmarks/  scripts/
```

## License

BSD 3-Clause. Copyright (c) 2025, David Al-Attar.
