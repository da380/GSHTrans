# Directions to look into

Five things raised as worth doing, expanded here into something that can be
planned against. Each keeps the original note at the top and then says what it
means concretely, where it overlaps what is already planned, and what would
have to be decided.

Nothing here is scheduled. `core-plan.md` and `field-algebra-plan.md` remain
the authorities on the numerical core and the field algebra respectively, and
`gshtrans-reference.tex` on what exists. Where a direction below contradicts a
decision already taken, it says so rather than quietly overriding it.

---

## 1. Tangential tensor fields

> *As raised:* tangential tensor fields. Should be easy to implement and so
> specialise operations.

### What it is

A tangential tensor has no radial slot: every index takes values in
`{-1, +1}` rather than `{-1, 0, +1}`. It is the `T^{ΩΩ}` piece of the
decomposition `T = r̂r̂ T_rr + r̂ T^{rΩ} + T^{Ωr} + T^{ΩΩ}` that D&T open
Appendix C with, and it is what most geophysical surface objects are: surface
strain and stress, the metric of the sphere, the second fundamental form, the
spin-2 fields of CMB and geodesy.

### What it costs and saves

| rank | components | tangential | real d.o.f. |
|---|---|---|---|
| 1 | 3 | 2 | 2 |
| 2 | 9 | 4 | 4 |
| 3 | 27 | 8 | 8 |
| 4 | 81 | 16 | 16 |

The saving grows as `(3/2)^p`, so a factor of five at rank 4.

### What it needs

`MultiIndex` hard-codes the alphabet `{-1, 0, +1}`. Generalising it to a
per-slot alphabet — a `Slots` policy with `AllSlots` and `TangentialSlots` —
is the whole of the change, because everything downstream is already generic
over *whatever indices exist*: the orbit enumeration walks the group it is
given, the symmetry policies permute slots without caring what they hold, and
the storage layout groups by upper index.

Three consequences fall out and are worth knowing before starting.

**The upper index is constrained.** `N = Σαᵢ` still, but now `N ≡ p (mod 2)`
and `|N| ≤ p`. A tangential rank-2 tensor has `N ∈ {-2, 0, +2}` only, and its
four components are `(--), (-+), (+-), (++)`.

**Reality gets simpler, except under symmetry.** The all-zero multi-index does
not exist for `p ≥ 1`, so negation has no fixed point and every orbit has size
two: `2^{p-1}` complex fields, `2^p` reals, and — for `NoSymmetry` — **no
pinned components and therefore no second buffer**, which is the one piece of
awkwardness phase 4 introduced. Under a permutation symmetry a self-paired
component can still appear: for a symmetric tangential rank-2 tensor, `(-+)`
is mapped to `(+-)` by negation and back by the symmetry, so it is pinned
real. The count then is one complex plus one real, three reals a point, which
is a real symmetric 2×2 matrix. That the machinery gets there unaided is the
check to write first.

**It is closed under the algebra but not under the gradient**, and this is the
important one. Products of tangential tensors are tangential; contraction of
two tangential slots against the metric stays tangential, since `g_{00}` needs
radial slots. But

```
(∇₁v)^{σ,0} = Ω v^0 − v^{0+σ} = −v^{σ}
```

for a tangential `v`, because the connection term reaches a component with a
radial slot. That term is the extrinsic curvature of the sphere — for a unit
sphere the normal part of `∇ₐv_b` is exactly `−v`, which is what this says —
so it is not an artefact to be designed away. **`SurfaceGradient` maps
tangential rank `p` to *general* rank `p+1`.** A type that claimed otherwise
would be wrong.

### What it does not express

Tracelessness. The traceless symmetric tangential rank-2 tensor — the spin-2
object — has two real degrees of freedom, but tracelessness is a linear
constraint and not an orbit of a group acting on indices, so the machinery of
`Orbits.h` cannot produce it. It would be a separate mechanism, and probably
not worth one: a caller can subtract the trace.

### Overlap

This reopens a decision. §16.3 of the field-algebra plan rejected a
"tangential tensor" type, on the grounds that it would complicate more than it
saved. That judgement was about one narrow question — whether to avoid storing
`SurfaceGradient`'s zero `σ = 0` block — and it is not an argument against the
type in general. The general case is stronger than the case I was answering,
and the storage table above is the reason.

---

## 2. Three-dimensional tensor fields

> *As raised:* 3D tensor fields. Built from a fixed angular grid and then a
> specified set of radii. These are just "dumb" storage objects in the library,
> but allow for things like transformations, and pointwise and tangential
> operations. Thought needs to be given here a bit as the aim would be to
> interact with, say, a finite-difference or finite-element discretisation in
> the radial direction that provides a means for doing radial differentiation.

### Overlap: most of this is already planned

`field-algebra-plan.md` §8, *Layered (3D) fields*, settles a good deal of it
and should be read first. In summary of what is already decided there:

- **2D is the primitive and 3D is a stack.** The angular field is
  first-class and never wrapped; a 3D field is one contiguous buffer viewed as
  `nR` angular slices, radius-major, exposed through a `Layered<…>` wrapper
  whose slice accessor returns an ordinary phase-1 node. There is no second
  expression system.
- **Bridges are explicit**: `Broadcast` lifts a 2D expression to every radius,
  `Slice` goes down, and there are no implicit conversions.
- **The wrapper is named neutrally** because time levels and ensembles want the
  same structure.
- **The radial mesh abstraction needs only nodes, quadrature weights and handle
  identity**; element connectivity stays in the application.
- **Parallelism is over slices**, on the thread-safety contract phase 1
  documented, with first touch matching the compute on NUMA machines.

### What the note adds, and what it forces

**The radial derivative is supplied, not owned.** The library should take an
operator rather than implement one: something that maps nodal values along `r`
to their derivative. For finite differences that is banded, for spectral
elements block-diagonal, and neither is the library's business. The interface
is the decision — a callable applied along the radial axis, or a sparse matrix
the library multiplies by — and the answer probably depends on whether the
caller wants to reuse a factorisation.

**Which representation it acts on.** The radial derivative commutes with the
angular transform, so it can be applied to samples or to coefficients. The
plan already notes that angular transforms want `[r][(l,m)]` and radial
operations want `[(l,m)][r]`, and that both are hot; that is exactly this
question, and the `Layout` policy plus an explicit repack is the answer it
reaches.

**It completes the gradient, and reinstates a factor.** With `∂⁰ = ∂/∂r`
available, `SurfaceGradient` becomes the angular half of the real `∇`. One
detail that must not be lost: D&T's `∂^±` carry `r^{-1}`, which the library
omits because an angular grid has no radius. In three dimensions it comes
back, and the full gradient is `[ê₀ ∂_r + r^{-1} ∇₁]`. A 3D gradient that
forgot the `r^{-1}` would be wrong by a factor that is invisible on the unit
sphere.

**Sizes are the design pressure.** A single complex scalar component at
`lMax = 256` and `nR = 100` is 210 MB; a real rank-2 tensor field at those
sizes is nine of them in reals. So laziness is worth real memory traffic here
in a way it is not in 2D, and the transform's batch axis and the radial axis
are the same axis — which is what step F was built for and what has not yet
had a consumer.

---

## 3. Specialisations for common objects

> *As raised:* Specialisations of the tensor classes for common objects.
> VectorFields, SymmetricSecondOrderTensorFields, things like that.

### What exists

Three aliases, in `TensorField.h`: `VectorField`, `Rank2TensorField`,
`ElasticTensorField`. They are thin and undocumented, and they were written to
make the tests readable rather than as a considered public surface.

### What a considered set would be

The objects that actually recur, with their tangential variants once §1 exists:

| name | rank | symmetry | typical use |
|---|---|---|---|
| `ScalarField` | 0 | — | already `SpinField<0>`; the alias is for uniformity |
| `VectorField` | 1 | none | displacement, flow |
| `SymmetricTensorField` | 2 | symmetric | strain, stress |
| `AntisymmetricTensorField` | 2 | antisymmetric | rotation, vorticity |
| `ElasticTensorField` | 4 | elastic | moduli |
| tangential forms of each | | | surface fields |

### The questions worth settling

**Aliases only, or extra API?** An alias costs nothing and reads better. Named
*accessors* — `v.Radial()`, `v.Plus()` — are a second vocabulary competing with
the multi-index one, and the library has been careful to have exactly one way
to name a component. My inclination is aliases plus a small number of named
*operations* that are meaningful only at a particular rank, where the name
carries real information: `Trace` (exists), `Deviatoric`, and — once §2 gives
`∂⁰` — `Divergence` and `Curl`, which are contractions of the gradient rather
than new operators.

**Defaults.** Every alias currently repeats `Reality` and `Grid`. Defaulting
`Reality` to `RealTensor` would match what applications mostly want, at the
cost of making the more surprising case the silent one.

This is small and mostly editorial, but it is the layer users see first, so it
is worth doing deliberately rather than by accretion. It should follow §1, so
that the tangential variants get named in the same pass.

---

## 4. Project structure

> *As raised:* Overall structure of the project, in terms of naming, file
> structure, CMake. Basically, this all probably needs to be updated and
> improved. Just surface level stuff, but of value.

Surface level, but there is one genuine bug in here and one genuine gap.

### CMake

- **`GaussQuad` is fetched over SSH**: `git@github.com:da380/GaussQuad.git`.
  Anyone without GitHub SSH keys cannot configure the project at all. This is
  a one-word fix to `https://` and it is the most consequential item on this
  list.
- **There are no install rules.** An `INTERFACE` library with no `install()`
  or `export()` cannot be consumed by `find_package`, so the only way to use
  GSHTrans is to vendor it or `FetchContent` it. Adding an install target and
  a package config file is standard and small.
- **Option names are leftovers**: `MY_PROJECT_BUILD_EXAMPLES`,
  `MY_PROJECT_BUILD_TESTS`, `MY_PROJECT_BUILD_BENCHMARKS`. They should carry
  the project's name.
- **`FetchContent` unconditionally**, where the convention is to try
  `find_package` first and fall back, so that a system Eigen or FFTW is used
  when present.
- **The version is `1.0` in `project()` and nowhere else** — no version header,
  nothing a consumer can test against.

### Layout and naming

- **`GSHTrans/src/`** holds headers, not sources. `include/` is the usual
  name, and `src/` will read oddly in an installed tree.
- **Extensionless umbrella headers** (`GSHTrans/All`, `Core`, `Field`,
  `Tensor`, `Expansion`) follow Eigen's convention. Fine, but it is a
  convention worth stating rather than inheriting.
- **`Concepts.h` has outgrown its name.** It now holds tag types, concepts,
  *and* four policy classes — `Execution`, `Batch`, `Chunking`,
  `WignerValues` — which are not concepts at all. Splitting the policies into
  their own header would make both halves easier to find.
- **Detail namespaces are inconsistent**: `WignerDetails`, `TensorDetails`,
  `EthDetails`, `ContravariantDetails`, `SymmetryDetails`, and a leftover
  `Internal`. One convention, applied everywhere.
- **`3j.h` is library code that nothing uses and nothing tests.** It is
  included in `GSHTrans/Core`, so it is public API, and it has no coverage at
  all. Either it earns a test or it should not be shipped. Note that
  `examples/wigner3j.hpp` and `wigner3j_tests.cpp` are a *different*,
  standalone implementation that does not use the library — two Wigner 3-j
  codes in one repository, neither exercised by the suite.

### Effort

All of it is a day at most, none of it blocks anything, and the SSH URL is
worth fixing on its own before anyone else tries to build this.

---

## 5. Dependencies

> *As raised:* Look into dependencies that are my own, getting plans formed
> that can be handed over for them to be updated. Both in general, and in light
> of needs linked to this project. GaussQuad and FFTWpp are the main ones. […]
> we should be open to replacing these bespoke codes with more standard ones if
> they are lightweight — prefer header only libraries beyond really key things
> like FFTW3 and maybe blas.

### The current set, and how deep each goes

| dependency | surface actually used | assessment |
|---|---|---|
| FFTW3 | via FFTWpp | external, required, keep |
| FFTWpp | ~12 names, incl. `Ranges::{Layout,Plan,View}` | load-bearing |
| GaussQuad | 6 entry points | trivially replaceable |
| NumericConcepts | 11 concepts | stable, low risk |
| Eigen | **none** | fetched and linked, never included |

**Eigen is the easy decision.** Nothing in the library includes it; the only
occurrence of the word is a comment about Eigen's operand-lifetime hazard. It
is fetched on every fresh configure for nothing. Either phase 3's deferred
"Eigen-backed pointwise linear algebra at ranks 2 and 4" gets built, or the
dependency should go.

### GaussQuad

The surface is six entry points: `LegendrePolynomial<Real>{}.GaussQuadrature(n)`,
then `Quadrature1D`'s `Transform`, `Points()`, `Weights()`, `X(i)`, `W(i)`.
That is small enough that replacing it is a contained job.

Measured, on this machine:

| nodes | build time | `Σw − 2` | `∫P²_{n-1}` relative error |
|---|---|---|---|
| 65 | 0.4 ms | `−1.1e−15` | `−3.0e−14` |
| 257 | 8.2 ms | `−5.7e−13` | `−2.8e−13` |
| 1025 | 145 ms | `3.5e−13` | `1.8e−13` |
| 2049 | 596 ms | `4.0e−12` | `2.0e−12` |

Two things to read from that. The cost is **O(n²)**, which is the signature of
Newton iteration with an O(n) polynomial evaluation per step; and the accuracy
**degrades with n**, from `10^{-15}` to `4×10^{-12}` by n = 2049, where a good
implementation holds near machine precision at any n.

Neither is a problem at the sizes used today — 8 ms and `10^{-13}` at
`lMax = 256` are invisible next to a 648 MB Wigner table — but both bite at
the sizes the library is being pointed at, and the O(n²) becomes the *dominant*
construction cost for a generating grid, which builds no table at all.

**The plan to hand over** would be: keep the `Quadrature1D` interface exactly,
replace the root-finding with Bogaert's method (explicit asymptotic expansions
for the nodes and weights, O(1) per node, accurate to near machine precision
for any n and used by most modern libraries), and add an accuracy test at
n = 10, 100, 1000, 10000 against `Σw = 2` and the orthogonality of the
highest-degree polynomial. Alternatively adopt a public-domain implementation
of the same method — it is a few hundred lines and header-only, which is
exactly the "lightweight standard replacement" criterion.

### FFTWpp

Harder to displace, and probably should not be: it wraps the advanced
interface, and it is the `Ranges::Layout` / `plan_many` support that makes the
batched transform's `(count, stride, dist)` descriptor cost a descriptor
rather than a repack. That is load-bearing and was confirmed by measurement,
not assumed.

Being a wrapper it should be performance-neutral, and the one thing worth
checking is whether it is: whether any path copies where it could plan. Two
things this project would ask of it if they are not there already:

- **Alignment control.** `core-plan.md` F3 was about FFTW's new-array execute
  being valid only for buffers in the same alignment class, which neither FFTW
  nor FFTWpp checks. Either a checked new-array execute or the ability to plan
  `FFTW_UNALIGNED` would let the transform skip the row copy it currently
  makes — measured invisible today, but it is the kind of thing that stops
  being invisible.
- **Wisdom import and export.** The grid used to pre-generate wisdom and set
  `WisdomOnly`, which step E removed because it could not anticipate every
  batched shape. A caller who wants persistent wisdom across runs currently
  has no route to it.

### NumericConcepts

Concepts only, and the reason it exists is sound: two projects agreeing on
what "real", "complex" and "precision" mean, rather than each defining a
subtly different set. `core-plan.md` records that the two *did* once differ,
over whether a const-qualified `std::complex` counts as complex. Low risk,
nothing to do.

---

## Suggested order

1. **The SSH URL in `CMakeLists.txt`** (§4). One word, and until it is fixed
   nobody else can build the project.
2. **Drop or use Eigen** (§5). One decision, and it removes a fetch.
3. **Tangential tensor fields** (§1). Self-contained, and the storage argument
   is strong at rank 4.
4. **Specialisations** (§3), immediately after, so the tangential forms are
   named in the same pass.
5. **The rest of the structural work** (§4), as a single deliberate tidy.
6. **The GaussQuad plan** (§5), written and handed over — it is someone else's
   work and does not block anything here.
7. **Three-dimensional fields** (§2), which is the largest and the one that
   completes the gradient, and which should start from `field-algebra-plan.md`
   §8 rather than from scratch.
