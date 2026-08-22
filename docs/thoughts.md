# Directions to look into

Things raised as worth doing, expanded here into something that can be planned
against. Each keeps the original note at the top and then says what it means
concretely, where it overlaps what is already planned, and what would have to
be decided.

Nothing here is scheduled. `core-plan.md` and `field-algebra-plan.md` remain
the authorities on the numerical core and the field algebra respectively, and
`gshtrans-reference.tex` on what exists. Where a direction below contradicts a
decision already taken, it says so rather than quietly overriding it.

**Status, 2026-08-22.** Sections 1 to 6 were the original set. Sections 2, 4
and 5 are largely done and say so at the top; 1, 3 and 6 are open. Sections 7
to 10 were raised later and are assessed here for the first time.

| | subject | status |
|---|---|---|
| 1 | tangential tensor fields | open, and the next candidate |
| 2 | three-dimensional tensor fields | **done** -- `field-algebra-plan.md` §17 |
| 3 | specialisations for common objects | open, follows §1 |
| 4 | project structure | **done**, except the `src/` rename, declined |
| 5 | dependencies | **done** upstream; GaussQuad and FFTWpp both refactored |
| 6 | Wigner 3-j symbols | open, and independent of everything else |
| 7 | the `Interpolation` library | assessed here |
| 8 | how three-dimensional the 3-D fields are | assessed here |
| 9 | interpolating a field, as a callable | assessed here |
| 10 | a wisdom mechanism for the computational options | assessed here |

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

### Two derivatives, not one

The thing that makes the type worth having is not the storage. It is that
**tangential tensors carry their own covariant derivative, and it is closed.**

There are two operators here and they must be kept apart.

- **The ambient surface gradient `∇₁`**, which is what D&T write and what
  `SurfaceGradient` computes: the angular part of the three-dimensional
  gradient. It differentiates the tensor *and its basis*, and the basis leaves
  the tangent plane.
- **The intrinsic covariant derivative `D`**, the Levi-Civita connection of the
  induced metric on the sphere. This is what surface differential geometry
  means by differentiation, and it is closed on tangential tensors by
  construction.

They are related by the Gauss formula: `∇₁` equals `D` plus a term in the
second fundamental form, which for the unit sphere is algebraic — no
derivatives at all. Concretely, for a tangential `T`, splitting `∇₁T` by
whether the inherited slot is radial:

```
(∇₁T)^{σ α₁…α_q}       = Ω^{∓N} T^{α₁…α_q}          all αᵢ tangential
(∇₁T)^{σ α₁…0ⱼ…α_q}    = −T^{α₁…σⱼ…α_q}             one slot radial
```

The first line says the tangential part of `∇₁T` has **no connection terms at
all**: every shift `αᵢ + σ` either leaves `{-1,0,1}` or lands on a radial slot,
which a tangential tensor does not have. So it is pure `Ω` multiplication —
which is `ð` up to `√2`. The second line is the extrinsic curvature: for a
unit sphere the normal part of `∇ₐv_b` is exactly `−v`, and that is what it
says.

Both were checked against the existing `SurfaceGradient`, applied to a rank-1
tensor with its radial component set to zero. Both hold **exactly**, to
`0.000e+00`.

So:

> **On tangential tensors the intrinsic covariant derivative is `ð` applied
> component by component, up to `√2`, and it is closed.**

That is worth stating plainly because it reconciles the two formalisms this
library carries. `ð` is not a poor relation of the contravariant derivative;
it is the *intrinsic* derivative, exact on exactly the objects — tangential,
spin-weighted — that it was invented for. The contravariant derivative is the
ambient one, and it is the right tool for general tensors with radial slots.
Neither subsumes the other, and the difference between them is one algebraic
term.

### What it needs

`MultiIndex` hard-codes the alphabet `{-1, 0, +1}`. Generalising it to a
per-slot alphabet — a `Slots` policy with `AllSlots` and `TangentialSlots` —
is the whole of the change, because everything downstream is already generic
over *whatever indices exist*: the orbit enumeration walks the group it is
given, the symmetry policies permute slots without caring what they hold, and
the storage layout groups by upper index.

On top of that, three small operators: the intrinsic derivative (which is
`ð` and therefore nearly free), the projection from a general tensor to its
tangential part, and the injection back. `SurfaceGradient` of a tangential
tensor then has an honest signature — it returns a general tensor — and a
caller who wants the intrinsic derivative asks for that instead of projecting
by hand.

Three consequences of the alphabet change fall out and are worth knowing
before starting.

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

**The ambient gradient still leaves the type**, as above. That is not a defect
to design around; it is the curvature, and a type that hid it would be lying.

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
type in general. The general case is much stronger than the case I was
answering: the storage table above, and above all the closed intrinsic
derivative, which the narrow question never raised.

---

## 2. Three-dimensional tensor fields

**Done.** Built as `field-algebra-plan.md` §17, in four steps recorded in
§§17.5-17.7: `RadialGrid`, `LayeredSpinField`/`LayeredSpinExpansion`, the
radial seam, layered tensors, `Gradient`, and `RadialMajor`. Everything this
section anticipated held, including the `r^{-1}` warning below. §8 is the
follow-on question -- *how much* the library knows about how the layers are
linked, which this section did not ask.

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

**Done, 2026-08-22**, except one item declined. The SSH URL, the option names,
`FetchContent` with `find_package` first, install and export rules, a version
header, CI, the policies split out of `Concepts.h` and the one stray detail
namespace are all in. What is *not* done is `GSHTrans/src/` -> `include/`:
GaussQuad and FFTWpp both kept `<Name>/src/` through their own refactors, so
the rename would be forty files of churn against the grain of the sibling
projects rather than towards a convention they share. The `3j.h` coverage gap
is §6's, not this section's.

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
- **`3j.h` has no test coverage at all**, although it is included in
  `GSHTrans/Core` and is therefore public API. It is not stray code — it is a
  documented facility and intended to stay — but it is untested, and §6 below
  says what testing it turns up. `examples/wigner3j.hpp` and
  `wigner3j_tests.cpp` are a *second*, standalone implementation that does not
  use the library, so there are two Wigner 3-j codes in the repository and
  neither is exercised by the suite.

### Effort

All of it is a day at most, none of it blocks anything, and the SSH URL is
worth fixing on its own before anyone else tries to build this.

---

## 5. Dependencies

**Done, and done upstream.** GaussQuad and FFTWpp were both refactored on the
strength of the plans this section produced; `GaussQuad-for-GSHTrans.md` and
`FFTWpp-for-GSHTrans.md` are the hand-over notes. Eigen is gone from the
dependency tree entirely -- GaussQuad replaced Golub-Welsch with an implicit-QL
eigensolver and the two `llt().solve()` calls with a Thomas solve -- and
`find_package(GSHTrans)` now works, which is what the whole exercise was for.
The one prediction below that did not survive is Bogaert: GaussQuad took
Glaser-Liu-Rokhlin as an option and kept Golub-Welsch as the default, because
GLR's weights drift as `O(n eps)` where the default's stay flat. §7 is the
same question asked about `Interpolation`, and the answer is now shorter
because these two have been through it.

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
| Eigen | none directly | **required transitively, by GaussQuad** |

**Eigen is not ours, and cannot simply be dropped.** No GSHTrans header
includes it — the only occurrence of the word in the library is a comment about
Eigen's operand-lifetime hazard — but *GaussQuad* includes
`Eigen/Cholesky`, `Eigen/Core` and `Eigen/Eigenvalues`, and fetches Eigen
itself. Removing `Eigen3::Eigen` from this project's link line is correct and
harmless, since the dependency belongs to GaussQuad and propagates from there;
removing Eigen from the build is not possible while GaussQuad needs it.

Why it needs it is the interesting part, and it changes the GaussQuad plan
below: `OrthogonalPolynomial.h` builds the Jacobi matrix and calls
`Eigen::SelfAdjointEigenSolver`. That is **Golub–Welsch**, and it explains
both measurements — the `O(n²)` cost and the accuracy that drifts with `n`,
since the nodes are eigenvalues of a matrix of size `n`.

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
replace Golub–Welsch with Bogaert's method (explicit asymptotic expansions for
the nodes and weights, O(1) per node, accurate to near machine precision for
any n, and used by most modern libraries), and add an accuracy test at
n = 10, 100, 1000, 10000 against `Σw = 2` and the orthogonality of the
highest-degree polynomial. Alternatively adopt a public-domain implementation
of the same method — a few hundred lines and header-only, which is exactly the
"lightweight standard replacement" criterion.

The case is stronger than "it could be better", because it is three things at
once: `O(n²)` becomes `O(n)`, `4×10⁻¹²` becomes `10⁻¹⁶`, and **the Eigen
dependency disappears** — Bogaert's method uses no linear algebra at all. That
would leave the whole dependency set header-only and small.

**And it would unblock installation.** Neither FFTWpp nor GaussQuad carries any
`install()` or `export()` rules, so neither can be found by `find_package`, so
GSHTrans cannot export a CMake package that references them: an `INTERFACE`
target can only be exported if everything it links is exported or imported.
GSHTrans therefore installs its headers and no package config, and consumption
is by `add_subdirectory` or `FetchContent`. Adding install and export rules to
those two is a small change in each and is the whole of what blocks a proper
`find_package(GSHTrans)`.

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

## 6. Wigner 3-j symbols

> *As raised:* The 3J stuff is intended to be built in at some point. The codes
> there were me messing about with a C++ translation of one of Woodhouse's old
> F77 routines. None will be ideal. From testing, the algorithm is pretty good
> for most practical cases, but not without edge issues, and it could be
> bettered using alternative approaches. But having some 3j functionality would
> be nice. 6j etc I don't personally care about but, if easy, why not.

### What is there

`GSHTrans/src/3j.h` is more than scaffolding: a documented C++20 port of
Woodhouse's `wig2.f`, offering `Wigner3jMatrix` (the full table over the
`(m₁, m₃)` plane at fixed degrees), `Wigner3jStack` (over `l₂` at fixed
`l₁, l₃`), a single-symbol convenience, and in-place kernels that fill
caller storage. Its stability argument is stated in the file: the corner value
is closed-form, two-term recursions run the edges, a three-term recursion
sweeps diagonals of constant `m₂` *from the classically forbidden corner
inward* — the direction in which the symbols grow — and the remaining half of
the plane comes from the reflection symmetry rather than from continuing into
an unstable regime. No factorials are formed.

### Where it stands up, and where it does not

I measured it, against exact values from `sympy` for small degrees and against
the orthogonality identity `Σ_{m₁m₃} (3j)² = 1` — which needs no reference at
all — for large ones.

**It is absolutely accurate over a wide range.** For triangles that are not
close to stretched, the error is at the noise floor: `10⁻¹⁶` to `10⁻¹⁴`
absolute, and the orthogonality sum holds to `10⁻¹⁵` for `(l,l,l)` up to
`l = 128`.

**The relative error on the smallest symbols is larger, and benignly so.** The
worst relative errors sit on the exponentially small symbols in the forbidden
corner — `2×10⁻⁸` relative on a symbol of size `3.5×10⁻⁶`, against a largest
symbol of order `0.1`. That is the same `10⁻¹⁴` absolute error seen from a
different angle. For anything that sums 3-j symbols, which is what coupling
coefficients do, it does not matter.

**It fails outright for stretched triangles**, and this is the real finding.
Taking `l₃ = l₁ + l₂`, the orthogonality sum departs from 1 by:

| `l` | `(l, l, 2l)`, double | long double | `(l, l, 3l/2)`, double | `(l, l, l)`, double |
|---|---|---|---|---|
| 20 | `4×10⁻¹⁴` | `10⁻¹⁸` | `2×10⁻¹⁶` | `10⁻¹⁶` |
| 30 | `4×10⁻³` | `8×10⁻⁹` | — | `4×10⁻¹⁵` |
| 40 | `2×10⁹` | `8×10²` | `8×10⁻¹⁵` | `6×10⁻¹⁵` |
| 60 | `2×10³²` | `5×10²⁵` | `3×10⁻⁴` | `2×10⁻¹⁵` |
| 100 | — | — | `5×10¹⁸` | `10⁻¹⁴` |
| 128 | `2×10¹¹²` | — | — | `2×10⁻¹⁰` |

The error grows **exponentially in the degree**, at roughly a decade per
degree, and long double buys about ten degrees before failing the same way.
That is the signature of a recursion being run in its unstable direction, not
of dynamic range: extra precision delays it and does not cure it. So the
stability argument in the file's header holds for fat triangles and breaks
down as the triangle approaches stretched, where the forbidden region is large
and the sweep must cross it.

**Why this is the regime that matters.** Coupling coefficients and Gaunt
integrals need `l₃` running all the way to `l₁ + l₂` — the stretched
configuration is not an exotic corner but the top of every coupling sum. For
band-limited fields at `lMax = 256` that means degrees to 512. The current
implementation is trustworthy to `l ≈ 30` there.

### What to do about it

Three routes, and they are not equally good.

**Fix the direction.** This is Schulten & Gordon's scheme: recurse inward from
*both* classically forbidden ends, match in the allowed region, and fix the
normalisation from the orthogonality sum rather than from a closed-form seed.
It keeps the existing interface and the no-factorials property, and it is what
the standard library routines do. Luscombe & Luban refined it.

**Compute exactly.** Prime-factorisation methods (Johansson & Forssén's
`wigxjpf` is the reference implementation) give full double precision at
arbitrary degree by doing the combinatorics in exact integer arithmetic, and
they cover 3-j, 6-j and 9-j together. The cost is a real dependency: C rather
than C++, a few thousand lines, and not header-only — which cuts against the
preference stated in §5, though the licence is permissive.

**Or hybridise, because the two classical methods fail in complementary
regimes.** Racah's closed form is a single alternating sum whose length is
roughly `min(l₁+l₂−l₃, …)`. At `l₃ = l₁ + l₂` that sum has exactly **one
term** and is therefore exact; near-stretched it has a handful. It is the fat
triangles, where the sum is long and alternating, that destroy it in floating
point — and those are precisely where the present recursion is at its best.
A dispatch on how close the triangle is to stretched would cover the whole
space with two simple methods and no new dependency. This is worth an
afternoon's investigation before committing to either of the others.

### Testing, whichever route

The orthogonality identity is the thing to build on: it needs no reference
implementation, it is cheap, and it caught this in one line. A test that
sweeps the triangle space — fat, intermediate and stretched, at several
degrees — plus exact comparison against rational arithmetic for small degrees
and the known closed forms (`(j j 0; m −m 0)`, and `l₃ = l₁ + l₂`) would have
made the boundary above visible from the start.

### 6-j

The answer depends entirely on the route. If the exact method is adopted, 6-j
and 9-j come with it and the question does not arise. If the recursion is
fixed by hand, 6-j is a separate implementation of the same Schulten–Gordon
idea — not hard, but not free either. That is an argument for deciding the
3-j route with 6-j in view, even though nothing needs 6-j today.

---

## 7. The `Interpolation` library

> *As raised:* another repo revisited, `da380/Interpolation`. Of potential use.
> Certainly they will link up downstream.

### What it is, as it stands

Header-only, C++20, `Interpolation/` with the same `<Name>/src/` layout as
FFTWpp and GaussQuad and the same extensionless umbrella headers. Six
facilities: `Linear`, `CubicSpline`, `Akima`, `Lagrange`,
`LagrangePolynomial`, `Polynomial1D`. Abscissae strictly increasing;
ordinates real or complex where noted.

**It is one-dimensional.** That is the single most important fact for how it
joins to this library, and it is a good thing rather than a limitation: the
radial axis is exactly where GSHTrans has deliberately stopped, so a
one-dimensional interpolator is the missing half rather than a competitor to
anything here. The angular half is not an interpolation problem of this kind
at all — see §9.

### Three things to settle before depending on it

These are not objections. They are the same three things GaussQuad and FFTWpp
had to fix, and both took about a fortnight, so the shape of the work is
known.

**It has no install or export rules.** `README.md` shows consumption by
`FetchContent` and a bare `Interpolation` target name with no namespaced
alias. That is precisely what blocked `find_package(GSHTrans)` until last
week: an INTERFACE target can only be exported if everything it links is
exported or imported. Taking a dependency on `Interpolation` as it stands
would *undo* §5's result. This is a small change and it is the one that has
to come first.

**Cubic splines pull in Eigen.** The dependency tree has just been cleared of
Eigen entirely, and this would bring it straight back — for a tridiagonal
solve. GaussQuad hit the identical case and its two `llt().solve()` calls
became a fifteen-line Thomas solve; a natural or clamped cubic spline is a
symmetric tridiagonal system and wants exactly the same treatment. So the ask
is one that has already been answered once in this family of libraries, and
the answer can be lifted.

**Every interpolator stores iterators rather than copying.** The README says
so: "Keep the sample containers alive and do not reallocate them while an
interpolator is in use." That is a lifetime contract of the kind this library
is careful about — `field-algebra-plan.md` §3.8 is entirely about operand
lifetime and aliasing, and §3.3 settled operand storage by value category
precisely so that a node cannot outlive what it reads. An interpolator handed
out by a GSHTrans field (§9) would inherit the hazard and hand it to a user
who never chose it. Either the interpolators gain an owning mode, or anything
built on them owns the samples itself and the borrowed form stays internal.

### Where it would actually be used

Four places, in decreasing order of how much they want it.

1. **Evaluating a layered field at an arbitrary radius.** The obvious one, and
   the one that needs nothing but `Linear` or `CubicSpline` along a gathered
   radial line. It composes with the seam: `ApplyRadially` already gathers a
   contiguous line and hands it to a callable, and this is the same gather
   with a callable that returns a value rather than a line.
2. **Remeshing between radial grids.** Model on one set of radii, computation
   on another. Common, tedious, and entirely one-dimensional.
3. **The spectral-element bases.** `LagrangePolynomial` gives the cardinal
   functions on the GLL nodes GaussQuad now places exactly at `±1`, and
   `Polynomial1D` differentiates them. That is the differentiation matrix a
   caller currently has to build by hand in order to pass a radial derivative
   to `Gradient`, and §8 argues it should be offered.
4. **`Polynomial1D` for the radial power laws** used throughout the layered
   tests. Minor, and no reason on its own.

### The recommendation

Worth depending on, after the three fixes, and worth writing them up as a
hand-over plan the way §5 did for the other two — that route has now worked
twice. Do not take the dependency before the install rules exist, because
doing so gives back the thing §5 spent its effort buying. And keep it
**optional**: a `GSHTRANS_WITH_INTERPOLATION` component, so that the core
transform never acquires a dependency it does not need. Nothing in the
angular library wants a one-dimensional interpolator; only the layered half
does.

---

## 8. How three-dimensional are the three-dimensional fields?

> *As raised:* there is a question over the extent to which this library
> builds proper 3D fields or just acts as a dumb container for sets of radial
> layers without knowing how they are linked. Flexibility here is likely to be
> important, but we can have minimal functionality in place.

### The honest answer

**It is a container plus one named seam, and everything about how the layers
are linked lives on the far side of the seam.** Precisely, the library knows
four things about the radial axis and no more:

| it knows | where |
|---|---|
| the radii, in increasing order | `RadialGrid` |
| quadrature weights over them, optionally | `RadialGrid`, `IntegrateRadially` |
| that two stacks are on the same radial grid | `RadialGrid::Identity` |
| that a radial line can be gathered contiguously and handed to a callable | `ApplyRadially` |

It does *not* know how to differentiate along `r`, how to interpolate between
layers, whether the nodes form elements, whether a radius is repeated at a
material interface, or what happens at either end. `Gradient` does not compute
`∂_r`; it takes a `RadialOperator` from the caller and applies it. `RadialGrid`'s
own comment states the position: element connectivity, the spectral-element
basis, differentiation matrices and any factorisation "belong to the
application that built them".

So: a dumb container, but not accidentally — deliberately, with the seam named
and documented, which is a different thing from not having thought about it.

### Why that was right, and where it now costs

It was right because the radial discretisations in view genuinely differ in
kind. A finite-difference derivative is banded, a spectral-element one is
block-diagonal, and a caller who has factorised an operator wants to apply the
factorisation rather than hand over a matrix. `field-algebra-plan.md` §17.2
chose a callable over a matrix for exactly that reason, and it has held.

What it costs is that **the library cannot offer anything that needs to know
how layers link.** Three consequences, and the first two are already live:

- **The interpolation of §9 cannot be done in `r`** without knowing whether it
  is legitimate to interpolate across a given pair of layers. A model with a
  discontinuity at the core-mantle boundary has two different values at the
  same radius, and a spline that smooths across it is silently wrong.
- **Every user supplies a radial derivative**, including for the ordinary
  cases. `Gradient` is unusable without one, and the tests reach for a
  hand-written `PowerDerivative` to have anything to pass.
- **There is no 3-D Laplacian, and cannot be one** at this level, because
  `∇²` needs `∂_r` twice.

### The middle position, which is what to build

Keep the seam. Add two things, neither of which takes a decision away from the
caller.

**A. A small library of ready-made radial operators.** Each is nothing but a
type satisfying `RadialOperator`, offered in a header of its own that nothing
in the core includes:

```
FiniteDifferenceDerivative<Order>(radii)   second and fourth order, one-sided at the ends
LagrangeDerivative(radii)                  the differentiation matrix on the given nodes
SplineDerivative(radii)                    from §7's cubic spline
ElementDerivative(radii, elements)         block-diagonal, GLL nodes per element
```

Offered, not imposed: a caller with their own operator passes their own, and
nothing here changes. This is "minimal functionality in place" almost exactly
as raised, it is a few hundred lines, and it makes `Gradient` usable out of
the box, which today it is not.

**B. Let `RadialGrid` optionally carry the element partition.** This is the
one piece of *structure* worth adding, and it is the smallest fact that
distinguishes a discretisation from a list of numbers: which radii belong to
which element, and therefore where a repeated radius is a genuine interface
rather than an error.

It is worth singling out because it is the fact that **more than one thing
needs and nothing can infer**. Interpolation must not cross an interface;
`ElementDerivative` needs the blocks; the constructor's existing
`is_sorted` check currently *rejects* a repeated radius, which is precisely
how a two-sided material interface is represented. Carrying it optionally
costs one `std::vector<Int>` and leaves a caller who has no elements exactly
where they are now.

Everything else — connectivity beyond that, basis functions, boundary
conditions, factorisations — stays the application's, and should.

### What this is not

It is not a proposal to make the library own a radial discretisation. The
distinction worth holding on to is between *knowing the mesh* and *owning the
method*: A and B give the library enough to know the mesh, and leave every
method a caller might disagree about on the far side of the seam where it is
now.

---

## 9. Interpolating a field, as a callable of the two angles

> *As raised:* providing the various fields with an interpolation method that
> returns a callable object of the two angles. We could pass a "method"
> variable which determines the scheme used, e.g. bilinear, bicubic (see
> Interpolation) or direct expansion.

### The shape of it

```cpp
auto f = SpinField<0, Grid>(grid, ...);
auto at = Interpolate(f, Scheme::Bicubic());
auto value = at(theta, phi);
```

`Interpolate` on any spin-weighted node, returning a callable of `(θ, φ)`. For
a tensor field it is per component, or a callable returning the component set;
that is a detail, and the rank-0 case decides everything else.

### The `method` variable is a policy, and there is already a house style

`Policies.h` holds four of these: `Execution`, `Batch`, `Chunking`,
`WignerValues`. Each is a value with named constructors —
`Chunking::ForCache(bytes)`, `WignerValues::Generated()` — rather than a
template parameter or an enum, and for a stated reason: it is a decision the
caller makes about the machine or the problem, not about the mathematics, and
putting it in a type would make the mathematics carry it. An interpolation
scheme is the same kind of thing. So `Scheme::Bilinear()`, `Scheme::Bicubic()`,
`Scheme::Spectral()`, in `Policies.h`, and the fifth member of a set that
already exists.

### The two axes are not alike, and that is the whole design

This is where the problem stops resembling one-dimensional interpolation twice
over.

**φ is uniform and periodic.** So interpolation in φ can be *exact* for a
band-limited field, by trigonometric interpolation — zero-pad the Fourier
coefficients and transform back, which is one FFT the library already owns.
Bilinear in φ would be throwing away accuracy that costs nothing to keep.

**θ is Gauss-Legendre and therefore non-uniform, and neither pole is a grid
point.** Both facts bite. Non-uniform spacing rules out the usual fixed-stencil
bicubic and calls for a genuine interpolant on given nodes — which is §7's
`CubicSpline` or `Lagrange` along a colatitude line. And the poles are outside
the convex hull of the nodes, so *every* local scheme extrapolates there. Near
the poles a spin-weighted field also has the `sin^{|m|}θ` behaviour that makes
polar truncation possible (`core-plan.md` §10), so an interpolant that ignores
it will be worst exactly where it is least defensible.

**Direct expansion is exact and pole-safe, and expensive per point.**
Evaluating `Σ f^N_{lm} Y^N_{lm}(θ, φ)` needs the whole `d^l_{Nm}(θ)` column,
which is `O(lMax²)` work for one point — against `O(1)` for a local scheme.
The recursion for it already exists and is already threaded: `WignerValues::Generated()`
runs it into per-thread scratch and was measured, so the machinery is in
place. The crossover is the thing to measure: at some number of evaluation
points, building a *second grid* and transforming onto it beats evaluating
point by point, and the answer is a straightforward benchmark rather than a
guess.

### What that suggests

Three schemes, and the middle one is the one to reach for by default:

| scheme | φ | θ | exact? | cost per point |
|---|---|---|---|---|
| `Bilinear` | linear | linear | no | `O(1)` |
| `Hybrid` | trigonometric, exact | spline on the GL nodes | in φ only | `O(1)` after setup |
| `Spectral` | exact | exact | yes | `O(lMax²)` |

`Spectral` is the reference the other two are *tested against*, which is worth
more than it sounds: it makes the accuracy of a cheap scheme measurable rather
than asserted, on any field, without an analytic answer to compare to.

### Three things to decide before building

**Lifetime.** The callable reads the field. Whether it borrows or owns is the
same question §7 raises about `Interpolation`'s iterators and the same
question `field-algebra-plan.md` §3.8 answered for expression nodes — and it
should be answered the same way, by value category, so that
`Interpolate(Materialise(...))` does not dangle. This is the one that will
cause a real bug if it is decided casually.

**What it is a callable *of*.** `(θ, φ)` as raised. Worth noting that a
`ScalarFunctionS2` concept already exists in `Concepts.h` and is what field
constructors take — so an interpolant that satisfies it can be fed straight
back into another grid's constructor, which is remeshing in one line and is
probably the commonest use.

**Setup versus evaluation.** A spline along every colatitude line is `O(lMax²)`
of setup that must not be redone per point, so the callable is a built object
and not a lambda over the field. That is an argument for `Interpolate` being a
named type rather than `auto`.

---

## 10. A wisdom mechanism for the computational options

> *As raised:* as we are building various computational options, do we want a
> simple "wisdom"-like method, whereby a user can trial different options for
> their problem? A poor man's version of what FFTW3 does, but maybe of value.

### The case for it is stronger here than it looks

The library already carries five choices whose right answer is
machine-dependent: the thread count, the chunk size, stored versus generated
Wigner values, the FFTW planner flag, and — since §17.7 — bulk repack versus
fused gather. Two more are coming: polar truncation and the transform-major
GEMM layout.

What makes this more than a convenience is that **the plans record, repeatedly,
that these cannot be settled by reasoning.** The chunk optimum at `lMax = 128`
"moves between runs". The batched forward transform gives 2.3× sequentially
and collapses to 1.12× at `lMax = 256` on eight threads. The generated Wigner
path lost on the laptop in every configuration and may still win on
`earth-tunya`, where the aggregate L3 is sixteen times larger. Direction-aware
chunking was worth 2.0× and was found only by measuring. Every one of those is
a case where a user of this library, on a machine neither of us has, will be
running with the wrong setting and have no way to know.

That is exactly the situation FFTW's wisdom exists for, and the analogy holds
in the way that matters: the object being tuned is built once and used many
times, so measurement at construction is cheap against the run.

### What it would look like

The key is a *problem shape*, not a call: `(lMax, nMax, precision, batch count,
direction, threads)`, plus a machine fingerprint. The value is a small tuple of
policy values. The store is a file the user may load at start and save at exit.

```cpp
auto wisdom = Wisdom::Load("gshtrans.wisdom");        // or an empty one
auto grid = Grid::Tuned(lMax, nMax, wisdom);          // measures what it must
wisdom.Save("gshtrans.wisdom");
```

Three properties it must have, each of which is a way it could go wrong:

- **The fingerprint is part of the key.** Wisdom carried to another machine is
  worse than none, and the failure is silent. `earth-tunya` and the laptop
  disagree about nearly everything measured so far.
- **Measurement repeats.** The laptop's noise floor is about 10%, and several
  of the differences at stake are smaller than that. A single timing is not a
  measurement, and back-to-back repetition is what the benchmark harness
  already does.
- **The candidate set stays small.** The knobs are not independent — chunking
  and threading interact, and direction changes the answer — but tuning them
  jointly is a combinatorial explosion for a gain that is mostly in one knob.
  One dimension at a time, in a fixed order, is enough and is honest about
  being a heuristic.

### Where to start, and it is small

**`Chunking::Tuned(...)`, and nothing else at first.** It is the knob that has
already been measured to matter most — getting it wrong cost 2.0× on the
batched inverse, and the right rule was not guessable — it is a single scalar,
and it is already "a property of the machine rather than of the call", which
is where the grid takes it. Timing three or four candidates at grid
construction costs a fraction of building the Wigner table it sits beside.

That gives the mechanism a first customer and a measurable answer before any
persistence, keying or fingerprinting exists. If it earns its place, the store
generalises around it; if it does not, nothing has been built that has to be
maintained.

### What it is not, and the honest caveat

It is not an autotuner, and it should not try to become one. FFTW searches a
space of plans it generates itself; this would time a handful of named
alternatives. The name "wisdom" is worth borrowing for the *persistence* — the
idea that a machine's answer is worth writing down — and not for the search.

And the caveat that applies to the whole idea: none of it substitutes for the
target-machine run that `core-plan.md` §10 item 2 is waiting on. A tuner
measures which of the options we have is best on a given machine. It does not
tell us whether the option set is the right one, and three of the open
performance questions are of the second kind.

---

## Suggested order

Rewritten 2026-08-22. The original list had the dependency work and the
project structure at the top, and both are done; three-dimensional fields were
last and are done too. What is left divides cleanly into things that are ready
and things that are waiting.

**Ready now, in order:**

1. **Tangential tensor fields** (§1). Self-contained, the storage argument is
   strong at rank 4, and the closed intrinsic derivative is the part that
   makes it more than a saving. The change is one generalisation --
   `MultiIndex`'s alphabet becoming a per-slot policy -- and everything
   downstream is already generic over whatever indices exist.
2. **Specialisations** (§3), immediately after, so that the tangential forms
   are named in the same pass rather than bolted on.
3. **Ready-made radial operators** (§8A), which is what makes `Gradient`
   usable without the caller writing a differentiation matrix first, and is
   the smallest useful answer to the question §8 asks.
4. **The element partition on `RadialGrid`** (§8B), which §9 and §8A both
   need and neither can infer.

**Waiting on something:**

- **`Interpolation`** (§7) waits on install rules upstream, exactly as
  GaussQuad and FFTWpp did. Write the hand-over plan now -- that route has
  worked twice -- and take the dependency after.
- **Field interpolation** (§9) waits on §7 for the θ interpolant, though
  `Scheme::Spectral()` needs nothing and could be built first as the reference
  the others are tested against.
- **The wisdom mechanism** (§10) waits on nothing technically, but its first
  customer should be `Chunking::Tuned` alone, and the case for the wider
  mechanism is better made after the target-machine run than before it.

**Independent of all of it:** the **3-j work** (§6), whose priority depends
entirely on when coupling coefficients are actually wanted. If it is wanted at
all, the orthogonality test goes in first: it is a few lines, it needs no
reference implementation, and it turns an unknown boundary into a known one.
There are still two Wigner 3-j codes in this repository and neither is
exercised by the suite.

**Not on the list, and deliberately:** `core-plan.md` §10's efficiency work --
the transform-major GEMM restructure and polar truncation. Those are the
subject of that document's own ordering, they are gated on measurements from
`earth-tunya`, and §10's standing judgement is that the field layer's own
problems are worth more per unit effort than another factor of two on the
transform.
