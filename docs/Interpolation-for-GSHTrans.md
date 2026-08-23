# Interpolation — what GSHTrans found using it

A hand-over note in the other direction from `GaussQuad-for-GSHTrans.md` and
`FFTWpp-for-GSHTrans.md`: those record what GSHTrans needed to know about a
dependency, and this records what a consumer learned from `da380/Interpolation`
at `main` (`56c2fe9`, 2026-08-22) and what would make it more useful. Everything
below was compiled or measured on this machine rather than read off the
interface.

**Short version: it is in, it is used, and nothing has to change. There is one
request worth considering, one cheap piece of CI insurance, and one convention
we now depend on and would like to keep depending on.**

---

## 1. How GSHTrans uses it now

An **optional dependency, on by default**: `GSHTRANS_WITH_INTERPOLATION`,
fetched-or-found in the same pattern as NumericConcepts, GaussQuad and FFTWpp,
tracking `main` as they do. A build with it off compiles, fetches nothing, and
runs 273 of the 280 tests.

It is used in exactly two places.

- **`Resample`** — carrying a three-dimensional field from one set of radii to
  another. `RadialInterpolation::Linear()`, `::CubicSpline()` and `::Akima()`
  are the policy, so the whole menu costs a policy argument rather than three
  implementations. This is the one runtime use.
- **As an oracle in the tests.** `SplineDerivative` is checked against
  `Interpolation::CubicSpline` on data with no polynomial structure, real and
  complex, to `1e-12`. An independent implementation of the same spline is a
  stronger check than any property GSHTrans could assert about its own.

It is deliberately **not** used for the three ready-made radial derivatives,
which are self-contained. §3 is why, and it is the substance of this note.

> **Postscript, 2026-08-23.** All four asks below were answered, and §3's
> request changed what GSHTrans built. `SplineDerivative` is now assembly over
> `CubicSplineSystem` — the sixty duplicated lines are gone — so it *is* a
> runtime user of this library, and the sentence above holds for the other two
> operators only. `ElementDerivative` joined them since, so the
> dependency-free set is three: finite differences, the differentiation
> matrix, and the element derivative. `docs/field-algebra-plan.md` §21 records
> what else changed, including the measurement in §3 being corrected downwards
> from an order of magnitude to 1.2×–2×.

---

## 2. What was checked and is fine

- **The interpolators take complex ordinates.** `InterpolationRanges` requires
  `RealOrComplexRange` for the ordinates, so a radial line of spherical-harmonic
  coefficients goes in directly with no real/imaginary split. Verified by
  compiling `CubicSpline` and `AkimaSpline` over `std::complex<double>`. This
  matters more than it sounds: the model application does all of its radial work
  in the spectral domain, where every line is complex.
- **The dependency set is `NumericConcepts` and nothing else**, which GSHTrans
  already has. So it adds no transitive dependency at all — the reason making it
  optional was almost not worth doing.
- **Install and export work.** `Interpolation::Interpolation` is exported and
  the install rules are unconditional, so a copy fetched into a GSHTrans build
  installs itself and `find_package(GSHTrans)` still resolves. Two of the other
  three dependencies had to be changed before they could do that; this one
  needed nothing.
- **The lifetime rule matches ours.** Borrow an lvalue, own an rvalue, is
  exactly the operand-storage rule GSHTrans chose for expression nodes and for
  the same reason. Arrived at independently, and it means a value can cross
  between the libraries under one rule rather than two.

---

## 3. The request: the node-dependent part is separable, and only you can
expose it

This is the one that changed what GSHTrans built.

**What happened.** GSHTrans needed `d/dr` of a spline through each radial line
of a field. A radial operator is called **once per line** — `nθ · nφ` in the
spatial domain, or about 66,000 coefficients at `lMax = 256` in the spectral
one — so the natural spelling is a `CubicSpline` per line. Two things make that
more expensive than it needs to be, and neither is a defect:

- construction computes the coefficients and holds them in a `std::vector`, so
  it costs about four allocations; and
- `Evaluate` locates its segment by binary search, so asking for the derivative
  at every node of the curve just fitted costs `nR log nR` lookups.

Both are the right design for what the interface is *for* — evaluating a curve
at a point. They are simply the wrong shape for "give me the derivative at
every node", which is what a differentiation operator is.

**The observation that matters:** for a fixed set of abscissae, the spline
system's **matrix depends on the nodes alone**. Only the right-hand side carries
the ordinates. So a caller with many ordinate sets on one grid — multi-component
data, an ensemble, a field of coefficients — can factorise once and solve many
times, and that is a common shape rather than a GSHTrans quirk.

**Measured**, over 66,049 lines, a spline constructed per line against an
operator that reuses its factorisation:

| `nR` | per line | reused factorisation | ratio |
|---:|---:|---:|---:|
| 9 | 13.9 ms | 6.8 ms | 2.0× |
| 17 | 23.0 ms | 13.8 ms | 1.7× |
| 33 | 58.4 ms | 27.8 ms | 2.1× |
| 65 | 73.4 ms | 61.2 ms | 1.2× |
| 129 | 193.6 ms | 137.4 ms | 1.4× |
| 257 | 401.9 ms | 296.9 ms | 1.4× |

**So it is 1.2× to 2×, and not more.** That is worth being honest about,
because the first version of this argument inside GSHTrans put it much higher
and was wrong. A 1.4× on an operation applied every iteration of a matrix-free
solve is still worth avoiding, but it is not the headline.

**The headline is the duplication.** Because there was no way to reach the
factorisation, GSHTrans now contains its own natural-cubic-spline system and
its own Thomas sweep. That is about sixty lines that exist only because
`Interpolation::Detail::SolveTridiagonal` is in `Detail` — a function that is
fully documented, typed so a real matrix can act on a complex right-hand side,
and carries its own no-pivoting argument from diagonal dominance. It reads like
considered public API that happens to be namespaced private. Two implementations
of one spline, in two repositories with one author, is the cost actually worth
avoiding.

**Two ways to fix it, smallest first.**

1. **Promote `SolveTridiagonal`.** One namespace change. GSHTrans would delete
   its Thomas sweep and keep its own assembly of the system.
2. **Expose the factorisation itself** — something like
   `CubicSplineSystem<Real>(abscissae)` with `Solve(ordinates, secondDerivatives)`,
   which `CubicSpline` would then use internally. GSHTrans would delete all
   sixty lines. This is the better shape for the library too: it names the thing
   that is reusable, and it makes "one grid, many datasets" a first-class case
   rather than something a caller has to notice.

A third, orthogonal and smaller: **an "evaluate at the nodes" path** that
returns the value or derivative at every abscissa without a binary search per
point. That is the other half of what a differentiation operator wants.

---

## 4. Cheap insurance: a GCC 13 leg in CI

CI covers GCC 14 and Clang 18. The deployment target for the codes GSHTrans
serves — `earth-tunya` — has **GCC 13.2**, and current `main` compiles there:
verified with `g++-13 -std=c++23`, including `CubicSpline`, `AkimaSpline` and
`LagrangeBasis` instantiated over `std::complex<double>`.

So this is locking in something that already works rather than asking for a
fix. The roadmap's decision to use deducing `this` — which GCC 13 lacks — is
the thing that would break it, and no use of it survived into the headers; a CI
leg is what would keep that true by accident rather than by inspection.

---

## 5. A convention we now depend on, deliberately

`Piecewise` is exactly what it says, checked against current `main`:
breakpoints strictly increasing and one longer than the pieces, the pieces
tiling with no gaps, continuity across a breakpoint **not** checked, evaluation
right-continuous so piece `k` owns `[b_k, b_{k+1})`, and `Limits(x)` returning
both one-sided values.

GSHTrans has just adopted all of that. Its `RadialGrid` can now carry an element
partition, and the representation was chosen to *be* breakpoints in that sense:
blocks are disjoint, two elements meet at a repeated radius, both one-sided
derivatives exist there, and continuity is the caller's business rather than
something the grid enforces. The reasoning is in
`docs/field-algebra-plan.md` §20.

The point of saying so: the two libraries should not be able to disagree about
what happens at the core–mantle boundary, and that agreement is now load-bearing
for a consumer rather than a coincidence. If the right-continuity convention or
`Limits` ever moves, GSHTrans wants to know.

---

## 6. Not your problem

Recorded so it is not mistaken for a request. `Interpolation` requires CMake
3.24. GSHTrans declares 3.20 while already using `FIND_PACKAGE_ARGS`, which is
itself a 3.24 feature, so GSHTrans's stated minimum has been wrong since before
this dependency existed. That is a one-line fix on the GSHTrans side and no
constraint on yours.

---

## 7. Summary of what is being asked

| | ask | size | why |
|---|---|---|---|
| 1 | Expose the spline factorisation, or at least promote `SolveTridiagonal` | small | deletes a duplicated implementation; 1.2–2× where it is used |
| 2 | An evaluate-at-the-nodes path | small | avoids a binary search per node for the case that has none to do |
| 3 | A GCC 13 leg in CI | trivial | a known deployment target, already passing |
| 4 | Keep `Piecewise`'s conventions | none | a consumer now depends on them |

Nothing here blocks GSHTrans. The dependency is in, on by default, and doing
useful work as it stands.
