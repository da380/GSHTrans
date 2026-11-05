# Known limitations

What the library deliberately does not do, listed so that it is not mistaken
for a bug. The reasons for the larger choices are in `design.md` and in the
reference's "Scope and limitations".

**Numerics and transforms**

- **Degree ceiling.** `lMax ≤ SafeDegree<Real>(mMax, nLow, nHigh)` — 1822 in
  double for every order at `n = 0`, higher for restricted orders, 913 at
  worst — set by underflow of the Wigner seed (`design.md` §5.1). Lifting it
  needs an exponent carried per column (Fukushima, *J. Geod.* 86, 2012), which
  is not built; long double reaches 30 740, and a table of a few orders goes
  much higher in double already.
- **6-j symbols are not implemented.** Schulten and Gordon give the recursion
  in the same paper, and the 3-j machinery would carry over.
- **The 3-j `Fill…` functions allocate one row of scratch per call.** Calling
  them in a hot threaded loop would want the scratch passed in.
- **ILP64 BLAS is not detected.** `Blas.hpp` assumes 32-bit integers; an ILP64
  BLAS would link and pass wrong arguments. No portable header-side check
  exists.
- **No persistence for tuning results.** Tuning is cheap at start-up up to
  lMax of a few hundred; around lMax = 512 building the tables to measure
  them stops being cheap, and a persisted store would become worth having.
- **Reduced-precision storage of a double table** is not offered; any
  experiment would first need a tighter round-trip oracle.

**Tensors**

- **One slot alphabet per tensor.** Mixed objects such as Dahlen & Tromp's
  T^{rΩ} cannot be one tensor.
- **No trace-free symmetric tensors.** The orbit machinery handles signed slot
  permutations only, so the trace of a tangential symmetric (spin-2) tensor
  cannot be removed from storage; the caller subtracts it.
- **`Materialise<Symmetry, Reality>` does not check** that the value has the
  symmetry or reality stated: only representatives are copied, so the
  difference is discarded. A check would need every derived component and a
  tolerance.
- **Tensor `Materialise` has no contiguous fast path**; it copies element by
  element even into a component-major target where `EvaluateInto` on a
  contiguous span would do.

**Layered fields**

- **Radial operators are conveniences.** They refuse what they cannot do
  correctly (narrow elements, Lagrange across interfaces) rather than gaining
  the capability (`design.md` §6).
- **Spline end conditions apply at every interface.** `SplineDerivative`'s
  `left` and `right` are imposed on every element, so the default `Natural`
  forces zero curvature on both sides of each interface. Interior pieces
  cannot be given different conditions.
- **`Gradient` assumes a radial operator with real weights.** Derived
  components are reconstructed by conjugation, which commutes with the
  operator only if its weights are real. Every ready-made operator's are; a
  caller's complex-weighted one would give wrong derived components silently.
- **The layered `Gradient` refuses an operand too short to hold its result.**
  Its radial half maps stacks to stacks of the same shape, so it cannot change
  the degree as the flat `SurfaceGradient` does.
- **The layered `Gradient` threads only its radial half**, and
  `SurfaceGradient` takes no `Execution` policy.
- **No automatic choice between `ExpandToLines` and `Expand` + `RadialMajor`.**
  The timing difference depends on `nR` (power-of-two strides) and threads in
  a way no simple rule captures.
