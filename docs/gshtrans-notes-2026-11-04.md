# GSHTrans in sl3dpp: how it is serving, and what would help

*Running notes, started 4 November 2026 while building the 3D elastic
extension of sl3dpp (lateral δμ, δκ in the static and stepped problems
and their adjoints). To be extended through the parallelisation and
optimisation pass. The question David asked: does GSHTrans do what this
application needs, and where could it be improved for it?*

## 1. What sl3dpp actually uses

Everything goes through `SurfaceGrid` (sl3dpp's handle) onto one
`GaussLegendreGrid` with spins to ±2 and 1.5× oversampling:

- **Batched scalar and spin transforms** (`SynthesiseBatch`,
  `ExpandBatch`, `SynthesiseSpinBatch`, `ExpandSpinBatch`) with the
  radial memory nodes as the batch: the whole pseudospectral pipeline
  (strain → grid → pointwise weight → coefficients → force) is three
  batched spin transforms each way per pass, plus one scalar pair for
  the dilatation under a lateral bulk modulus. This is the dominant
  cost of a 3D step and of the static lateral solve (at lmax 64 the
  lateral force is the whole of the 2.5 s per conjugate-gradient
  iteration; the banded solves are 0.07 s).
- The plan cache, threading, the loop and matrix kernels (the matrix
  kernel is opt-in in sl3dpp: a wash at lmax 64 after the upstream
  livelock fix, as recorded in `docs/dev/lateral-coupling.md` §6).
- `Integrate`, `FromFunction`, `Expand`/`Synthesise` for scalar
  surface fields (the sea level engine's grid is a separate, band-exact
  scalar grid).
- Several grids at reduced degree for depth-truncated node blocks
  (`transform_degrees`): the batched call cannot mix degrees, so
  sl3dpp partitions the nodes into blocks, one grid each.

Not used: the lazy spin-field algebra and `TensorField` (the pipeline
batches across ~10³ radial nodes, which the field layer's per-tensor
batching does not do; the reality reduction and the canonical-component
algebra are sl3dpp's own, fixed by the 1D limit), the layered fields,
the 3-j symbols.

## 2. What works well

- The batched transform is exactly the right primitive: one call per
  spin component with every node in the batch, bitwise equal to the
  per-field calls, threaded internally. The 10× from naive loops to
  batched calls (lateral-coupling.md §6) came entirely from this.
- The shared immutable grid: one Wigner table per resolution, cheap
  handles, so the stepper, the static solver, the adjoints and the
  tests share it.
- Determinism: the degree loops are bitwise independent of the thread
  count, which the reproducibility gates rely on; the transforms agree
  to round-off (each thread plans its own FFTs).

## 3. Where it costs us, and what would help

Observations so far (to be quantified in the optimisation pass):

1. **Real fields through complex transforms.** Every field sl3dpp
   transforms is real (the m ≥ 0 reduction and the ±N reality pairing
   are done on sl3dpp's side). A spin transform that takes a real
   tensor's stored components (spins −2, −1 complex, 0 real) and
   produces the m ≥ 0 coefficients directly would halve the work of the
   spin-0 path and remove the reconstruction loops at the force exit.
   GSHTrans has the reality reduction in its field layer; what the
   pipeline needs is that reduction *on the batched raw transform*.
2. **Depth-variable degree in one batch.** The truncation rule gives
   each radial node a maximum degree; sl3dpp blocks the nodes by degree
   and transforms each block on its own grid. A batched transform that
   accepts a per-item degree (or a ragged batch) on one grid would
   remove the block bookkeeping and the duplicate grids.
3. **Memory traffic, not flops.** The pass is bound by moving the grid
   tensors (0.9 GB each at lmax 64, 4.9 GB at 128): synthesise into a
   buffer, multiply in place, analyse out. A transform entry point
   that fuses the pointwise weight (synthesise → multiply by a real
   field → analyse, never materialising the grid tensor beyond one
   latitude ring) would cut the traffic by about 3× and the peak memory
   by one grid tensor per branch. This is the single most valuable
   change for this application.
4. **The Wigner table** is the gigabytes at high lmax (2.6–5.2 GB at
   lmax 256). On-the-fly generation exists as a policy; the trade-off at
   the batch sizes sl3dpp uses (10³ nodes per call, so the table is
   read once per node) has not been measured here yet.
5. **Distribution.** GSHTrans is single-process. sl3dpp distributes
   over radial segments (each rank holds its nodes' grids, and the
   transforms are local), so the transform never needs MPI; what the
   static solve and the stepper under lateral moduli need is only that
   the batched calls work per rank on its own nodes, which they do.
   Distribution over the sphere (latitude rings) would be a different
   design and is not needed while the radial partition balances.
6. **Oversampling control.** The 3/2 rule is applied globally
   (`ForBand`); the products are quadratic so this is right, but a grid
   whose longitude count is chosen for FFT efficiency at a given
   oversampling (not just ≥ 2 lmax + 1) would be a small, free win.

## 3a. First numbers (4 November 2026, laptop, 16 threads)

One pass of a 3D step at lmax 64 on prem8 with a Maxwell mantle
(`TestLateralModuli3D.DISABLED_StepTiming`):

| what the pass carries | memory nodes | s per pass |
|---|---|---|
| lateral viscosity (3 spin pairs of transforms) | 1108 | 0.66 |
| + lateral shear over the whole solid | 1700 | 0.96 |
| + a four-companion scalar pipeline (the lateral density, since removed) | 1700 | 1.8–3.3 (varies run to run on the laptop) |

At lmax 32 the same three are 0.075, 0.10 and 0.19 s: the pass scales
as ~9× for a doubling of lmax (nodes ×2, grid ×4), transform-bound. The
static lateral solve at lmax 64 is 2.5 s per conjugate-gradient
iteration, all of it the lateral force (the banded solves are 0.07 s).

What these say for GSHTrans:

- A scalar transform costs about as much as a spin transform per field
  here (the density's eight scalar transforms cost more than the six
  spin ones), so there is no cheap scalar path to lean on; batching the
  four companions into one call would save call overhead, not flops.
- The fused product transform (item 3 above) is the lever: every
  companion and every strain component goes grid → weight → coefficients
  with the full grid tensor materialised in between.
- The node count is the other lever, and that is sl3dpp's (the active
  elements rule, built the same day); a ragged batch (item 2) would let
  the depth truncation cut further.

## 3b. The conjugate-gradient step (5 November 2026)

Closing a step by conjugate gradients instead of the fixed point does
not change what GSHTrans is asked to do: an iteration is one pipeline
(synthesis, pointwise weight, analysis) like a pass, so the transform
cost per iteration is unchanged and the iteration count falls (the
ladder in sl3dpp/docs/dev/lateral-coupling.md §5a''). What it does
change is the balance: at small lmax the band products and the Krylov
algebra show (an iteration ~1.8 passes at lmax 16), at lmax 64 the
transforms dominate both, so everything in §3 about the pipeline
applies as before, and a fused product transform would cut the CG
iteration exactly as it cuts a pass.

## 3c. The field layers: would they serve, and what would make them (5 November 2026)

David's question: the pseudospectral products — do we use the lazy
spin-field algebra and the tensor machinery, and if not, why not, and
could they be made useful? Read against the headers at v1.0.0+8
(SpinField/, Tensor/, Layered/, Expansion/).

**What sl3dpp hand-codes today.** `GridCoupling` holds the deviatoric
strain as its own five-real layout per grid point (two complex
components at N = ±2 and ±1, one real at N = 0, the reality reduction
applied by hand), stacked radius-major over ~10³ memory nodes, and the
step's product is a loop: `Multiply(field, weight)` scales all five
arrays by the pointwise weight, `MultiplyScalar` a companion. The
transforms are the grid's batched calls on raw spans with the radial
nodes as the batch (`SynthesiseSpinBatch`/`ExpandSpinBatch`, three
spin pairs a pass; one scalar pair per companion). Nothing of the
field layers is touched.

**What the layers are, against that.**

- *The lazy spin-field algebra* (`SpinField`, the `Binary`/`Unary`
  nodes, the overloads) is pointwise and index-preserving by
  construction: a product of an N = 2 field with an N = 0 weight is a
  node of N = 2 that reads each operand at the same point, and the
  reality and index rules are checked at compile time. That is exactly
  our product. And `Expand(node)` does NOT materialise: it calls the
  node's `EvaluateInto` straight into the transform's input vector and
  transforms that — one pass over memory, the fused product-transform
  of §3 at the single-field level, already there. What it lacks is the
  batch: `Expand` takes one angular field, allocates its own value
  vector and its own expansion per call, and runs one transform of
  count 1. Our pipeline is 10³ such fields per component; calling it
  10³ times forfeits the batched Wigner evaluation, which is the whole
  cost structure (§3). So the algebra is right and unusable at our
  scale as it stands.
- *The layered fields* (`LayeredSpinField`, `LayeredTensorField`,
  radius-major stacks, each component its own batch) are precisely our
  layout — the header says as much, "the layout the applications this
  exists for already use" — and their `Expand`/`Evaluate` run the
  batched transform with the stack as the batch, which is our call.
  But they take an OWNING stack, not a node: there is no layered
  expression type, so a product has to be written into a stack before
  it can be transformed — the extra traversal we would want to remove,
  the same two-pass shape our hand loop has now. The lazy algebra and
  the batching are each present and do not compose.
- *The tensor layer* (`TensorField`, `Orbits`, `TensorExpansion`)
  derives the stored components of a real symmetric rank-2 tensor
  from the symmetry and the reality reduction — the same five reals
  per point we lay out by hand, with the pinned-real N = 0 component
  and the conjugate pairs made a type property rather than a
  convention, and the (q, j) assembly a compile-time traversal. For a
  rank-2 deviatoric tensor the gain is clarity and the compile-time
  check of what we have fixed by the 1D limit; it buys no speed. Where
  it would matter is a rank-4 object: the relabelling plan
  (~/dev/planning/relabelling-sl3dpp-2026-11-04.md) weights the strain
  by a 21-coefficient tensor pointwise, and the layer's `PointMajor`
  layout and compile-time contraction over stored orbits is built for
  "a rank-4 tensor applied pointwise to a rank-2 one — 81 numbers at
  one point and none at any other". Hand-coding that is where the
  bugs would be; the layer would also make the fluid/solid alphabet
  (`BundleMaps`) explicit.
- *`LayeredGradient` / `RadialOperator`*: the 3D gradient with the
  angular half (grad_1, spectral) from the library and the radial
  derivative from the application's discretisation. Our strain is not
  formed that way — the degree operators give the strain components of
  each (l, m) solution directly from the radial basis (the canonical
  components of the 1D form) and the grid pipeline only transforms
  them — so the gradient machinery is not on our path. It would be on
  the path of a relabelling-style referential formulation, where
  fields are differentiated on the grid.
- *`RadialMajor`* (the [(l, m)][r] layout for radial solves touched
  many times, with the repack) is what the preconditioner of the
  conjugate-gradient step does by hand (`StrainCoefficients` per
  degree, the banded solves per (l, m)); we gather per degree from the
  `SpectralTensor` and it works, but it is the same repack written
  twice.
- *The 3-j symbols*: the mode-coupling route we do not take (products
  are formed on the grid); no use.

**What would make them useful here, in order of payoff.**

1. A layered `Expand` that accepts a layered NODE — a pointwise
   product of a `LayeredSpinField` with a `LayeredSpinField<0, Real>`
   weight, or a `Binary` over slices — and evaluates it slice by slice
   into the batch buffer before one batched transform. That is the
   fused product-transform of §3 and it is a small addition: the
   two-dimensional `Expand` already evaluates a node into its buffer;
   the layered one needs the same loop over `nR` slices into a batch
   buffer the caller may supply. With it the step's pass is
   `Expand(weight * Evaluate(coefficients))` per component, one
   traversal fewer, and sl3dpp's `Multiply` loops go away.
2. Caller-supplied buffers on the layered transforms (an `Into` form):
   the per-pass allocation is already a measurable cost at lmax 64
   (§3a), and our pipeline keeps its buffers across passes.
3. The ragged batch of §3 (nodes at different truncation degrees in
   one call), which the layered type could carry as a per-slice degree.
4. Then, and only then, moving `GridCoupling`'s tensor onto
   `LayeredTensorField<2, Symmetric, Real>`: the stored-orbit set is
   ours, the transforms are ours, and the five hand-written arrays
   become one typed object — worth it when the rank-4 weight of the
   relabelling work arrives, not before.

Not worth pursuing: the per-field (unbatched) algebra for the step
(10³ transforms of count 1), and the mode-coupling 3-j route.

## 4. Questions to settle in the optimisation pass

- Measure the pipeline's split between transform flops and memory
  traffic at lmax 64 and 128 on the laptop and on the server, with the
  loop and matrix kernels; decide whether (3) is worth proposing
  upstream.
- Try the on-the-fly Wigner policy at lmax 128 with 10³-node batches.
- Whether the conjugate-gradient step (replacing Picard/Anderson in the
  stepper) changes the balance: fewer passes but each a full transform
  pair still — the transform remains the unit of cost either way.
