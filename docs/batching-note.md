# Design note: general batching

For Phase 11 of `fix-plan.md`. New work, not a fix; to be agreed before
anything is built. Every figure was measured on 2026-09-19 on the development
laptop (16 threads, OpenBLAS held to one thread), from an optimised build, with
the prototype in a scratch copy. Nothing in the repository has changed.

There are three decisions, **B1–B3**, each with a recommendation.

---

## 1. What is already there

`Batch` in `Policies.h` is FFTW's advanced interface: a count, a stride and a
dist, given independently for the input and the output of

    grid.ForwardTransformation(lMax, n, in, inBatch, out, outBatch, policy);

and its inverse. Element `j` of field `k` is at `j * stride + k * dist`. The
radial case already sits on it — `LayeredSpinField::Batch()` is
`Batch::Contiguous(nR, FieldSize())` — and so does the flat tensor, with
`Interleaved` for its point-major layout. So any *affine* layout works today,
including a caller's own `[r][point][component]` or `[point][r][component]`
array, by taking a subspan at the component's base.

Two things are not there. One is a layout that is not affine. The other is
already expressible and merely unspelt.

## 2. Fields at known offsets

### 2.1 The proposal

    static Batch Batch::At(std::vector<Int> offsets, Int stride, Int size);

A batch whose `k`-th field starts at `offsets[k]`, with `stride` between its
successive elements, each field `size` elements long. And, since the commonest
reason to want one is "these fields out of that batch",

    Batch Batch::Subset(std::span<const Int> which, Int size) const;

so that `field.Batch().Subset(solidRadii, field.FieldSize())` is the solid
regions of a layered field and nothing else.

Uses: a subset of radii (the solid regions only; a mask; every radius a
preconditioner touches); finite-element storage with padding between elements
or duplicated interface nodes; the several same-`N` components of a tensor
across radii in a caller's own array.

### 2.2 Why it is cheap

Because of how the kernels were written. They touch a layout *only* through
`Count`, `Stride`, `Offset(j, k)`, `Span` and `Disjoint` — fourteen call sites
— and always gather into chunk-local scratch before FFTW or the BLAS sees
anything, so neither ever meets the caller's layout. `Offset` becomes
`j * stride + base(k)`, and **no kernel changes**.

Prototyped and timed, 32 complex fields at lMax = 256, forward / inverse in ms,
best of nine, two interleaved rounds (the second in brackets):

| | loop, 1 thr | loop, 8 thr | matrix, 1 thr | matrix, 8 thr |
|---|---|---|---|---|
| today | 231 / 198 (235 / 195) | 146 / 67 (152 / 70) | 69 / 85 (69 / 86) | 21.5 / 23.3 (21.9 / 23.4) |
| prototype, affine batch | 201 / 200 (206 / 205) | 140 / 70 (140 / 74) | 67 / 86 (69 / 88) | 21.2 / 23.5 (24.0 / 24.9) |
| prototype, `Batch::At` | 201 / 198 (204 / 202) | 141 / 71 (144 / 72) | 67 / 90 (68 / 87) | 21.5 / 23.6 (22.8 / 25.8) |

The affine path is unchanged — the branch in `Offset` is on a value that is
invariant over every loop it sits in, and the compiler hoists it — and an
offset batch is within a few per cent of a contiguous one. (The 231 → 201 on
the sequential loop forward is the code-placement effect `lessons.md` records
for exactly that row, and not a gain.)

### 2.3 What has to be decided

**Who owns the offsets.** Two honest answers.

- *The batch does*, as a shared immutable table. A `Batch` stays a value that
  can be copied, stored and returned, and a dangling table is impossible. The
  cost is an atomic count when a `Batch` is copied, which happens a dozen times
  a transform and never per element — the rule `lessons.md` now states.
- *The caller does*, and the batch holds a span — FFTW's guru interface works
  this way. Free, and one more lifetime to get wrong; and what going wrong looks
  like here is a transform that *writes* through stale offsets, in silence.

**Disjointness** stops being arithmetic and becomes a sort: two fields overlap
iff their offsets agree modulo the stride and differ by less than
`size * stride`. It is checked once, in `At`, which is why `At` takes the size.
A call then checks only that the fields it is given are no longer than the
batch was built for.

**Equality** compares offsets by value, so that two batches describing the same
layout are equal whichever table they hold.

> **B1. Build `Batch::At` and `Batch::Subset`?** Recommended: yes.
>
> **B2. Who owns the offsets?** Recommended: the batch. This review was mostly
> about things that went wrong in silence, and a span's failure mode is
> exactly that; the cost of avoiding it cannot be measured.

### 2.4 Tests

- A transform through `At` with affine offsets equals the affine transform,
  bit for bit, in every kernel and under threads.
- A shuffled, gapped layout round-trips, and the gaps are untouched.
- `Subset` of a layered field's batch transforms those radii and no others.
- Overlap — for stride one and for an interleaved stride — a negative offset,
  and a field longer than the batch was built for, all throw.
- A `Batch` outlives the vector it was made from.

## 3. Transforming straight into radial lines

### 3.1 It already works

`[(l, m)][r]` — the layout `RadialMajor` holds, in which a radial line is
contiguous — *is* `Batch::Interleaved(nR, nR)`. So

    grid.ForwardTransformation(lMax, n, field.Data(), field.Batch(),
                               lines, Batch::Interleaved(nR, nR), policy);

writes radial-major coefficients directly, today, through the public
interface, and the result is **bit-identical** to transforming and then
transposing (checked: the difference is exactly zero). Nothing has to be built
for this to be possible. The question is whether it is worth spelling, and that
is a measurement.

### 3.2 The measurement

A: transform to radius-major, then the tiled transpose (today's
`Expand` + `RadialMajor`). B: the direct transform. No allocation inside
either timed region. B / A, forward and inverse:

| | nR = 64 | 100 | 128 | 200 |
|---|---|---|---|---|
| matrix, lMax 128, 1 thread | 0.83 / 0.84 | 0.83 / 0.87 | 0.85 / 0.87 | 0.94 / 0.89 |
| matrix, lMax 256, 1 thread | 0.96 / 0.98 | 0.93 / 0.97 | 1.01 / 0.99 | 0.94 / 0.99 |
| matrix, lMax 256, 8 threads | **1.08 / 1.10** | 0.79 / 0.78 | **1.18 / 0.99** | 0.82 / 0.81 |
| loop, lMax 128, 1 thread | 0.94 / 1.06 | 0.96 / 0.98 | 0.95 / 1.10 | 0.96 / 1.02 |
| loop, lMax 256, 8 threads | **1.22 / 1.07** | 1.08 / 1.05 | **1.13 / 1.12** | 1.05 / 1.03 |

Three things to read from it.

- **The prize is small.** The transpose is five to ten per cent of a transform,
  so that is the most B can win sequentially, and it does: a few per cent at
  lMax = 256, fifteen at 128.
- **The power-of-two hazard is real and is where it was predicted.** The
  scatter into radial-major order has stride nR, and at nR = 64 and 128 under
  threads B *loses* by eight to twenty per cent — the same cache-set collision
  the Fourier stage guards against, and the reason `RadialMajor`'s transpose
  is tiled. At nR = 100 and 200 B wins by twenty.
- **With the loop kernel under threads B loses everywhere**, by three to
  twenty per cent.

So on time alone this is a wash with a trap in it, and "direct is faster"
would have been the wrong thing to build on.

### 3.3 What it is actually good for

Memory. Route A holds the radius-major expansion *and* the radial-major copy;
route B never makes the first. That is `nR × coefficients × 16` bytes a field:
211 MB at lMax = 256 with nR = 200, and 2.1 GB at lMax = 512 with nR = 500, per
scalar field — which on a production run is the difference that matters, where
ten per cent of a transform is not.

### 3.4 The proposal

Spell it, as the memory-saving route, and leave the choice to the caller:

    auto lines = ExpandToLines(field, lMax, policy);   // a RadialMajor
    auto back  = EvaluateLines(lines, policy);         // a layered field

with `RadialMajor` gaining a `Batch()` (`Interleaved(nR, nR)`) and a way to be
made empty at a given shape, which it has privately already. No automatic
switching between the routes: the table above does not support a rule simple
enough to be trusted, and it was measured on a laptop. The documentation says
what was measured, names the power-of-two case, and says that the route to
take when time matters more than memory is the one a measurement on the
machine in question prefers. `run-server-benchmark.sh` gains the comparison,
so that the 64-core answer exists before anyone needs it.

> **B3. Build `ExpandToLines` / `EvaluateLines` as the memory-saving route,
> with no automatic choice?** Recommended: yes. The alternative is to build
> nothing and document the one-line idiom of §3.1, which is a fair answer too:
> the capability is there either way.

### 3.5 Tests

- `ExpandToLines` equals `RadialMajor(Expand(...))` bit for bit, for a complex
  and a real-valued field, both kernels, sequential and threaded.
- `EvaluateLines(ExpandToLines(f))` round-trips.
- The lines it returns go through `ApplyToLines` and the operator-grid check
  exactly as a transposed buffer does.

## 4. What I would not build

- **A two-level affine form** (`countA, distA, countB, distB`): `At` subsumes
  it.
- **Mixed upper indices in one call.** The Wigner block and the GEMM matrix are
  per-`n`, so the kernel has nothing to share; `TensorField` already is that
  loop.
- **An FFTW-style plan object.** Its attractions are real — validate once, make
  workspaces outside the parallel region, a home for tuning results — but it
  cuts across the recorded decision to keep policies on the grid handle with
  `With()`, and `At`'s construct-time validation takes most of the benefit.
  The thing to revisit if `At` grows.

## 5. Order of work

1. `Batch::At`, `Batch::Subset`, their tests, and the README's batching
   paragraph (B1, B2).
2. `RadialMajor::Batch()`, `ExpandToLines`, `EvaluateLines`, their tests, an
   example, and the benchmark row (B3).

Each is one stage with the full gate run after it.
