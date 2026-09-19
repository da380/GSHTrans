# Design note: 3-j symbols and Wigner tables at large degree

For Phase 6 of `fix-plan.md` (review findings W1, W2, W3). To be agreed before
anything is built. Every figure below was measured on 2026-09-19, on the tree
at Phase 5, with harnesses that are described where they are used; the
prototypes were made in a scratch copy and nothing in the repository has
changed.

There are four decisions, **Q1–Q4**, each with a recommendation.

---

## 1. 3-j symbols (W1, and the 3-j half of W3)

### 1.1 What is wrong

One cause, two symptoms. `SchultenGordonRow` recurses inward from both ends of
a row and joins the halves. Where to stop the forward half is decided by
SLATEC's rule — *stop when |c1| = |B/A| first rises* — which is a proxy for
"the values have stopped growing". In a row with a classically allowed region
the proxy is good: |c1| is smallest near the middle of that region. In a
near-stretched row there is no such region, the row is a single hump, and the
minimum of |c1| is not at the hump. The forward half then stops short, or
overshoots, and one of the halves runs *downhill* for many steps — the
direction in which the unwanted solution grows.

- **Symptom 1, accuracy.** Tens of steps downhill cost digits: errors of
  1e-10 where 1e-14 is available.
- **Symptom 2, the throw.** When the two halves end up more than 1e154 apart
  in scale, `ratio * ratio` overflows, the row normalises to zero, and the
  residual check throws.

### 1.2 How bad it is

Against `tests/RacahReference.h`, which is exact to rounding where its sum is
short — and a short sum *is* the stretched region, reached by a cyclic
permutation of the columns. Errors are relative to the largest entry of the
table. 77 tables per sweep, l = 2 … 1000.

| sweep | throws today | worst error today |
|---|---|---|
| (l, 2l, l) | **44 of 77** | 1.4e-11 |
| (l, 2l−1, l) | 8 of 77 | 2.0e-11 |
| (l, l, 2l) | **44 of 77** | 9.7e-12 |
| (l1, l1+l3, l3), 150…400 | 17 of 121 | 3.9e-10 |
| `float`, (l, 2l, l), l ≤ 200 | 7 of 29 | 1.4e-4 |

So more than half of all stretched tables up to l = 1000 cannot be built, and
`Wigner3jStack(l1, l3)` includes the stretched l2 = l1 + l3 by default.

Away from the stretched region, double against long double on 300 random
triples with l ≤ 500: 12 exceed 2000 ε, worst 8.2e-10 at (10, 430, 423).

`float` is worse than the review said. On 300 random triples, `float` against
`double`: 1 throws and **7 are wrong by more than 1e-3 of the table's scale,
the worst by 0.60**, all passing the residual check — for instance
(74, 466, 392) is off by 0.33.

### 1.3 Proposal A, prototyped

Two changes to `SchultenGordonRow`, about fifteen lines.

1. **Stop on the quantity itself.** Run forward while |g| grows and stop at the
   first decrease, `k ≥ 2 && |g[k]| < |g[k−1]|`. That is the edge of the
   allowed region when there is one and the top of the hump when there is not,
   and in both cases the three matching points sit where the values are
   largest, which is where a least-squares match is best conditioned. The
   backward half then climbs all the way to the join.
2. **Join without squaring the ratio.** The total is formed as
   `hypot(ratio·√sumForward, √sumBackward)`, which cannot overflow where the
   entries themselves do not.

The phase tracking is untouched — `lessons.md` records what it cost to learn.

| | today | **A** |
|---|---|---|
| throws, all sweeps above (double and float) | 124 | **0** |
| worst error vs Racah, (l, 2l, l) | 1.4e-11 | 3.6e-12 |
| worst error vs Racah, (l1, l1+l3, l3) | 3.9e-10 | 7.1e-13 |
| random triples exceeding 2000 ε | 12 | 3 |
| worst random triple | 8.2e-10 | 8.9e-12 |
| fat triangles (l, l, l), l ≤ 500, vs long double | 3.2e-14 | 1.0e-14 |
| `float` vs `double`, random: wrong by > 1e-3 | 7, worst 0.60 | **0**, worst 1.3e-4 |
| `float` at (74, 466, 392) | 0.33 | 5.6e-6 |
| cost, fixed set of tables | 0.183 s | 0.181 s |

What is left at the 1e-11 level is not the algorithm's doing, and I checked
both residues. The 3.4e-10 on (l, 2l−2, l) is identical in every version and is
the *oracle*: there the library's double and long double agree to 1.3e-14, and
the Racah sum has three cancelling terms built on a double-precision `lgamma`.
The 8.9e-12 at (0, 478, 478) is a flat row of 957 equal entries, where any
three-term recurrence accumulates about n² ε; it is the same before and after.

I also prototyped a combined rule — SLATEC's, overridden only while the values
still climb in a forbidden region, or once they fall in one. It removes the
throws equally, and is marginally *worse* than A on fat triangles and in
`float`. A is simpler and at least as good everywhere measured, so the
combined rule is dropped.

**Proposal B (Luscombe–Luban ratio recursion) is not needed.** It was the
fallback if A failed the table. A passes it.

> **Q1. Adopt proposal A?** Recommended: yes.

### 1.4 The residual tolerance

`ResidualTolerance` is 1000 (n+1) ε relative to the row's scale — 0.094 for a
long `float` row, which is how a row wrong by a third passes it. With A the
rows are right, so the check has nothing to catch in the cases measured; but a
check that cannot fail is the theme of this whole review. I propose measuring
the residuals A actually leaves, over the sweeps above in all three
precisions, and setting the tolerance at ten times the worst of them — a
number that comes from the algorithm rather than from a guess. If that turns
out to be within a small factor of today's, it stays as it is and the note
says so.

> **Q2. Tighten `ResidualTolerance` to a measured bound?** Recommended: yes,
> as a second commit after A, so that the two can be told apart if anything
> trips.

### 1.5 Tests (written first, failing today)

- No throw on the three stretched sweeps to l = 1000 (a thinned set in the
  default suite, the full one labelled slow), nor on the (l1, l1+l3, l3) grid.
- Agreement with `RacahSymbol` to 1e-11 of the table's scale wherever
  `RacahSumLength ≤ 2`, over the same sweeps. (Two, not three: three is where
  the oracle's own cancellation starts to show.)
- `double` against `long double` to 1e-10 on 300 seeded triples, l ≤ 500.
- `float` against `double` to 1e-3 on the same triples, and no throw.
- The existing cyclic-permutation tests at high degree, unchanged — they are
  the only check that sees a negated row.

---

## 2. Wigner tables (W2, and the rest of W3)

### 2.1 The limit is a formula, and it checks

The seed of each (m, θ) column, at l = |m|, is about (sin θ)^m, and the column
grows back to O(1) only once l sin θ ≳ m. A seed that underflows *and* matters
therefore needs m |ln sin θ| > |ln min| with m < lMax sin θ, that is
lMax · s |ln s| > |ln min|, and s |ln s| is largest at s = 1/e. So the plain
recursion is safe iff

    lMax  <  e · |ln min|.

Measured, as the worst defect of Σ_m |d^l_{nm}|² = 1 over five colatitudes
including asin(1/e), for |n| = 0, 1, 2, 4 (the four agree, so one column is
shown):

| `double`: e·\|ln min\| = **1926** | lMax 1700 | 1830 | 1900 | 1930 | 1960 | 2000 |
|---|---|---|---|---|---|---|
| defect | 3e-13 | 3e-13 | 2e-12 | 4e-9 | 4e-6 | 3e-3 |

| `float`: e·\|ln min\| = **237** | lMax 190 | 230 | 240 | 250 | 270 |
|---|---|---|---|---|---|
| defect | 3e-5 | 4e-5 | 3e-5 | 6e-5…3e-4 | 1e-2…3e-2 |

The onset is where the formula puts it, to within a per cent, in both
precisions and at every upper index tried. The small rise just *before* it
(1900 in double) is the seed going denormal and losing bits, which the formula
accounts for if the seed is required to stay above min/ε:

    MaxSafeDegree<Real>  =  ⌊ e · (|ln min| − |ln ε|) ⌋  =  1827 (double), 194 (float), 30 747 (long double).

At 1830 the double table is indistinguishable from one at 1700.

### 2.2 Proposal

**Enforce it.** A `constexpr`-evaluable `MaxSafeDegree<Real>()`, and a throw
above it from `Wigner`, `WignerMatrices` and `SphericalGrid::Impl` (the last
because a generating grid has no table to refuse). The message names the
precision and the limit. This turns a silent O(1) error into an exception,
which is the whole of the fix that D1 asked for; the extended-exponent
recursion that would lift the limit is not planned.

> **Q3. Enforce `MaxSafeDegree` = ⌊e (|ln min| − |ln ε|)⌋, i.e. 1827 in
> double?** Recommended: yes. The alternative is the bare 1926, which admits
> tables that have started to lose digits.

### 2.3 `float`

Enforced as it stands, a `float` grid stops at lMax = 194. Whether that is
acceptable is the question, and there are two honest answers.

**Option F1 — recursion in `double`, storage in `float`.** `float` earns its
keep in the *transform*: half the table traffic, half the field memory. It
earns nothing in the recursion. So for `Real = float` the recursion runs in
`double` into per-thread scratch and is narrowed on store. The limit becomes
the double one, 1827.

- `WignerMatrices` already computes each block into scratch and scatters it,
  so there it is a change of the scratch's type.
- `Wigner` writes in place today and would gain a scratch block per thread.
- A *generating* grid (`WignerValues::Generated`) runs the recursion inside
  every transform, so there the cost is paid per transform and not once.
  It must change too, and not only for the limit: the library promises, and
  tests, that the stored and generated paths agree **bit for bit**, which they
  can only do if both narrow from the same `double` values.
- Cost, measured: the double recursion takes **1.47×** the float one on one
  thread and **1.75×** on eight (lMax 128–512, nMax 2). For a stored or matrix
  grid that is construction time only — 213 ms against 123 at lMax = 256 on
  eight threads. For a generating `float` grid it is up to that factor on the
  recursion's share of every transform.
- `long double` and `double` are untouched, bit for bit.

**Option F2 — leave `float` in `float` and enforce 194.** No new code. `float`
is then an option for small problems only, and says so when asked for more.

The 3-j side was expected to need neither, proposal A leaving `float` 3-j
symbols correct to a few 1e-6 on the triples first tried. *Building it showed
otherwise for flat rows:* (0, 448, 448) in `float` is off by 5.7e-3 before and
after A, because a row of 900 values loses n² ε however good the algorithm.
So the 3-j rows are computed in at least double and narrowed as they are
scattered — the same principle as F1, and five lines, the row having always
been a scratch buffer.

> **Q4. F1 or F2?** Recommended: **F1**, since you asked for `float` to be a
> real option and 194 is not much of one — with the generating path included,
> because the bit-for-bit property is worth more than the generating `float`
> grid's speed. If you would rather not touch the generating path, F2 is the
> coherent alternative; a mixture (F1 for stored, float for generated) is the
> one thing I would not do.

### 2.4 The seed-row binomial

`Wigner.h` forms C(2|n|, ·) in `Real` before taking its root, which overflows
for |n| ≳ 520 in double and |n| ≳ 65 in float — beyond anything the library's
own use reaches (|n| ≤ 4), but the API admits nMax = lMax. The running product
is kept where it fits, because it is exact there and changing it would change
every existing table in the last bit; above that the entry is formed as
exp(½ ln C + (l−m) ln s + (l+m) ln c), with ln C accumulated as a sum of logs
(no `lgamma`, so ThreadSanitizer's reason for the present form still holds)
and the two logarithms `Arguments` already carries. Under F1 the `float` case
uses the double thresholds. *Test:* finite and unitary at n = 600, lMax = 700
in double, where two thirds of the table is non-finite today (87 062 entries
of 131 401; at n = 500 none are, at n = 520 nearly half).

### 2.5 Tests (written first, failing today)

- Above `MaxSafeDegree` the three constructors throw; at it they do not.
- Unitarity at `MaxSafeDegree<double>` for |n| ≤ 4 to 1e-11, at five
  colatitudes including asin(1/e). Slow — several seconds — so labelled.
- Under F1: a `float` table equals the narrowed `double` table exactly; a
  `float` grid at lMax = 512 round-trips through the loop kernel, the matrix
  kernel and the generating path to ~1e-4, and stored and generated agree bit
  for bit. None of these can be built correctly today.
- The large-|n| seed row.

---

## 3. Order of work

1. 3-j tests, then proposal A (Q1).
2. `ResidualTolerance`, measured (Q2).
3. `MaxSafeDegree`, its tests and the three throws (Q3).
4. The seed-row binomial.
5. F1 or F2 (Q4).

Each is a commit-sized step with the full gate run after it. Steps 1–2 and 3–5
touch different files and could be two separate stages if you would like to
commit them separately.
