# Wigner 3-j symbols: the plan

Written 2026-08-24, when `thoughts.md` §6 was picked up. That section is the
assessment — what is in `3j.h`, where its stability argument holds and where it
breaks, and the three routes out — and is not repeated. This is the work order
and the decisions.

Its own document rather than a section of `core-plan.md` or
`field-algebra-plan.md`, because it is independent of both: `3j.h` touches no
grid, no transform and no field, and `thoughts.md` says so in as many words.
Nothing here blocks or is blocked by anything there.

---

## 1. What was measured, which sharpens §6 rather than repeating it

`thoughts.md` §6 gives a table of the failure. It was reproduced against
current `main` so that this plan starts from numbers it has checked, and three
things came out sharper than that table states them.

**The identity is exactly one, over the plane the table already holds.**
Summing the squared symbols over the whole `(m₁, m₃)` grid at fixed degrees,

```
Σ_{m₁ m₃}  ( l₁ l₂ l₃ ; m₁ , −(m₁+m₃) , m₃ )²  =  1
```

which is the standard `Σ_{m₁m₂}(·)² = 1/(2l₃+1)` at fixed `m₃`, summed over
the `2l₃+1` values of `m₃`. Measured at `(l,l,l)` for `l = 2, 5, 10, 20`: 1 to
the last digit each time. **It needs no reference implementation, it is one
loop over a table that has just been built, and it is the whole basis of what
follows.**

**The boundary is where §6 puts it, and now has a sharp edge.** Taking
`l₃ = l₁ + l₂` and reporting `(2l₃+1) × Σ`, which should be `2l₃+1` exactly:

| `l` in `(l, l, 2l)` | expected | double | long double |
|---:|---:|---:|---:|
| 10 | 41 | 41 | 41 |
| 20 | 81 | 81 | 81 |
| 25 | 101 | 101 | 101 |
| 30 | 121 | **121.49** | 121 |
| 35 | 141 | `3.7 × 10⁶` | **141.02** |
| 40 | 161 | `3.3 × 10¹¹` | `1.3 × 10⁵` |
| 60 | 241 | `5.7 × 10³⁴` | `1.2 × 10²⁸` |

So double is exact to `l = 25`, first departs at 30, and is catastrophic by
35; long double buys **five to ten degrees** and then fails identically. That
is the signature §6 identified — a recursion run in its unstable direction,
which extra precision delays and does not cure.

*Refined while writing T1's tests, and it matters for where a tolerance
goes.* At the resolution of that table the failure looks like a cliff between
25 and 30. It is not: the departure from one is `1e-16` out to `l = 16`,
`4.4e-14` at 20, **`1.8e-8` at 25**, `1.2e-5` at 28 and `4.1e-3` at 30. So the
loss is **exponential from about `l = 20`**, and there is a band around 22 to
28 where the values are degrading but still usable. A check has to decide
what to do about that band, which is what [J1]'s tolerance is for.

Away from stretched the boundary is much further out. At `l₃ = 3l₁/2` the sum
is exact to `l = 40`, departs at 60 (181.06 against 181) and is catastrophic
by 80. At `(l, l, l)` it holds to `l = 20` and, per §6, to 128.

**And the closed-form seed is exact, not merely accurate.**
`(j j 0 ; m −m 0) = (−1)^{j−m}/√(2j+1)` reproduces to **zero** absolute
difference at `j = 1, 4, 30`. Worth knowing because it means a test on that
family discriminates a wrong phase or a wrong index and nothing else — there
is no tolerance in which a bug could hide.

**One more fact, and it changes an item on §4's list.** `examples/wigner3j.hpp`
is **not** an independent implementation. It is the earlier standalone port of
the same Woodhouse `wig2` routine that `GSHTrans/src/3j.h` was derived from —
same recursion, same convention, same `(-1)^{m₁}` phase — and its own header
comment names the same accuracy caveat and even proposes the same completeness
relation as "a cheap self-diagnostic". So the repository does not hold two
Wigner 3-j codes worth keeping; it holds one code twice. That removes the
possibility of using one as an oracle for the other, which was the attractive
reading of §4's note, and it settles what to do with it.

---

## 2. Decisions taken

**[J1] The completeness relation becomes a runtime self-check, not only a
test.** `Wigner3jMatrix` sums its own squares at construction and throws if the
result differs from one by more than a tolerance.

This is the decision the rest of the plan is built on, and the argument for it
is that **the failure mode is silent and catastrophic**. A library that returns
`10¹¹²` where it should return a number of order `0.1` is worse than one that
refuses, because the caller has no way to know: the values are finite, the
table has the right shape, and a coupling sum built on them produces a
plausible-looking wrong answer. §6 found this with one line, and the same line
can stand between every caller and the same surprise.

*It is affordable, which is why it is possible at all.* The table is
`(2l₁+1)(2l₃+1)` entries and the check is one pass over it — the same order as
building it, and a small multiple less work. Nothing else in this library gets
a check that strong for that price; the reason it can here is that the identity
is exact and needs no reference.

*It detects rather than locates*, which is enough: [J3] recomputes the whole
table when it fires, so knowing which entry went wrong would buy nothing.

*The tolerance is a decision in itself and is deliberately loose.* The point is
to catch `10⁶` and above, not to police the last bits — §6's own measurements
put the honest floor at `10⁻¹⁴` for fat triangles at `l = 128`, and a tight
tolerance would reject good tables. Something near `10⁻⁸` separates the two
populations by twenty orders of magnitude, so the choice is not delicate.

**[J2] The duplicate goes, and it is a deletion rather than a promotion.**
`examples/wigner3j.hpp` and `examples/wigner3j_tests.cpp` are the same
algorithm as the library's, so keeping them cannot buy the independent-oracle
check that `Interpolation::CubicSpline` bought `SplineDerivative` (§19.5) or
that the loop kernel buys the matrix kernel ([C12]). What they would buy is two
copies of one recursion drifting apart. `examples/wig.cpp` moves onto the
library's own type, which is what an example should have been doing.

*This is the answer to `thoughts.md` §4's "there are two Wigner 3-j codes in
the repository and neither is exercised by the suite".* After [J1] and step
T1 there is one, and it is.

**[J3] The route is the hybrid, and the self-check is what dispatches it.**
§6 offers three: fix the recursion's direction (Schulten–Gordon), compute
exactly in integer arithmetic (`wigxjpf`), or hybridise the existing recursion
with Racah's closed form, and it recommends investigating the third "before
committing to either of the others". Taken, and the reason is stronger than §6
had:

Racah's formula is a single alternating sum whose length is
`min(l₁+l₂−l₃, …) + 1`. **At `l₃ = l₁ + l₂` that sum has exactly one term**,
so there is no cancellation at all and the only error is the rounding of a
product of factorials — and one term is the case where the recursion is at its
worst. Going the other way, the fat triangles that destroy Racah's long
alternating sum are exactly where the recursion is at the noise floor. The two
methods fail in complementary regimes, which is the observation §6 makes.

What §6 does not have is the dispatcher. **[J1] supplies it.** The rule is not
a boundary formula in `(l₁, l₂, l₃)` that someone has to derive, calibrate and
then re-derive for `long double`:

> Build the table by the existing recursion. Check the identity. If it fails,
> rebuild by Racah and check again.

That is a dispatch on *the thing that actually went wrong*, it is free of
tuning constants, it adapts to precision without being told, and it degrades
correctly on a triple neither method handles — because the second check fires
too, and the library throws rather than returning nonsense.

*The known cost, so that it is not discovered:* Racah's factorials must be
formed in log space and exponentiated, and `core-plan.md` T11 records what
that does — a logarithm of size `O(l log 4)` loses bits in proportion, giving
about `1000 ε` at the sizes this library works at. For a one-term sum that is
the whole error, so near-stretched symbols would carry a relative error near
`10⁻¹³`. Against `10¹¹²`, that is the trade being made, and it should be stated
in those terms rather than as "exact".

*And what would change the decision:* if the measurement of step T2 shows a gap
— triples where the recursion fails and Racah's sum is already too long to be
trusted — then the hybrid does not cover the space and Schulten–Gordon becomes
the answer rather than the fallback. T2 exists to find out, and it is an
afternoon.

**[J4] `wigxjpf` is rejected for now, and the condition for revisiting it is
6-j.** *T2 adds a second condition, and it is now the stronger of the two: the
gap.* The classical pair does not cover intermediate shapes above `l ≈ 80` or
fat triangles above `l ≈ 160`, and closing that means either Schulten–Gordon
or exact integer arithmetic. If the gap ever needs closing, the choice between
those two should be made together with 6-j rather than separately. It is C rather than C++, a few thousand lines, and not header-only,
which cuts against the preference `thoughts.md` §5 states and which this
library has honoured everywhere except BLAS — and BLAS earned its exception by
being a thing every target machine already has, which `wigxjpf` is not.

What would earn the exception is wanting 6-j and 9-j as well, since the exact
route delivers all three and the recursion route means implementing
Schulten–Gordon again for each. So the decision is made **with 6-j in view**,
per §6, and the position is: nothing needs 6-j today, the hybrid does not
foreclose it, and if 6-j is ever wanted the right move is to reconsider
`wigxjpf` for the whole family rather than to hand-write a second recursion.

**[J5] The stack keeps looping, and the `l₂` recursion is not built.**
`Wigner3jStack(l₁, l₃)` builds one `Wigner3jMatrix` per middle degree and does
not use the three-term recursion in `l₂` at all. That is `O(l³)` work where the
recursion would be `O(l²)`, and it is left alone: the `l₂` recursion has its
own stability direction to get right, it would need its own self-check, and
nothing in this library consumes a stack yet. Recorded so that the loop reads
as a decision rather than as an oversight, and so that whoever wants coupling
sums at scale knows where the factor of `l` is.

---

## 3. The steps

**T1 — the test family, and the boundary as a known quantity.** *Done*, as
`tests/TestThreeJ.cpp`. Before any algorithm changes. Four groups, none
needing a reference implementation:

- **the identity**, swept over the triangle space — fat, intermediate and
  stretched, at several degrees — asserting `Σ = 1` where §1's table says it
  should hold, and asserting that it *fails* where §1 says it fails. The
  second half matters as much as the first: it is what stops a later change
  silently narrowing the working range, and it is the negative test that makes
  the positive one meaningful;
- **the closed forms**: `(j j 0; m −m 0) = (−1)^{j−m}/√(2j+1)`, exact to zero
  per §1, and the stretched form `l₃ = l₁ + l₂`, which is Racah's one-term
  case and therefore also closed;
- **exact comparison at small degrees**, against values computed in rational
  arithmetic and written into the test as literals, which pins the convention
  and the phase rather than the accuracy;
- **the symmetries**: the table's own reflection, and invariance under an even
  permutation of the columns.

*This is the step `thoughts.md` says goes in first*, and it stands on its own:
after it, `3j.h` has coverage where it had none, and the boundary is a
documented property rather than something someone rediscovers.

*Two things this document had wrong, both found by writing the tests against
it.* The gradual decay above, and the closed form at the stretched corner:
§3's sketch of it as `(−1)^{l1−l2} sqrt((2l1)!(2l2)!/(2l1+2l2+1)!)` is the
formula for a different corner. At the **fully** stretched symbol every
factorial cancels and it is simply

```
(l1  l2  l1+l2 ; l1  l2  −(l1+l2))  =  1 / sqrt(2 l3 + 1)
```

which the library reproduces to zero absolute difference. That is a much
stronger test than the one intended, because there is no arithmetic in the
expected value to be wrong in the same way the code is — and it is the case
T4's Racah path must reproduce, since the sum there has exactly one term.

**T2 — price the hybrid, per [J3].** An afternoon, and it is a measurement
rather than an implementation. Racah in log space, written for the experiment
and not for keeps; compare against the recursion across the triangle space and
across degrees; and answer the one question that decides the route — **is there
a gap?** Plot, or tabulate, the region where the recursion's identity fails
against the region where Racah's alternating sum has lost too much, and see
whether they overlap or leave a band uncovered.

If they overlap, [J3] stands and T3 is small. If they do not, the answer is
Schulten–Gordon and this plan's §2 is amended rather than followed.

#### T2's answer: they do not overlap, and [J3] is amended

*Done.* Racah first agrees with the recursion where the recursion is good —
`1.3e-15` at `(2,2,2)`, `2.6e-13` at `(12,12,12)`, `1.6e-11` at `(20,20,20)`,
degrading with the alternating sum's length exactly as expected — so the two
are implementations of the same thing and the comparison below is a
comparison of accuracy rather than of convention.

Then the question. Completeness departure from one, by each method, with the
longest alternating sum Racah forms for that triple:

| triple | recursion | Racah | sum length | covered by |
|---|---:|---:|---:|---|
| (35,35,70) | 2.6e4 | **1.000000** | 1 | Racah |
| (64,64,128) | 8.0e35 | **1.000000** | 1 | Racah |
| (128,128,250) | 2.1e102 | **1.000000** | 7 | Racah |
| (256,256,500) | 1.1e241 | **1.000009** | 13 | Racah |
| (60,60,90) | 1.000312 | **1.000002** | 31 | both, marginally |
| (70,70,105) | 8.4e2 | **1.000103** | 36 | Racah, marginally |
| (80,80,120) | 1.7e7 | 1.146 | 41 | **neither** |
| (90,90,135) | 5.9e9 | 3.8e2 | 46 | **neither** |
| (100,100,150) | 4.5e18 | 1.6e6 | 51 | **neither** |
| (128,128,192) | 2.1e34 | 2.1e16 | 65 | **neither** |
| (128,128,128) | **1.000000** | 5.8e18 | 129 | recursion |
| (160,160,160) | 1.056 | 7.5e30 | 161 | **neither** |
| (200,200,200) | 1.4e10 | 1.3e46 | 201 | **neither** |

**There is a gap, and it is not small.** Two regions fall outside both
methods: **intermediate shapes** — around `l₃ ≈ 1.5 l₁` — from about
`l = 80`, and **fat triangles** from about `l = 160`. The second is news
beyond §6, which had `(l,l,l)` holding to 128: it does, and 160 is already
5.6% wrong.

*So [J3]'s hope is refuted and its own contingency applies* — "if they do
not, the answer is Schulten–Gordon and this plan's §2 is amended rather than
followed". §2 is amended below rather than abandoned, because the measurement
also shows the hybrid is worth building anyway.

**T3 — [J1]'s self-check.** *Done.* The identity in `Wigner3jMatrix`'s
constructor, with the tolerance and the throw. It lands before T4
deliberately: a caller who picks the library up between the two gets a refusal
where they used to get `10¹¹²`, which is a strict improvement even with no
second algorithm behind it.

*The test is that it fires* — a stretched triple at `l = 40` must throw — *and
that it does not fire* on the whole of the range T1 established as good, which
is what says the tolerance is not doing damage. Both are asserted, and the
second names every family T1 measured as accurate, plus a triple that is not
a triangle at all: its table is identically zero by design, so its
completeness sum is zero and correctly so, and the check must skip it rather
than refuse it.

*The tolerance is `100 · sqrt(ε)`*, so it follows the precision rather than
assuming double — `1.5e-6` there, `3.4e-2` in single, `3.3e-8` in long
double. That accepts the degrading band around `l = 22` to 25 and refuses
from about 28, which is the line the refined measurement above puts it on.
The populations it separates are twenty orders apart, so the choice is not
delicate; where it sits inside the degrading band is the only judgement in
it.

*What T1's boundary tests became.* They asserted that the completeness sum
**is not one** at stretched high degree. They now assert that construction
**throws**, which is the same fact one layer up, and they carry a note saying
that T4 should invert them.

**T4 — Racah, and the dispatch.** *Done, and built despite T2*, because what
T2 refutes is [J3]'s claim that the hybrid *covers the space*, not its claim
that the hybrid is worth having. The union is a large gain over the recursion
alone: everything near-stretched — which §6 rightly says "is not an exotic
corner but the top of every coupling sum" — becomes exact where it used to be
`10¹⁰²` wrong. What T2 changes is the honest description of the result, from
"covered" to "covered except a named band", and the check is what stands
between a caller and that band.

The dispatch is [J3]'s and it works as designed: recurse, check, fall back to
Racah, check again, refuse. No boundary formula, no tuning constants, adapts
to precision without being told.

*The test that earns it* is the one T1 wrote as a negative, inverted: the
stretched triples that used to be refused now construct, satisfy the identity
to `1e-12`, and reproduce `1/sqrt(2l₃+1)` at the corner to `1e-13`. A second
test asserts the gap — that `(90,90,135)`, `(100,100,150)`, `(128,128,192)`,
`(128,128,160)` and `(200,200,200)` are **refused**, since being answered
wrongly there is the thing this whole plan exists to prevent.

#### [J1]'s tolerance was wrong, and the measurement is why

That decision says the two populations are "twenty orders of magnitude apart,
so the choice is not delicate". **They are not cleanly separated.** In double
the departures fall into three groups, not two:

```
good      1e-16 to 1e-8   fat triangles at moderate degree, and anything the
                          closed form answers
marginal  2e-6 to 1e-4    (60,60,90), (256,256,500), (70,70,105)
lost      1.1  to 1e112   everything past the boundary
```

A tolerance set just above the noise — the first attempt used
`100·sqrt(ε)`, about `1.5e-6` — **refuses the marginal band**, whose values
carry perhaps `5e-5` relative error and are worth having. The real gap is
between `1e-4` and `1.1`, so the tolerance goes there instead: `1e-3`, with
the epsilon-scaled floor kept for single precision, which cannot resolve the
marginal band at all. Three orders of margin below the lost population is
ample.

**T5 — [J2]'s deletion.** *Done*, and there were **three** copies of the
Woodhouse routine in this repository rather than the two §4 counted:
`GSHTrans/src/3j.h`, `examples/wigner3j.hpp`, and a third written out inline
inside `examples/wig.cpp`. All three standalone files are gone, replaced by
`examples/22-wigner-3j.cpp` — a numbered example against the library's own
types, registered as a test like every other, which is what an example should
have been doing. It ends by demonstrating the refusal, since that is the
thing a caller most needs to know about.

---

## 4. What this does not touch

The grid, the transform, the Wigner `d`-functions and the field algebra.
`3j.h` shares no code with `Wigner.h` despite the name — one is rotation-matrix
elements on a colatitude grid, the other is coupling coefficients of three
degrees — and nothing in this plan brings them together.

**And what it deliberately does not build:** the Gaunt integrals and coupling
coefficients that are the reason to want 3-j at all. Those are a consumer, they
belong wherever the physics does, and building them now would fix an interface
before there is a caller to fix it for. What T1 to T4 deliver is a 3-j
implementation that is either right or says it is not, which is the
precondition for any of it.
