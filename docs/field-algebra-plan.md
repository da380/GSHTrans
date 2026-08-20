# Field algebra: unified plan

Single planning document for the rebuild of the field-algebra layer. It
supersedes the earlier `field-algebra-plan.md` and
`spin-component-expression-plan.md`, whose content is absorbed here.

Companion: `canonical-components.tex`, which fixes the conventions and states
the facts the type system is meant to enforce. That note is the authority on
mathematics; this one is the authority on design.

Companion: `core-plan.md`, the **step −1** plan for the numerical core. Three
of the decisions below (Q3, Q4, Q5) turned out to be changes to the core rather
than to the field layer, so they were moved there; §1 records which core items
must land before which phase-1 steps.

The **[Qn]** questions have been answered. §9 records each decision, its
rationale, and what it changed in the body — including the two places where the
answer went against the recommendation. The body text above §9 is written as
settled design. §10 recorded two items as still open; both are now closed, and
§10 says how. **Nothing in this document is waiting on an answer.** The work
order is `core-plan.md` §8.

---

## 1. Scope and roadmap

Rebuild the field-algebra layer from scratch. The numerical core — `Wigner`,
`GaussLegendreGrid`, `Indexing`, `Views`, `3j` — is not rewritten, but it is no
longer treated as frozen: the decisions of §9 imply a bounded set of changes to
it, planned separately as **step −1** in `core-plan.md`. Nothing in phase 1
touches the core's algorithms; step −1 touches its interface, its ownership
model, and one grid size.

Superseded: `CanonicalComponentField*`, `CanonicalComponentExpansion*`, and the
dormant `ScalarField/`, `VectorField/`, `MatrixField/`, `MatrixFieldOld/`,
`ScalarFieldExpansion/` trees.

*Status.* All of those are now **deleted**, together with `FieldBase.h` and
`ExpansionBase.h` — which settles the `CheckPointIndices` question step 7
raises below: the buggy comparison was deleted rather than fixed.

`CanonicalComponentExpansion*` was kept back at step 7, on the grounds that
deleting the spectral side would remove working, tested code with nothing to
put in its place until phase 5. **It went after step F**, because that
argument stopped holding in two ways that step −1 created:

- **It could still represent what Q5 removed from the library.**
  `Traits<CanonicalComponentExpansion>` set `MRange = NonNegative` whenever
  `Value = RealValued`, at *any* `_N`, and its test instantiated
  `RealCanonicalComponentExpansion<1, Grid>`. That is the reduced `m ≥ 0`
  storage scheme core step A deleted, surviving on the spectral side of a
  library whose spatial side is now type-prevented from expressing it (§9,
  Q5).
- **It held the grid as `_Grid&`**, taken from a non-const reference, which is
  the ownership model core step B removed everywhere else (§3.7). Its test
  pinned those semantics with an address comparison, so the tree could not be
  moved onto the handle without rewriting the test that justified keeping it.

Deleting it cost no coverage that phase 5 will want: it touched no transform,
and `GSHIndices` is exercised far more heavily by `TestGaussLegendreGrid`.
`Traits.h` went with it, the CRTP `Internal::Traits` scaffolding having had no
other user, and the `GSHTrans/Expansion` umbrella went with that. Phase 5
therefore starts from a clean spectral side rather than from a superseded
one — which is what it was always going to do.

| phase | content |
|---|---|
| **−1** | **Numerical core: ownership, transform interface, grid sizing, plan caching, batching. Separate document, `core-plan.md`. Its steps A–D precede phase 1.** |
| **1** | **Spin fields — a single field of definite upper index, and its lazy algebra. This document's main subject.** |
| 2 | Tensor storage: multi-index machinery, orbit enumeration, `TensorField` over one buffer handing out phase-1 views |
| 3 | Tensor algebra: product, metric contraction, trace, transpose, (anti)symmetrisation, rank-4 on rank-2 |
| 4 | Reality reduction: store one component per negation orbit, derive the rest on the fly |
| 5 | Spectral side: tensor-valued expansions, batched transforms, raising/lowering |

Ranks 0, 1, 2 and 4 are exposed. The index machinery is generic in the rank
because that is cheaper than writing each rank by hand, but the API surface
stays bounded.

Phases 2–5 appear here only in §8, and only to the extent that they constrain
phase 1.

### Why a rewrite rather than repair

Three defects in the current expression layer are structural. All three are
confirmed against the code:

1. **Operand lifetime and slicing.** Nodes hold
   `const CanonicalComponentFieldBase<N, Derived>&` — a *base* reference
   (`CanonicalComponentFieldUnary.h:150`, `:199`;
   `CanonicalComponentFieldBinary.h` likewise). Nested expressions bound to
   `auto` dangle; confirmed under AddressSanitizer (`stack-use-after-scope`,
   `CanonicalComponentFieldAdd::Grid()`). The rvalue overloads do not help —
   `operator+(Base&& u0, Base&& u1) { return u0 + u1; }` rebinds to an lvalue
   in the body and dispatches to the `const&` form.
2. **Conjugation does not reverse the upper index.**
   `CanonicalComponentFieldConj` derives from `...Base<_N, ...>`
   (`CanonicalComponentFieldUnary.h:107`). It must be `-N` (theory note §5).
3. **`real` and `imag` are offered at every upper index.**
   `CanonicalComponentFieldOverloads.h:46–62` are unconstrained in `N`; the
   `requires(N == 0)` clauses only begin at `:68`. They are covariant only at
   `N = 0` (theory note §5).

Fixing 2 and 3 changes the signatures of the affected nodes, and fixing 1
changes how every node stores its operands, so little of the present code
survives.

---

## 2. Conventions: what is blocked and what is not

The theory note flags three convention-dependent points: the sign in
`e_± = ∓(θ̂ ± iφ̂)/√2`, the index order and value convention of the Wigner
`d`-functions, and the overall signs in the ð relations.

**The index algebra is convention-free.** Nothing in phase 1's types depends
on any of the three. What depends on them is the *reality phase* `(-1)^N`
(phase 4), the spectral relations (phase 5), and any phase-1 test whose oracle
runs through a transform.

**The `d`-function convention is now settled and verified.** It was the one
item gating step 6 of §7 and the transform-based oracles of test family 3;
that gate is lifted. The statement of record:

> `Wigner`'s stored value at upper index `N`, degree `l`, order `m` is
> `sqrt((2l+1)/(4π)) · d^l_{Nm}(θ)`, where `d^l_{Nm} = P^N_{lm}(cos θ)` is the
> generalised Legendre function of **Dahlen & Tromp (1998) eq. (C.115)** —
> **upper index first**, consistent with the Q12 choice of orbit
> representatives.

Verified against the code, not inferred from it. The `l = 1` values of
D&T (C.115),

```
P^{-1}_{1,-1} =  (1 + cos θ)/2      P^{0}_{1,-1} =  sin θ/√2      P^{1}_{1,-1} = (1 - cos θ)/2
P^{-1}_{1, 0} = -sin θ/√2           P^{0}_{1, 0} =  cos θ         P^{1}_{1, 0} = sin θ/√2
P^{-1}_{1, 1} =  (1 - cos θ)/2      P^{0}_{1, 1} = -sin θ/√2      P^{1}_{1, 1} = (1 + cos θ)/2
```

reproduce all nine stored values at `θ = 0.7` after dividing by
`sqrt(3/(4π))`. The check discriminates between the two candidate conventions:
the four entries with `N − m` odd change sign under transposition, since
`d^l_{mN} = (−1)^{N−m} d^l_{Nm}`, and the code matches `d^l_{Nm}`. Two existing facts corroborate it —
the recursion seed `WignerMinOrder(l, n)` evaluates to `d^l_{n,−l}`, and
`tests/CheckLegendre.h` already pins the `N = 0` row against
`std::sph_legendre`, which fixes the orthonormal scaling and the
Condon–Shortley phase.

The earlier text in this section read the convention off the *indexing* as
`d^l_{mN}` and flagged the value convention as unchecked. That reading was
wrong in the index order and right to flag the values; the table above settles
both. Because the fact is load-bearing for phase 4's stored components, it
becomes a permanent test rather than a note — `core-plan.md` §8 task T1.

The `Ortho` normalisation of `eq:gsh` is separately pinned by the one-point
grid: `GaussLegendreGrid.h:141` sets `out[0] = in[0] * 2 / inv_sqrtpi`, i.e.
`f⁰₀₀ = 2√π f`, exactly as the note derives. After `core-plan.md` step A2 it
is the *only* normalisation the library offers.

Still convention-dependent and still unchecked, because nothing before phase 4
observes them: the sign in `e_± = ∓(θ̂ ± iφ̂)/√2` and the overall signs in the
ð relations. Neither blocks any step of §7.

## 3. Phase 1: the object

The unit is a single field of definite upper index `N` on the sphere. Tensors
come later and are built from it.

### 3.1 Naming

The level-0 object is renamed. "Canonical component" is reserved for the tensor
layer, where a component is identified by a *multi-index*, not by an upper
index — theory note §2 is explicit that for rank ≥ 2 the two are different
things and that a collection labelled only by `N` does not determine the
tensor. Continuing to call the level-0 object a canonical component would
entrench that confusion.

Settled: `SpinField<N, Grid, Value>` for the owning terminal, concept
`SpinWeighted` for any node, `SpinFieldExpr` informally. Parameter order follows
the existing `CanonicalComponentField<N, Grid, Value>` for continuity.

### 3.2 The node concept

One concept, `SpinWeighted<T>`, defines what every node — terminal, view or
expression — provides:

- `T::UpperIndex` — `static constexpr Int`, the spin weight `N`.
- `T::Value` — `RealValued` or `ComplexValued` (the existing tags in
  `Concepts.h`), **subject to `Value == RealValued ⟹ N == 0`**. This is theory
  note §7 item 5 promoted into the concept: it is a constraint on
  `SpinWeighted` itself, so no node — terminal, view or expression — can claim
  to be real-valued at nonzero upper index. §3.4 shows the constraint is closed
  under the whole algebra, so it costs no expressiveness; what it forbids is a
  real-valued *terminal* at `N ≠ 0`, and that is a thing the library should not
  be able to represent.
- `T::Real`, `T::Complex`, `T::Scalar` — `Scalar` is `Real` when
  `Value` is `RealValued`, else `Complex`.
- `T::GridType`, and `const GridType& Grid() const`. The grid is a
  value-semantic handle after step −1 (§3.7, `core-plan.md` step B): copying it
  is a pointer copy, and two grids are "the same grid" when
  `Grid().Identity() == other.Grid().Identity()`. Terminals hold one by value;
  expression nodes forward their left operand's.
- `Scalar operator[](Int iTheta, Int iPhi) const` — **by value, always**,
  including on terminals. Uniform value return keeps expression, view and
  terminal interchangeable; a mutable reference accessor exists on terminals
  and mutable views only, outside the concept.
- `template <typename S> requires std::convertible_to<Scalar, S>`
  `void EvaluateInto(std::span<S>) const` — writes the field in the canonical
  layout. Default implementation loops over `operator[]`; terminals override
  with a contiguous copy. Templated on the output scalar, not fixed at
  `span<Scalar>`, so that a `RealValued` expression can be written into a
  `Complex` destination (§3.8); the `requires` clause makes the reverse a
  compile error rather than a truncation.

No virtual functions. No CRTP base is *required*; a small CRTP convenience base
supplying derived helpers (point-index iteration, `Size()`, `Integrate()`
forwarding) is acceptable if it stays free of data members and of the lifetime
bugs the old `FieldBase`/`...FieldBase` pair introduced. Everything dispatches
through the concept.

**Iteration order contract.** `(iTheta, iPhi)` with phi fastest, flat index
`iTheta * nPhi + iPhi`. Verified as the transform's layout:
`GaussLegendreGrid.h:169` computes `offset = iTheta * nPhi`, and
`GridBase::Points()` / `PointIndices()` are `cartesian_product(colatitudes,
longitudes)`, which is theta-major. State this once, on the concept, and make
`EvaluateInto`'s default the definition of record.

**Thread-safety contract.** Evaluation is `const` and stateless: `operator[]`
and `EvaluateInto` mutate nothing, and callables are invoked as `const`.
Concurrent evaluation of the same node from several threads into disjoint
outputs is therefore safe by construction. This is a documented guarantee, not
an accident — the layered layer of §8 parallelises over slices on it.

**Slice targets.** The span passed to `EvaluateInto` may be a slice of a larger
buffer (a radial slab, a tensor-component plane), contiguous in the canonical
order, with **no alignment requirement beyond `Scalar`'s**. That contract is
made honest by step −1 (`core-plan.md` step C), which stops the transform
executing FFTW plans on caller storage; until it lands, the contract is a
promise the core does not yet keep.

**Views satisfy the concept.** A non-owning view — a pointer into someone
else's buffer plus a grid handle — must be admissible as a `SpinWeighted` node.
This is why `operator[]` returns by value everywhere and why the grid handle is
shared: the tensor-component views of phase 2 and the radial-slice views of §8
depend on it. Phase 1 ships a minimal view fixture to pin it down (test family
4). Note that the existing `Views.h` is spectral-side only (`GSHView`,
`GSHSubView` over coefficient storage); the spatial view is new work.

### 3.3 Operand storage and value categories

This is the load-bearing decision; defect 1 lives here.

Rule, applied per operand at operator-call time via forwarding references:

| operand kind | value category | stored as |
|---|---|---|
| terminal or view | lvalue | `const T&` |
| terminal or view | rvalue | `T` (moved in) |
| expression node | lvalue | `T` (copied) |
| expression node | rvalue | `T` (moved) |

Mechanism: `operator+(L&& l, R&& r)` deduces `L`, `R` as forwarding references;
the node's template arguments are the *deduced* (possibly reference) types; the
stored member type is

```cpp
template <typename A>  // A as deduced; may be an lvalue reference
using OperandStorage =
    std::conditional_t<IsTerminal<std::remove_cvref_t<A>>::value &&
                           std::is_lvalue_reference_v<A>,
                       const std::remove_cvref_t<A>&,
                       std::remove_cvref_t<A>>;
```

Consequences to record explicitly:

- `auto e = MakeField(...) + v;` is safe: the rvalue terminal is moved into the
  node. This closes the hole a plain `IsTerminal` trait leaves open (an rvalue
  temporary held by `const&`).
- Expression nodes are always held by value. They are small (a few references
  and a functor), so copying is cheap, and it makes `auto f = e * w;` safe when
  `e` is a named expression that later goes out of scope before `f`.
- Residual, undetectable-in-C++ hazard: an *lvalue terminal* destroyed while an
  expression referencing it is alive. Same residual hazard as Eigen. Document
  it in one place; do not attempt `shared_ptr` ownership of terminals to paper
  over it — that changes the cost model of every field.
- Because operands are stored as `const T&` or `T`, never as a base reference,
  every node is copyable and `auto` never slices. Add
  `static_assert(std::copy_constructible<Node>)` in the node templates.

Callables (§3.6) are always decayed and stored by value.

`IsTerminal` is a trait specialised by the owning field and the view types, not
a property inferred from the interface.

### 3.4 The index algebra as traits

Two node templates only:

```
Binary<Op, IndexRule, L, R>   // L, R as deduced; see §3.3
Unary<Op, IndexRule, A>
```

`IndexRule` is a tag carrying the compile-time computation and the
admissibility constraint:

| rule | result `N` | constraint | used by |
|---|---|---|---|
| `Equal` | `N_L` | `N_L == N_R` | `+`, `-` |
| `Sum` | `N_L + N_R` | — | `*` |
| `FirstOnly` | `N_L` | `N_R == 0` | `/` |
| `Same` | `N_A` | — | unary `-`, scalar `*`, `/` |
| `Negate` | `-N_A` | — | `conj` |
| `Zero` | `0` | per-op | `abs`, `abs2`, `real`, `imag`, `Map` |

Constraints are expressed as `requires` clauses **on the free operators**, so
that an unlawful combination is a clean overload-resolution failure at the call
site rather than an error inside a node, and negative tests can be written as
`static_assert(!requires { u + v; });`. A `static_assert` with a readable
message inside each node duplicates the constraint as a backstop for anyone
constructing nodes directly.

**Value propagation.** `CombinedValue<L, R>`: `RealValued` iff both operands
are `RealValued`. Scalar multiplication by a `Complex` scalar promotes
`RealValued` to `ComplexValued`; by a `Real` scalar it preserves. `conj` and
unary `-` preserve `Value`. `real`, `imag`, `abs`, `abs2` produce `RealValued`.

**Closure lemma for the reality constraint.** The concept's
`RealValued ⟹ N == 0` (§3.2) is preserved by every node in §3.5, so it never
has to be re-checked at a node and never blocks a lawful expression:

| rule | argument |
|---|---|
| `Equal` (`+`, `-`) | result real iff both operands real ⟹ both `N = 0` ⟹ result `N = 0` |
| `Sum` (`*`) | result real iff both real ⟹ `N_L = N_R = 0` ⟹ `N_L + N_R = 0` |
| `FirstOnly` (`/`) | result real iff both real ⟹ `N_L = 0` |
| `Same` (unary `-`, real scalar `*`) | `N` and `Value` both unchanged |
| `Same` with a complex scalar | promotes to `ComplexValued`; premise discharged |
| `Negate` (`conj`) | real ⟹ `N = 0` ⟹ `-N = 0` |
| `Zero` (`abs`, `abs2`, `real`, `imag`, `Map`) | result `N = 0` unconditionally |

So the only construction the constraint rejects is a `RealValued` terminal or
view at `N ≠ 0`. Write the table's negative cases as `static_assert`s in test
family 1; the lemma is what makes them exhaustive.

**Precision.** One `Real` per expression tree: constrain
`std::same_as<typename L::Real, typename R::Real>`. No mixed-precision
promotion.

### 3.5 Node inventory

All pointwise **and index-preserving**: node `(iTheta, iPhi)` reads its
operands only at `(iTheta, iPhi)`. This is an invariant, not an accident. It is
what makes the aliasing argument of §3.8 work, and it is the property that
gradients, raising/lowering and radial derivatives will *not* have, which is
why they live in the spectral layer instead. The stronger "index-preserving"
form matters because a re-indexing view (a phi shift, a transpose) would be
pointwise in the loose sense and would break the argument.

Binary, via `Binary` + functor + rule:

- `f + g`, `f - g` — `Equal`.
- `f * g` — `Sum`. Note in the header that the product of two band-limited
  fields exceeds the grid's truncation; the type system permits it (it must),
  and the dealiasing question belongs to the grid (§8).
- `f / g` — `FirstOnly`. Runtime hazard (zeros of `g`) is the user's.

Unary, via `Unary` + functor + rule:

- unary `-` — `Same`.
- `s * f`, `f * s`, `f / s` for scalar `s` — `Same`, scalar captured by value.
  `s` convertible to `Real` (preserving) or `Complex` (promoting).
- `s / f` — dedicated unary with constraint `N_A == 0`.
- `conj(f)` — `Negate`. At `N == 0` on a `RealValued` operand it is the
  identity; permitted, not special-cased.
- `abs(f)` — `Zero`, any `N`, result `RealValued`.
- `abs2(f)` — `Zero`, any `N`, result `RealValued`. Theory note §5 lists
  `|f|²` as admissible in its own right, and it is what `Norm` actually needs;
  routing that through `real(conj(f) * f)` costs a complex multiply per point
  and a `Complex` accumulator for a quantity known to be real.
- `real(f)`, `imag(f)` — `Zero`, constraint `N == 0`.
- `Map(f, F)` for arbitrary callables — §3.6. `Map`, not `Transform`, which
  would collide with `ForwardTransformation`/`InverseTransformation` and with
  "spherical harmonic transform" throughout the codebase's vocabulary.

Document (do not encode) that `|f|` is not band-limited and not smooth at zeros
of `f`, so the spectral layer must not assume that an `N == 0` expression is
truncatable at the grid's `lMax`.

Nothing else. In particular no `pow`, `exp`, `log` as named nodes — they are
all `Map(f, F)` at `N == 0`, and adding named sugar later is trivial.

### 3.6 Callable nodes

`Map(f, F)` with constraint `N == 0`:

- `F` decayed, stored by value, moved in; `std::invocable<F, Scalar>` required.
- Result `Value`: `RealValued` iff `f` is `RealValued` and
  `std::invoke_result_t<F, Real>` is real; otherwise `ComplexValued`.
  Implemented as a trait on `invoke_result_t`, not guessed.
- Invoked as `F(f[iTheta, iPhi])`. No point-dependent callables
  (`F(theta, phi, value)`) in phase 1 — if wanted later that is a second
  overload, not a change to this one.
- Regression to carry forward: an lvalue callable with observable copy
  semantics, proving the node owns a copy.

### 3.7 Grids

**Handle.** The grid is made **value-semantic and cheap to copy** in step −1:
`GaussLegendreGrid` becomes a handle holding `std::shared_ptr<const Impl>`,
where `Impl` owns the quadrature, the Wigner table and the FFTW plan cache
(`core-plan.md` step B). Terminals then hold a grid **by value**, and
expressions forward their left operand's copy. No `shared_ptr<const Grid>` in
the field layer, no `Grid&` member, and no hand-written assignment operators:
the whole class of lifetime and constness problems around the present
`_Grid& _grid` member disappears at the source rather than being wrapped.

Until step −1 lands, the field layer must not be written against
`GaussLegendreGrid` copies — the current class copies its entire Wigner table.
This is why `core-plan.md` step B is sequenced before phase-1 step 2.

**Identity.** `Grid::Identity()` returns `_impl.get()`. Binary nodes compare it
for equality at construction, in **all build modes**; failure throws
`std::invalid_argument`. Structural comparison is rejected: equal-parameter
grids with different wisdom or quadrature objects are different grids for our
purposes, and handle identity is the cheap, sufficient test. Currently `Grid()`
is taken from the left operand alone and a mismatch is silent.

**Support for `N`.** The check splits, and the split was previously mis-stated
as purely runtime:

- *Compile-time.* `NRange` is a template parameter of `GaussLegendreGrid`
  (`All`, `NonNegative`, `Single`). Whether the grid admits negative `N` at all
  is therefore a `static_assert`, not a runtime test.
- *Runtime.* `|N| <= grid->MaxUpperIndex()` (`GaussLegendreGrid.h:88`), in all
  build modes.

**Both checks live on terminals and views only.** Expression nodes range
freely over `N`. This is forced: on an `NRange = NonNegative` grid, a node-level
check would make `conj(f)` unconstructible at `N > 0`, killing
`Integrate(conj(f) * g)` (§3.9) even though it is pure pointwise arithmetic, and
breaking phase 4 outright, where the derived components of a real tensor are
`(-1)^N conj(stored)` over a scheme that stores `N ≥ 0`. On a `Single` grid any
`Sum`-rule node would be unconstructible. Terminals and views are exactly the
things that get stored and transformed, so they are exactly where the grid's
support matters. The transform boundary validates independently:
`ValidateTransformRequest` (`GaussLegendreGrid.h:327`) throws
`std::invalid_argument` when `|n| > lMax` or `n` is outside `UpperIndices()`, so
an unsupported `N` cannot escape into a transform regardless.

### 3.8 Materialisation, assignment and aliasing

`SpinField` gains:

- A constructor from any `SpinWeighted` expression with (a) equal
  `UpperIndex`, (b) `Value` convertible — `RealValued` expression into a
  `ComplexValued` field allowed, the reverse a compile error, (c) same `Real`.
  Grid handle taken from the expression. Evaluation goes through
  `EvaluateInto`.
- `operator=` from an expression, same constraints, plus runtime checks in all
  build modes: grid pointer equality with the destination, and size agreement.
  Assignment **never rebinds the destination grid**; write the test before the
  code.
- `+=`, `-=`, `*=`, `/=` from expressions and scalars, with the same index
  constraints as the binary forms (`+=` requires equal `N`; `*=` and `/=` by an
  expression require the right-hand side have `N == 0`, since the
  destination's `N` cannot change).

  *Implementation note (step 4).* Each compound operator is constrained by
  "can the binary expression be formed, and is the result assignable here",
  rather than by restating the index rules. That gets them right by
  construction — including that a complex scalar cannot multiply a real field
  in place, because the destination's value kind cannot change — and it keeps
  the failure at the call site, which is what makes
  `static_assert(!requires { u *= v; })` a usable negative test. Putting the
  check in the body instead makes every unlawful use a hard error, and the
  first draft did exactly that.

  **Small discrepancy, for the author.** "from expressions and scalars" cannot
  be met in full: `+=` and `-=` by a scalar would need an `f + s` node, and
  §3.5's inventory does not have one — it lists `s * f`, `f * s`, `f / s` and
  `s / f`, and says "nothing else". As implemented, `u += 2.0` simply does not
  compile, which is consistent with §3.5 and inconsistent with the sentence
  above. Adding `f + s` and `f - s` at `N == 0` would be a ten-line change and
  is arguably natural (adding a constant to a scalar field), but §3.5 is
  emphatic and this layer is deliberately small, so nothing was added. Worth a
  decision either way; it is not blocking.
- A free `Materialise(expr)` returning `SpinField`, for deliberately breaking a
  lazy chain — e.g. before feeding a product into a transform twice.

**Aliasing.** Because every node is pointwise and index-preserving (§3.5),
`u = expr` where `u` occurs in `expr` reads element `(iTheta, iPhi)` of `u`
only when writing that same element. In-place evaluation is therefore safe and
no temporary is needed. State this as a theorem-with-proof-sketch in the header,
tie it explicitly to the invariant, and add a `u = conj(u) * v + u` style
regression so that any future non-pointwise or re-indexing node that sneaks
into this layer fails a test rather than corrupting data silently.

### 3.9 Integration

**`Integrate(f)` is the only reduction phase 1 ships.** Free function,
constraint `N == 0` (theory note §5: `∫ f dΩ` vanishes identically unless
`N = 0`), returns `Scalar`.

`GridBase` supplies `Weights()` as `cartesian_product(CoLatitudeWeights,
LongitudeWeights) | transform(multiply)` in the canonical order, so a weighted
accumulation is directly available. But the longitude weights are
`views::repeat(dPhi, nPhi)` — *uniform* — so the quadrature factorises exactly
as `dPhi * Σ_θ w_θ (Σ_φ f)`, which is `nTheta` multiplies instead of
`nTheta·nPhi` and avoids iterating a `cartesian_product` view. Prefer the
factored form; keep `Weights()` for the generic case.

**No `InnerProduct`, no `Norm` at this level.** Both are deferred to the tensor
layer. The reasoning is not that they are hard but that at level 0 they are the
wrong objects: the pairing that matters is the metric one of theory note
`eq:metric`, carrying the `(-1)^α` factors, and it is a duality product between
a tensor and its dual rather than an L² inner product on a single component.
Naming the level-0 L² form `InnerProduct` would claim the name the tensor layer
wants, exactly as "canonical component" was claimed in §3.1.

What phase 1 keeps is the *fact* that made the level-0 version attractive:
`Integrate(conj(f) * g)` type-checks for **any** common `N`, because `conj(f)`
carries `-N` and the product lands at zero. That remains the showcase for the
index algebra and stays a test (family 3) — written out as the expression, not
hidden behind a name. Likewise `sqrt(Integrate(abs2(f)))` is available to any
caller who wants it; `abs2` survives on its own merits (§3.5), not as `Norm`'s
implementation detail.

One caveat to carry to whoever does write those reductions: `Integrate(abs2(f))`
for band-limited `f` of degree `l` integrates a degree-`2l` quantity, so on a
grid sized exactly for `lMax` it is where missing quadrature headroom first
bites. This is a grid-sizing question, answered in `core-plan.md` step D
(`ForBand`), not something the field layer should paper over by silently
refining.

### 3.10 Errors and diagnostics

- Unlawful index combinations: overload-resolution failure on constrained
  operators, backstopped by `static_assert` naming the rule violated
  ("addition requires equal upper indices; got N = +1 and N = -1").
- Grid mismatch, size mismatch, unsupported `N`: `std::invalid_argument`, all
  build modes. This class of check has been promoted out of assert-only twice
  in review; start there.
- Everything else (division by zero, NaNs): not checked; document.

---

## 4. What the existing code actually provides

Verified, so the plan does not have to guess:

| fact | where |
|---|---|
| Spatial layout is theta-major, phi fastest, flat index `iTheta*nPhi + iPhi` | `GaussLegendreGrid.h:169`, `GridBase::Points()` |
| `nTheta = lMax + 1` (Gauss–Legendre); `nPhi` = the least fast FFT length `≥ 2·lMax + 1` (`520` at `lMax = 256`). *Was `max(1, 2*lMax)`; changed in `core-plan.md` T3* | `GaussLegendreGrid.h`, `Utility.h` |
| `FieldSize()`, `Weights()`, `CoLatitudeWeights()`, `LongitudeWeights()`, `PointIndices()`, `ProjectFunction()` all exist on `GridBase` | `GridBase.h` |
| Longitude weights are uniform (`repeat(dPhi, nPhi)`) | `GaussLegendreGrid.h:103` |
| `MaxUpperIndex()` runtime; `NRange` compile-time (`All`/`NonNegative`/`Single`); `MinUpperIndex()`/`UpperIndices()` derived | `GaussLegendreGrid.h:88`, `GridBase.h:17` |
| Transform validates `lMax`/`n` by **throw**, sizes by **assert** | `GaussLegendreGrid.h:327` |
| Transforms take generic `InRange&&`/`OutRange&`, use `.size()` and iterators, and `std::copy` into an aligned work buffer when the input is not writable | `GaussLegendreGrid.h:124`, `:169–178` |
| Real-valued input is supported **only at `n = 0`**, with reduced `m ≥ 0` coefficient storage; a real transform at `n ≠ 0` throws. *Q5, done in `core-plan.md` T2* | `GaussLegendreGrid.h`, `tests/TestRealFieldSymmetry.cpp` |
| Complex-valued transforms require `MRange = All` | `GaussLegendreGrid.h:114` |
| Terminals allocate `FFTWpp::vector<Scalar>` — `fftw_malloc`, SIMD-aligned | `CanonicalComponentField.h:82`, FFTWpp `Core.h:116` |
| `Wigner` indexed `[n][iTheta][l][m]`, holding `sqrt((2l+1)/(4π)) · d^l_{nm}` — upper index **first**, per Dahlen & Tromp (C.115). Verified, §2 | `GaussLegendreGrid.h:181`, `Wigner.h:474–486` |
| ~~The `(lMax, lMax)` coefficient is explicitly zeroed on complex forward transforms, and again in the random-coefficient generators.~~ *Both workarounds deleted in `core-plan.md` T3; the orders `m = ±lMax` are now resolved.* | — |
| `FieldBase` (with the `CheckPointIndices` bug) is used only by superseded and dormant trees | `grep`; `GSHTrans/Field` |
| C++23, gtest via FetchContent, Eigen already a dependency | `CMakeLists.txt` |

Two rows carry consequences beyond their own line.

*The `(lMax, lMax)` zeroing — gone.* This row used to carry a caveat: a
round-trip through the complex forward transform at `lMax` would not reproduce
that mode, because with `nPhi = 2·lMax` the orders `m = ±lMax` were the same
discrete mode and the transform zeroed the coefficient rather than return a
wrong one. `core-plan.md` T3 sized `nPhi` at the least fast FFT length
`≥ 2·lMax + 1`, so both orders are now resolved and both workarounds are
deleted. **Test family 4's round-trip can be written at `lMax` as one would
naively expect**, which is what the caveat was there to prevent.

*Real input at any `n`.* This was a working, tested path that Q5 deliberately
removed (`core-plan.md` T2). Recorded here so that the deletion reads as a
decision rather than as an oversight when someone later finds the commit. What
was kept is the mathematics: `tests/TestRealFieldSymmetry.cpp` still exercises
`eq:basiclevel` at `n = 2` through complex transforms, which is the relation
phase 4 needs; what went is the storage scheme, which assumed the *self*-
relation that holds only at `N = 0`.

---

## 5. Consuming expressions without materialising

Fable's design nominates `EvaluateInto` as the seam through which phase 5 will
transform an expression directly. That is right, but the code offers a second
route that costs nothing to keep open: `ForwardTransformation` already accepts
an arbitrary `InRange&&`, requiring only `.size()`, iterators, and a value type
of exactly `Real` or `Complex`. When the input is not an `output_range` it
copies into its own aligned work buffer (`GaussLegendreGrid.h:174`). So a node
that models a sized random-access range of `Scalar` could be transformed today,
with no change to the grid.

Recommendation: do **not** put a range interface in the concept. It can be
added later as a free adaptor `Values(node)` built from `operator[]` and
`PointIndices()`, which touches no node. `EvaluateInto` stays the primary seam,
because it is the one that can write into a destination of a *different* scalar
type (§3.2) and into a slice of a larger buffer, which the range route cannot.

---

## 6. Tests

Four families, each with named acceptance criteria for the CLI.

1. **Compile-time algebra.** `static_assert` matrix over `N ∈ {-2,-1,0,1,2}`
   for every operator: lawful combinations satisfy `SpinWeighted` with the
   predicted `UpperIndex` and `Value`; unlawful ones fail `requires`. Must
   cover: `u + v` at unequal `N`; `real`/`imag` at `N ≠ 0`; division by
   `N ≠ 0`; complex-scalar promotion; `conj` double-negation returning the
   original `N`; `RealValued` claimed at `N ≠ 0` failing `SpinWeighted`; and
   the closure table of §3.4 discharged case by case, so that no lawful
   expression is collateral damage of that constraint.
2. **Lifetime.** Under ASan/UBSan: nested expressions bound to `auto`;
   expressions returned from functions; expressions stored in containers;
   **rvalue terminals in expressions** (`auto e = MakeField(...) + v;` — the
   case a plain `IsTerminal` trait misses); a named expression copied into a
   longer-lived expression; the lvalue-callable ownership regression of §3.6.
3. **Numerics.** Lazy vs materialised agreement to tight tolerance for a deep
   mixed tree at several `N`; the in-place aliasing cases of §3.8;
   `Integrate(conj(f) * g)` at `N ≠ 0` — first that it type-checks at all, then
   against the same quantity computed through explicit conjugate-component
   transforms (this doubles as the first oracle reuse for the phase-4 symmetry
   tests); `Integrate` of a constant against the sphere-area normalisation
   already pinned by the one-point grid.
4. **Interface stability.** `EvaluateInto` default equals the element loop;
   terminal override equals default; `EvaluateInto` into a `Complex`
   destination from a `RealValued` expression; layout matches transform input
   order by round-tripping through the scalar transform at `N == 0`; a minimal
   non-owning view over external storage satisfies `SpinWeighted` and
   participates in expressions interchangeably with an owning field. The
   `(lMax, lMax)` caveat of §4 no longer applies: `core-plan.md` T3 has landed,
   so the round-trip is written at `lMax` with the `m = ±lMax` modes included.

Sanitiser discipline as established: every family runs Debug and
Release+ASan/UBSan, leak detection off in the ptraced environment, and the
size/grid checks are exercised with assertions disabled.

---

## 7. Work order

Each step compiles and passes its tests before the next begins. Steps A–D of
`core-plan.md` precede step 1 here; step E onwards is sequenced against phase 5
and gates nothing below.

**All seven steps are complete.** Phase 1 is implemented, and the four test
families are in `tests/TestSpinField.cpp`, run in Debug, in Release, and under
ASan and UBSan. What phase 1 does *not* yet have is recorded in §11.

1. Concepts, `IsTerminal`, `OperandStorage`, index-rule tags, value traits.
   Pure headers plus compile-time tests (family 1 skeleton).
2. `SpinField` terminal: storage, shared grid handle, construction checks,
   `operator[]`, `EvaluateInto` override, mutable access. No algebra yet.
3. `Unary` and `Binary` templates plus the §3.5 operator set, constraints on
   the free operators, grid identity checks. Families 1 and 2 complete.
4. Materialisation, assignment, compound assignment, aliasing regressions.
   Family 3 (lazy-vs-materialised, aliasing) complete.
5. Callable nodes and their ownership regression.
6. `Integrate`, the minimal view fixture, and the transform round-trip.
   Families 3 and 4 complete. *(No longer gated: the `d`-function convention is
   settled and verified, §2.)*
7. Delete the superseded `CanonicalComponentField*` expression paths, or
   quarantine them behind a deprecation header if anything downstream still
   includes them; run the full pre-existing suite to confirm nothing outside
   the superseded layer regressed. Decide the fate of `FieldBase` here — its
   `CheckPointIndices` uses `<=` where it needs `<` (`FieldBase.h:29–32`), but
   it is used only by superseded and dormant code, so it is a delete-or-fix
   at this step rather than a phase-1 task.

   *Done.* Deleted outright, along with the five dormant trees, `FieldBase.h`
   and `ExpansionBase.h`; nothing needed quarantining, because nothing outside
   those trees included them. `FieldBase` was deleted rather than fixed.
   `CanonicalComponentExpansion*` stayed at this step, and needed one
   `#include "../Traits.h"` that it had been getting transitively through
   `ExpansionBase.h` — the same class of defect as core F10. It was deleted
   later, after core step F; §1 records why the argument for keeping it
   stopped holding.

   `examples/FieldExample.cpp` was rewritten against `SpinField` rather than
   deleted, so the public API is exercised from outside the test tree. Writing
   it turned up the index algebra doing its job: `scaled += conj(scaled) * (u *
   conj(u))` does not compile, because the right-hand side lands at `-N` rather
   than `N`.

---

## 8. Later phases, and the constraints they place on phase 1

Recorded so phase 1 does not foreclose them. None is implemented now.

### Phase 2 — tensor storage

*Prerequisite:* `core-plan.md` step F. Phase 2's `Layout` policy is designed
around what the transform can consume, and [C9] dropped the `PointMajor` repack
requirement on the strength of the batched, strided transform existing. Phase 1
never batches and is unblocked; phase 2 should not start until F has landed.


- `MultiIndex<Rank>` with constexpr conversion to and from a flat index, and
  `constexpr UpperIndex(MultiIndex)` implementing the signed sum of theory
  note `eq:N`.
- Constexpr orbit enumeration under index negation and under a permutation
  symmetry group, generating the stored-component table of theory note §3.1.
- `TensorField<Rank, Symmetry, Reality, Grid, Layout>` owning one contiguous
  buffer and handing out phase-1 component *views* at the correct upper index,
  rather than a tuple of separately allocated fields. The parameter is called
  `Reality` (`RealTensor`/`ComplexTensor`), *not* `Value`: phase 1's
  `Value ∈ {RealValued, ComplexValued}` is a statement about one component's
  scalar type, whereas a real tensor has one real component and many complex
  ones (theory note §7 item 6). Same word, different concept — keep them apart
  in the vocabulary.
- `Layout` is a policy, not a per-rank constant. `ComponentMajor` keeps each
  component contiguous, which is what the transforms need; `PointMajor` keeps
  all components at a point together, which is what a rank-4 tensor applied
  pointwise to a strain field wants. Rank-4 objects are needed both as local
  operators and, as inversion parameters, as transformable fields, so both
  layouts must exist at rank 4. `ComponentMajor` is the default. Transforms
  do **not** require it: the batched transform primitive is described by
  `(count, stride, dist)`, so `PointMajor` components are transformable in
  place with `stride = nComponents, dist = 1` (`core-plan.md` [C9]). A repack
  to `ComponentMajor` remains available and may well be faster; it is a
  performance choice to be benchmarked, not a precondition.

*Constraint on phase 1:* the component view must satisfy `SpinWeighted` while
pointing into a sub-range of someone else's buffer — hence §3.2's view
requirement and its test fixture. Per-component strides are phase 2's to choose
freely, which is true only because `core-plan.md` step C removes the transform's
alignment expectations from caller storage.

### Phase 3 — tensor algebra

Tensor product, contraction against the metric of theory note `eq:metric`, trace,
transpose, index permutation, (anti)symmetrisation, and rank-4 on rank-2 for
stress. Each delegates pointwise work to phase-1 nodes, so tensor expressions
stay lazy: a lazy tensor's component accessor returns a lazy component
expression.

*Constraint on phase 1:* index permutation permutes *tensor slots*, not grid
points, so it does not violate the index-preserving invariant of §3.5. Keep it
that way — the invariant is about `(iTheta, iPhi)`.

### Phase 4 — reality reduction

Store one component per negation orbit (theory note §3.1–3.2); expose the
derived components through read-only views applying `T^{-α} = (-1)^N conj(T^α)`
on the fly, so generic code can traverse all `3^Rank` components without
knowing which are stored. Halves storage for every real tensor: 81 reals rather
than 162 at rank 4.

*Done, and §14 records it.* The reduction is exact rather than approximate —
`3^p` reals per point at every rank and symmetry — which took one thing the
phase-2 plan did not anticipate: **the storage stops being uniform**.

*Constraint on phase 1:* a derived component is a `Unary<Conj, Negate, View>`
composed with a sign — that is, a phase-1 expression node over a phase-1 view,
which the design already supports. Two consequences worth recording now:

- Traversal over "all components" is necessarily a *compile-time* traversal
  (a template for-loop), because stored and derived components have different
  types. Runtime loops over components are not available and should not be
  designed for.
- The derived component of a stored `N ≥ 0` component has `N ≤ 0`. This is the
  concrete case that settled §3.7: expression nodes carry no grid `N`-support
  check, or phase 4 breaks on any grid with `NRange = NonNegative`.
- Orbit representatives follow **Dahlen & Tromp (1998)**, with orthonormalised
  harmonics — taken to agree with Phinney & Burridge up to that normalisation.
  Pin the correspondence explicitly in a comment when phase 4 is written: it is
  the one place where the convention becomes observable in *stored data* rather
  than in an intermediate.

### Phase 5 — spectral side

Tensor-valued expansions, per-component transforms, then the raising and
lowering operators of theory note §6. The node interface of §3.2 is chosen so
this does not require revisiting it.

**Batched transforms.** The transform primitive is "transform a batch of `k`
same-spin fields", with the single field as `k = 1` — not a scalar primitive
looped from outside. The batch is public API, described by
`(count, stride, dist)`, and threads only when explicitly asked
(`core-plan.md` [C9]); a batch shares grid, `lMax` and `n`, since the Wigner
block is what is being amortised. The FFT stage batches across `nR × nTheta` rows
via FFTW's advanced interface, and the Wigner stage at fixed `(N, θ)` applies
the same `d`-values across all radii, turning a bandwidth-bound matrix–vector
contraction into a matrix–matrix one whose arithmetic intensity grows with
`nR`. This is likely the largest performance lever in the library. Benchmark it
jointly with any exploitation of the ±n Wigner symmetry, not separately.

Known prerequisites in the existing code, now planned as `core-plan.md` steps
C, E and F: per-call plan creation and work-buffer allocation move to a
grid-owned cache; the transform stops executing plans on caller storage; and the
transform primitive becomes "transform `k` contiguous same-spin slices" rather
than a scalar primitive looped from outside. None of that is required before
phase 1, but the field-layer interface of §3.2 is chosen so it can land
underneath without being revisited.

### Layered (3D) fields

Applications build 3D fields as the product of the angular grid with a set of
radial nodes (radial SEM in practice, but only the outer index matters here).

**2D is the primitive; 3D is a stack.** The angular field is first-class and
never wrapped. A 3D field is one contiguous buffer viewed as `nR` angular
slices, radius-major (`[r][iTheta][iPhi]`, each slice contiguous in the
canonical order), exposed through a `Layered<...>` wrapper whose slice accessor
returns a phase-1 node — a non-owning view for terminals, a slice expression
for lazy nodes. The index algebra, value traits and constraints lift unchanged,
since `UpperIndex` and `Value` are properties of the type, not the slice. There
is no second expression system. 2D is *not* represented as `nR == 1`: that
idiom stays available to application code, but the library treats a one-slice
stack and an angular field as distinct types.

**Bridges are explicit.** `Broadcast(expr2d)` lifts a 2D expression to a
layered node returning the same expression at every radius (zero-copy; this is
how `field3d * field2d` works — surface masks, r-independent coefficients).
`Slice(field3d, r)` is the explicit way down. No implicit conversions in either
direction.

**Parallelism.** OpenMP over slices (later over `(component, r)` pairs), each
thread calling `EvaluateInto` on its slice of the destination — legitimised by
the §3.2 thread-safety contract. Pointwise 3D work is memory-bandwidth bound at
production sizes (a single complex scalar component at `lMax = 256`, `nR = 100`
is ~210 MB), so lazy evaluation is worth real memory traffic, and on NUMA
machines the destination should be first-touched in the same
parallel-over-slices structure as the compute. Nested parallelism with
`Wigner::ComputeAll`'s internal parallel loop must be explicitly excluded (one
parallel region at a time).

**Spectral-side layout.** Angular transforms want `[r][(l,m)]`; radial
operations (SEM derivatives, mode-wise radial solves) want `[(l,m)][r]`. Both
are hot, so the 3D expansion takes the same `Layout` policy plus explicit
repack planned for tensors. Radial derivatives and interpolation are linear in
`r` and non-pointwise, hence excluded from the lazy layer by the §3.5
invariant — same status as raising/lowering. The radial mesh abstraction needs
only nodes, quadrature weights and handle identity; SEM connectivity stays in
the application layer. Name the wrapper neutrally (`Layered`, not `Radial`):
time levels and ensembles want the identical structure.

### Deliberately deferred, with the seam named

- **Transforms of expressions** — phase 5 consumes `EvaluateInto`; nothing here
  may assume materialised input.
- **Dealiasing** — the product of two band-limited fields exceeds the grid's
  truncation. This is a grid-level question, flagged against
  `GaussLegendreGrid`, not something the field algebra decides.
- **Vectorised evaluation** — `EvaluateInto` is the hook; no SIMD work now.
- **`3j`, Eigen-backed pointwise linear algebra at ranks 2 and 4** — Eigen is
  already a dependency; phase 3 may use it.

---

## 9. Decisions taken

All twelve questions are answered. Two answers went **against** the
recommendation (Q5, Q10) and one went **beyond** it (Q3); those three are the
ones that changed the shape of the work, and between them they are the reason
`core-plan.md` exists.

| | question | decision | where it lands |
|---|---|---|---|
| Q1 | naming | `SpinField`, concept `SpinWeighted` | §3.1 |
| Q2 | `EvaluateInto` output scalar | templated on `S` | §3.2 |
| Q3 | grid ownership | **make the grid itself lightweight** | §3.7, core step B |
| Q4 | FFTW alignment | stop new-array execute on caller storage; keep the field layer agnostic | §3.2, core step C |
| Q5 | `RealValued` at `N ≠ 0` | **forbid — and cut it from the library** | §3.2, §3.4, core step A |
| Q6 | grid `N`-support check | terminals and views only | §3.7 |
| Q7 | pointwise invariant | "pointwise and index-preserving" | §3.5 |
| Q8 | `abs2`, factored quadrature | both adopted | §3.5, §3.9 |
| Q9 | callable node name | `Map` | §3.5, §3.6 |
| Q10 | `InnerProduct` name | **no `InnerProduct`/`Norm` at level 0 at all** | §3.9 |
| Q11 | tensor-level `Value` | `Reality` | §8, phase 2 |
| Q12 | orbit representatives | Dahlen & Tromp (1998), orthonormalised | §8, phase 4 |

### Q5 — `RealValued` at `N ≠ 0`: forbidden, and removed from the library

*Answer:* "there is a good argument for cutting this feature, eventually from
the whole library. I would rather have it be correct mathematically, and then
consistent, than avoid a rewrite."

This overrides the earlier recommendation to allow it at level 0, and it is the
right call: the recommendation was an argument from what the code happens to
support, against a statement the theory note already makes globally (§7 item 5),
not merely about the tensor layer. Real-valuedness of a spin-weighted field is
not preserved by the frame rotation `e_± ↦ e^{∓iψ} e_±`, so it is not a property
any component of any tensor can have at `N ≠ 0`. A library that can represent it
can represent something that does not exist.

Three consequences, in increasing order of reach.

1. **In the concept** (§3.2): `Value == RealValued ⟹ N == 0`, checked on
   `SpinWeighted` itself.
2. **In the algebra** (§3.4): the constraint is closed under every node, so it
   is enforced once at the concept and never re-derived. The closure table is
   the evidence that this costs no expressiveness. Worth noticing that the
   *reverse* implication is what makes the design pleasant: `abs`, `abs2`,
   `real`, `imag` and `Map` all land at `N = 0` precisely because that is the
   only place their results could be real.
3. **In the core** (`core-plan.md` step A): the reduced `m ≥ 0` coefficient
   path at `n ≠ 0` becomes unreachable from the field layer, and should be
   removed rather than left as dead weight. This is a deliberate deletion of
   working, tested code — `tests/TestRealFieldSymmetry.cpp` at `n = 2` goes with
   it, or is rewritten at `n = 0`.

The one thing worth stating carefully, so the deletion does not overreach: the
relation that test exercises, theory note `eq:basiclevel`, is **true and
useful**. What it says is that the coefficients of `conj(f)` at `-N` are
determined by those of `f` at `+N`. That is a relation between two *different*
fields, and it is the engine of phase 4's reality reduction. What is *not*
implied by it, and what the reduced `m ≥ 0` storage at `n ≠ 0` silently assumes,
is the self-relation `f^N_{l,-m} = (-1)^{m-N} conj(f^N_{lm})` — which holds only
if `f` is its own conjugate, i.e. only at `N = 0`. Delete the storage scheme;
keep the relation, and reuse the test as an oracle for phase 4.

### Q3 — grid ownership: make the grid lightweight

*Answer:* "I suspect it will be cleaner to make the grids lightweight as you are
suggesting."

Taken up in full: `GaussLegendreGrid` becomes a handle over a shared immutable
`Impl`, so grids are values, copying one is a pointer copy, and "same grid" is
handle identity (§3.7, `core-plan.md` step B). The two commented-out members at
`GaussLegendreGrid.h:346–347` show this was already the intended direction.

This is better than the recommended `shared_ptr<const Grid>` in the field layer
for a reason worth recording: with `shared_ptr` in the field layer, *every*
consumer of a grid — fields, expansions, the layered wrapper, application code —
independently decides how to hold it, and they will not all decide the same way.
Putting the indirection inside the grid makes the question unaskable. It also
gives the FFTW plan cache of `core-plan.md` step E somewhere to live that is
shared by construction rather than by convention.

### Q10 — no `InnerProduct` or `Norm` at level 0

*Answer:* "we keep this stuff at the Tensor level. It is not essential for the
spin fields. And here it's more of a duality product anyway."

Sharper than the recommendation, which was to keep the level-0 one and rename
the tensor-level one. §3.9 is rewritten accordingly: `Integrate` is the only
reduction phase 1 ships, and `Integrate(conj(f) * g)` stays as a *test*, which is
where its value was anyway — it demonstrates the index algebra, and demonstration
does not need an API. The observation that the level-0 pairing is a duality
product rather than an inner product is the substantive point: with the metric
of theory note `eq:metric` the natural pairing at the tensor level contracts a
tensor with a dual tensor and carries `(-1)^α`, and calling the level-0 L² form
`InnerProduct` would have made that harder to name honestly later.

### Q4 — alignment

*Answer:* "Yes, I think I agree with your recommendation. But we can also try to
make this layer somewhat agnostic, and later refine the FFT or GSHT steps to
optimise if that seems worthwhile."

Adopted with the qualification, which is the important half: the *contract* is
that the field layer imposes no alignment on slice targets, and the
*implementation* satisfies it by copying rows into the plan's own aligned
buffers. Whether that copy is later avoided — via `FFTW_UNALIGNED` plans, via a
second plan variant selected on the caller's actual alignment, or by folding the
copy into the batched pack stage where it is free — is a question for the core
and does not reach the field layer. `core-plan.md` step C does the safe thing;
step F is where the copy stops being a cost at all.

### Q1, Q2, Q6, Q7, Q8, Q9, Q11, Q12 — recommendations adopted

- **Q1** `SpinField` / `SpinWeighted`, per §3.1.
- **Q2** `EvaluateInto` templated on the output scalar, per §3.2. Decided now
  because adding to that member later touches every node.
- **Q6** the grid's `N`-support check lives on terminals and views only, per
  §3.7; phase 4 is the case that forces it.
- **Q7** the invariant is stated as "pointwise **and index-preserving**", per
  §3.5, so that a future phi-shift or transpose view cannot silently invalidate
  the aliasing argument of §3.8.
- **Q8** `abs2` is a node in its own right and `Integrate` uses the factored
  quadrature, per §3.5 and §3.9.
- **Q9** `Map`, per §3.5.
- **Q11** the tensor-level parameter is `Reality`, per §8.
- **Q12** Dahlen & Tromp (1998) orbit representatives with orthonormalisation,
  per §8. Recorded now, needed at phase 4; the correspondence with Phinney &
  Burridge is to be pinned in a comment at the point of use.

---

## 10. What was open, and how it closed

Both items are now closed. Neither leaves residual work in this document.

1. **The Wigner value convention — closed, verified.** The stored values are
   `sqrt((2l+1)/(4π)) · d^l_{Nm}` with the upper index first, per Dahlen &
   Tromp (1998) eq. (C.115), confirmed against all nine `l = 1` values. §2
   carries the statement, the table and the corroborating evidence;
   `core-plan.md` §6 [C8] carries the same fact from the core's side, and §8
   task T1 turns it into a permanent test.

   This unblocks §7 step 6 and the transform-based oracles of test family 3,
   which were the only things it gated.

2. **How far the Q5 deletion goes into the grid's template parameters —
   closed in the core plan.** `MRange` survives, documented as "real scalar
   grid" rather than as an order range, with a constructor-time rejection of
   `MRange = NonNegative` combined with `nMax != 0`, and a scheduled revisit at
   step G where the storage question is opened anyway. `core-plan.md` §6 [C1].

   Decided alongside it, and worth knowing here because it changes a public
   template signature: the `Normalisation` axis is **deleted** — `Ortho` is the
   only normalisation, `FourPi` goes (`core-plan.md` §6 [C7], step A2). Any
   phase-1 test that instantiates `Wigner` directly should not name a
   normalisation.

What remains genuinely undecided in the project is all in `core-plan.md` §7,
and all of it is performance work deferred to measurement: Wigner storage
layout and symmetry reduction, a runtime `NRange` set, stored-value precision,
and the fate of `RowMajor`. None of it reaches this document — the phase-1
interface of §3.2 was chosen so that it can all land underneath without being
revisited.

---

## 11. What phase 1 has, and what it does not

Recorded at the end of phase 1 so that phase 2 starts from facts rather than
from this document's intentions.

**Implemented and tested.** The `SpinWeighted` concept with the reality
constraint; `SpinField` and `SpinFieldView` as terminals; `Binary` and `Unary`
with the six index rules; the §3.5 operator set; `Map`; construction,
assignment and compound assignment from expressions; `Materialise`;
`Integrate`. Four test families in `tests/TestSpinField.cpp`, run in Debug, in
Release, and under ASan and UBSan.

**Deliberately absent, per the plan.** `InnerProduct` and `Norm` (§3.9);
`pow`, `exp`, `log` as named nodes, which are all `Map` (§3.5); a range
interface on the concept, which stays a possible free adaptor `Values(node)`
(§5); transforms consuming expressions directly, which is phase 5's use of
`EvaluateInto`.

**Absent and worth knowing.**

- *No CRTP convenience base.* §3.2 permitted one for `Size()`, `Integrate`
  forwarding and point iteration. None was written: the only shared code was
  the default evaluation loop, which is a free function
  (`EvaluateNodeInto`) that each node calls. Adding a base later is possible;
  nothing depends on its absence.
- *The thread-safety contract is documented but not tested.* §3.2 guarantees
  that concurrent evaluation of one node into disjoint targets is safe, and
  the code has no mutable state to violate it, but no test exercises it. The
  layered layer is the first thing that will rely on it, and that is where a
  test belongs.
- *`operator[]` does not bounds-check in release builds.* Point indices are
  `assert`-only, deliberately: §3.10 lists grid, size and upper-index checks as
  throwing in all build modes and everything else as unchecked.
- *`f + s` and `f - s` do not exist*, so `u += 2.0` does not compile. See the
  note in §3.8; it is a decision for the author, not an oversight.

**Constraint on phase 2, restated.** `core-plan.md` step F must land first:
phase 2's `Layout` policy was designed around the batched, strided transform
when [C9] dropped the `PointMajor` repack requirement.

---

## 12. Phase 2 in detail

Written when phase 2 was started, from §8's sketch plus the theory note and the
phase-1 code as they actually are. §8 remains the statement of intent; this is
the working plan.

### 12.1 What the theory note requires

Three facts drive everything here, and all three are in
`canonical-components.tex`:

- **A component is labelled by a multi-index, and its upper index is the signed
  sum** `N = Σ αᵢ` (`eq:N`). For rank ≥ 2 the two are different things: a
  rank-2 tensor has three components at `N = 0`, so a collection labelled only
  by `N` does not determine a tensor. This is why phase 1's unit is called
  `SpinField` and not `CanonicalComponent`.
- **Permutation symmetry reduces the stored set** to one component per orbit of
  the permutation group, with a sign when the group is antisymmetric.
- **Reality reduces it further** — one per orbit of index negation, by
  `T^{-α} = (-1)^N conj(T^α)` (`eq:reality`) — but that is phase 4, and §8
  keeps it there. Phase 2 must not foreclose it.

### 12.2 The four steps

**Step 1 — the index machinery.** Pure compile-time, no storage, no grid.

- `MultiIndex<Rank>`: an array of `αᵢ ∈ {-1, 0, 1}`, with constexpr conversion
  to and from a flat index in `[0, 3^Rank)`, `UpperIndex()` implementing
  `eq:N`, `Negated()`, and application of a slot permutation.
- `SlotPermutation<Rank>`: a permutation of the slots together with a sign,
  `+1` for a symmetric generator and `-1` for an antisymmetric one.
- A `Symmetry` policy supplying a constexpr list of generators.
  `NoSymmetry`, `Symmetric<Rank>`, `Antisymmetric<Rank>`, and `GeneratedBy<…>`
  for anything else — the elastic-tensor symmetry
  `c_{ijkl} = c_{jikl} = c_{ijlk} = c_{klij}` is three generators and is the
  case rank 4 exists for.
- Constexpr **orbit enumeration** over the group generated by those, producing
  for each of the `3^Rank` components: which component is stored for it, and
  the sign relating them. Written so that negation can join the generating set
  without changing the algorithm, which is what phase 4 will do.

The acceptance tests are the theory note's own tables: 5 stored components for
a real rank-2 tensor under negation alone, 41 at rank 4, and the worked
symmetric rank-2 example whose four orbits give 6 reals per point. Those
numbers are computed by the enumeration, not asserted against a hand-written
list.

**Step 2 — `TensorField` in `ComponentMajor`.** One contiguous buffer of
`nStored × FieldSize`, handing out phase-1 nodes per component:

- a **stored** component is a `SpinFieldView` at the right upper index,
  pointing into the buffer;
- a component related to a stored one by a `+1` permutation **is** that same
  view;
- one related by `-1` is a phase-1 expression, `-1 × view`.

So `Component<MultiIndex>()` returns different types for different components,
and traversal over all `3^Rank` of them is a compile-time loop. §8 already
records this for phase 4; it starts at phase 2, because antisymmetry has the
same consequence as reality does.

**Step 3 — `PointMajor`.** The same buffer indexed the other way, which makes
each component's samples *strided*. `SpinFieldView` holds a `std::span` and
assumes contiguity, so this step is where that changes — a stride defaulting
to one, with the contiguous fast path kept in `EvaluateInto`. Deferred to its
own step so that the decision is taken with the rest of phase 2 built, and
because `core-plan.md` [C9] means a strided component is transformable without
a repack, so nothing forces the layout.

**Step 4 — the transform bridge.** A batched transform over the components
sharing an upper index, using the `(count, stride, dist)` descriptor: `stride =
1, dist = FieldSize` in `ComponentMajor` and `stride = nStored, dist = 1` in
`PointMajor`. This is what step F was built for, and it is the first thing to
consume it.

*Steps 3 and 4 were taken in the other order*, since step 4 is what makes a
tensor field useful and step 3's motivation — a rank-4 tensor applied
pointwise — belongs to phase 3. Nothing depended on the order.

### 12.4 What phase 2 landed, and the one thing that was not foreseen

All four steps are in, with 127 tests running in Debug and under ASan and
UBSan. The index machinery reproduces every table in the theory note and, as a
check nobody wrote it for, gives **21 independent components for the elastic
symmetry** `c_{ijkl} = c_{jikl} = c_{ijlk} = c_{klij}`.

**The buffer is ordered by upper index, not by flat multi-index, and that was
forced rather than chosen.** A batch is described by `(count, stride, dist)`,
so its members must be uniformly spaced. In flat order the stored components
sharing one upper index are scattered at no fixed spacing as soon as there is
any symmetry — so no single descriptor covers them, and the batching that step
F exists for would have been unreachable from the layer it was built for.
Ordering by upper index makes each group a contiguous run on both sides.

This is worth recording because §8 wrote the batch descriptor and the tensor
layout as independent decisions, and they are not: **the descriptor's
uniformity requirement is a constraint on tensor storage.** `Orbits.h` stays in
flat order, which is pure combinatorics; the layout belongs to the field.

Two smaller things learned:

- **A component's accessor cannot return one type.** Stored components are
  views, symmetric-permutation relatives are the same view, antisymmetric
  relatives are an expression, and a vanishing orbit has nothing to return.
  §8 recorded this as a phase-4 consequence of reality; it arrives at phase 2,
  because antisymmetry has it too.
- **Constraints must be `requires`-clauses, not `static_assert`s**, or the
  negative tests are vacuous — and the pack size has to be checked before the
  multi-index is formed, or a wrong-length component is a hard error inside
  `std::array` rather than an unsatisfied constraint.

### 12.3 Decisions taken here

**`Reality` is carried but inert until phase 4.** §8 names it as a phase-2
template parameter and puts the reduction in phase 4. So phase 2 stores one
component per *permutation* orbit and `RealTensor` stores the same set as
`ComplexTensor`; phase 4 adds the negation orbits to the generating set, and
the enumeration of step 1 is written so that this is the only change. The
alternative — reducing on reality now, since the storage layout is being
decided anyway — was rejected as a unilateral reordering of agreed phases, but
it is the one place where following the plan costs some rework, and it is
recorded here so the choice is visible rather than accidental.

**Components are addressed by multi-index at compile time**, i.e.
`t.Component<-1, 0>()`, not by a runtime tuple. The upper index is a template
parameter of the returned node, so it cannot be otherwise.

**No lazy tensor expressions in phase 2.** That is phase 3. Phase 2 hands out
phase-1 nodes and nothing more, so every operation on a tensor in phase 2 is
written by the caller as operations on components.

---

## 13. Phase 3 in detail

The tensor algebra: product, contraction, trace, permutation, symmetrisation,
and rank-4 on rank-2. Written when phase 3 was started, on the same footing as
§12.

### 13.1 The shape of it

Phase 2 gave a tensor a component accessor returning a phase-1 node. Phase 3's
whole content is that **an operation on tensors is an operation on that
accessor**. A tensor expression is any object with a rank, a grid and a
`Component<Alphas...>()`; the operations build new ones whose accessors compose
their operands', and no arithmetic happens until a component is evaluated.

So the layer is lazy for free, and it needs no second expression system:
every node it returns is a phase-1 node, and the index algebra that governs
them is already enforced. The tensor layer's own correctness is entirely about
*which* components it asks for.

| operation | rank | component |
|---|---|---|
| `Permute<π>(T)` | `p` | `T^{π(α)}` — a relabelling, no arithmetic |
| `S ⊗ T` | `p + q` | `S^α T^β`, a phase-1 product; upper indices add by `eq:N` |
| `Contract<j,k>(T)` | `p − 2` | `Σ_α (−1)^α T^{…α…−α…}`, a phase-1 sum of three terms |
| `Trace(T)` | `p − 2` | `Contract<0,1>` |
| `Symmetrise(T)` | `p` | `(1/\|G\|) Σ_π sign(π) T^{π(α)}` |

Two things fall out rather than being arranged. **Contraction lands at the
right upper index by construction**: the contracted pair contributes
`α + (−α) = 0`, so all three terms share an upper index and phase 1's `Equal`
rule admits their sum — a contraction that paired the slots wrongly would fail
to compile. And **permutation cannot break phase 1's aliasing theorem**,
because it permutes tensor slots and not grid points; §8 records this and it is
what keeps the lazy layer sound.

### 13.2 The steps

1. **The concept, and `Permute`.** `TensorExpr`; `IsTerminalTrait` for
   `TensorField` so that operand storage follows phase 1's rules exactly;
   permutation as pure relabelling, with `Transpose` for rank 2. No arithmetic,
   so it is the step that pins the machinery for expanding a compile-time
   index array back into a template pack.
2. **The tensor product.** Splitting a component's multi-index between two
   operands.
3. **Contraction and trace**, against the metric `g_{αβ} = (−1)^α δ_{α+β,0}`.
4. **(Anti)symmetrisation**, as a sum over the symmetry group's elements.
5. **Materialisation** into a `TensorField`, which is where a lazy tensor stops
   being lazy.

*All five are in, with 145 tests.* The composite they were for —
`Contract<2,3>(Contract<2,4>(c ⊗ e))`, an elastic tensor applied to a strain —
is written in the test file in three lines, and every intermediate is checked
by the compiler for rank and upper index.

Nothing in phase 3 needed a change to phase 1 or phase 2. That was the bet §8
made when it said "each delegates pointwise work to phase-1 nodes", and it paid
off: the only new machinery is turning a computed multi-index back into
template arguments, which `TensorDetails::ComponentOf` does once and everything
else uses.

### 13.3 Decisions taken here

**Operand storage follows phase 1 exactly** — `OperandStorage`, `IsTerminal`,
and the rule that an lvalue terminal is held by reference and everything else
by value. The hazard that leaves is the same one Eigen has and phase 1
documented; a second, different rule at the tensor level would be worse than
the hazard.

**A materialised expression has no symmetry unless it is asked for.** The
product of two symmetric tensors is not symmetric, and inferring symmetry from
an expression tree is a research problem rather than a design. `Materialise`
therefore produces `NoSymmetry` by default and takes the symmetry as an
explicit template argument when the caller knows better — at which point it is
the caller's assertion, checked by construction rather than by inspection.

**Vanishing and derived components propagate as they must.** A component of a
product is derived if either factor's is; a component of a contraction is a sum
whose terms may individually be stored, derived or vanishing. Nothing here
needs to know: it asks the operand for a component and gets a node, and the
`Vanishes` trait is what a traversal consults before asking at all.

---

## 14. Phase 4, and the thing phase 2 did not anticipate

Turning the switch on in `Orbits.h` was, as §12.3 promised, a one-line change:
`RealTensor::ReducesOnReality` becomes true, negation joins the generating set,
and the orbits get larger and fewer. Everything else followed from what phase 2
had already built — except one thing.

### The storage stops being uniform

An orbit that contains its own negation pins its component to a **single real
number**: real, or purely imaginary when a permutation sign gets in the way.
Phase 2's buffer is `Complex`, and storing a pinned component there would waste
half of each *and*, worse, let a caller write an imaginary part into a
component that cannot have one — which is exactly the class of thing this
library's type system exists to prevent.

So a real tensor has **two buffers**: complex components grouped by upper index
as before, and pinned components in a buffer of `Real`. That is what makes the
reduction exact rather than approximate:

| | components | reals per point | naive complex |
|---|---|---|---|
| rank 1 | 1 complex + 1 real | 3 | 6 |
| rank 2 | 4 complex + 1 real | 9 | 18 |
| rank 2, symmetric | 2 complex + 2 real | 6 | 18 |
| rank 2, antisymmetric | 1 complex + 1 imaginary | 3 | 18 |
| rank 3, symmetric | 4 complex + 2 real | 10 | 54 |
| rank 4 | 40 complex + 1 real | 81 | 162 |
| rank 4, elastic | 8 complex + 5 real | 21 | 162 |

Every one of those is the count of independent real degrees of freedom, and
none of them is special-cased.

**The pinned components are always at upper index zero**, which is what makes
this work at all: permutation preserves the slot sum and negation reverses it,
so a component fixed by the combination satisfies `N = −N`. Phase 1 forbids
`RealValued` anywhere but `N = 0`, so the two constraints agree exactly — and
the real buffer is therefore *one* transform group rather than several,
running through the transform's real path and its reduced `m ≥ 0` coefficient
storage. The saving in the spectral domain is the same saving as in the
spatial one.

### What phase 1 and phase 2 got right in advance

- **`conj` reverses the upper index.** A derived component is
  `sign · conj(view)` at `−N`, so getting this wrong — one of the three
  structural defects in the layer phase 1 replaced — would have made phase 4
  impossible rather than merely wrong.
- **The grid's `N`-support check is on terminals and views only** (§3.7). The
  derived partner of a stored component at `N` sits at `−N`, so on a grid
  carrying only non-negative upper indices half of them could not be
  terminals. They are expressions, and expressions are unchecked.
- **A component accessor returning different types** was already true at phase
  2 through antisymmetry, so phase 4 added a third and fourth case rather than
  a new idea.

### One decision taken here

**`Materialise` takes the reality as well as the symmetry**, defaulting to
`ComplexTensor`. Where the target component is pinned and the expression is
complex, the surviving part is taken — the real part, or the imaginary part
for an `Imaginary` orbit. That is the caller's assertion that the expression
really is a real tensor, of exactly the same kind as asking for a symmetry and
unprovable here for the same reason. Without it the algebra could not produce
a real tensor at all, which would have made the reduction useful only for
inputs.

---

## 15. Phase 5 in detail

The spectral side: expansions as objects rather than raw coefficient buffers,
and the operators that connect different upper indices.

### 15.1 What is actually needed

Phase 2's transform bridge writes into a `std::span<Complex>` the caller
allocates from `CoefficientSize(lMax)`. That works and is untyped: nothing
stops a caller indexing the wrong block, and there is nowhere for `ð` to live.
Phase 5 gives that span a type.

**The storage does not split the way phase 4's did.** A pinned component is a
real *field*, but its coefficients are complex numbers in the reduced `m ≥ 0`
storage — so the expansion side is one complex buffer throughout, with blocks
of differing lengths. That is what `CoefficientSize(lMax)` already computes.

### 15.2 The steps

1. **`SpinExpansion`** — the spectral counterpart of `SpinField`: a degree, an
   upper index, and `operator[](l, m)`. Built on `GSHView`, which `Views.h`
   already provides over coefficient storage, wrapped so that the upper index
   is a compile-time property and the reduced `m ≥ 0` storage of a real field
   is a type distinction rather than a convention. Plus the transforms between
   it and `SpinField`.
2. **`TensorExpansion`** — one buffer, per-component expansion views, and the
   batched transform to and from `TensorField`. Mirrors phase 2 but without
   the second buffer.
3. **`ð` and `ð̄`** — raising and lowering, which in the spectral domain are a
   multiplication by an `l`-dependent factor with no coupling between different
   `(l, m)` (theory note `eq:eth`).

### 15.3 The convention question, and what can be checked

The theory note flags the **overall signs** in `eq:eth` as needing to be fixed
against Phinney & Burridge, and says the magnitudes are not in doubt. §2 of
this document has carried that as one of the two remaining unchecked
conventions.

Implementing the note's stated signs and testing what they imply is the honest
position, and two identities are checkable without settling the convention:

- **`ð̄ð` is the surface Laplacian on a scalar.** The factors multiply to
  `−l(l+1)`, which is the eigenvalue of `∇²` on `Y_{lm}`. This pins the
  *relative* sign of the two operators and both magnitudes.
- **`[ð, ð̄] = −2N`.** A structural identity, independent of the overall sign.

Flipping the sign of *both* operators leaves both identities intact, so the
overall sign remains a convention this document cannot settle from inside.
It is recorded at the point of use, and it is the one thing in phase 5 that
wants an answer from the literature rather than from a test.

### 15.4 What phase 5 landed

All three steps, with 170 tests. `SpinExpansion` and `TensorExpansion` on one
side, `Expand`/`Evaluate` between the domains, and `Raise`/`Lower`.

**The expansion side really is one buffer**, as §15.1 predicted: a pinned
component is a real field but its coefficients are complex, so what varies is
the block *length* and not the buffer's type. The reduction carries across —
a real rank-2 tensor stores five blocks rather than nine here too, and a
pinned block is about half the length of a complex one at the same degree.

**A derived component is not offered in the spectral domain.** On the spatial
side deriving one is pointwise, `T^{-α} = (-1)^N conj(T^α)`, so a view plus an
expression node does it. On the spectral side it is `eq:complevel`,
`T^{-N}_{l,-m} = (-1)^m conj(T^N_{lm})`, which *reverses the order index* —
not a view over anything. So `TensorExpansion` exposes only stored components,
and deriving is done in the spatial domain or by applying the relation to the
representative. This asymmetry between the two domains was not foreseen and is
the one place the mirror between `TensorField` and `TensorExpansion` breaks.

**The `ð` sign convention is the one open item in this document.** §2 listed
two unchecked conventions; the `e_±` sign remains untested because nothing yet
observes it, and the `ð` signs are now implemented as the theory note states
them, with `ð̄ð = ∇²` and `[ð, ð̄] = −2N` tested. Both survive flipping the pair,
so the common sign wants an answer from Phinney & Burridge rather than from a
test.

---

## 16. The contravariant derivative

Added after phase 5, on the observation that `ð` is not the operator a
geophysical user of tensor fields actually wants. `core-plan.md` is unaffected.

### 16.1 What is being built

D&T (C.151)–(C.153), without the radial part. The surface gradient `∇₁` has no
`e₀` component at all (D&T C.145), so an angular library can implement it
completely; the `σ = 0` component is `∂/∂r`, which a user supplies externally
and which folds into 3D fields later. That is a scope decision, not a gap.

For a rank-`q` tensor expansion, `SurfaceGradient` produces a rank-`(q+1)` one
with

```
(∇₁T)^{σ α₁…α_q}_{lm} = Ω^{∓N}_l T^{α₁…α_q}_{lm} − Σ_i T^{α₁…(α_i+σ)…α_q}_{lm}
```

for `σ = ±1`, and zero for `σ = 0`, where `Ω^{±N}_l = √(½(l±N)(l∓N+1))` and
the upper sign goes with `σ = −1`. Coefficients whose shifted index leaves
`{−1,0,1}` are zero, as are those at `l` below the component's own `|N|`.

### 16.2 The one thing that needs building first

Every term above is a coefficient of some component at fixed `(l, m)` — but
the shifted components need not be *stored*. Phase 5 exposed only stored
components in the spectral domain, precisely because deriving one there is
`eq:complevel`, which reverses the order index rather than acting pointwise.
The gradient needs them, so that gap closes first:

**`expansion.Coefficient<α…>(l, m)`**, defined for every representable
component, resolving in turn:

| the component is | the coefficient is |
|---|---|
| stored | read directly |
| a permutation relative | `± ` the representative's |
| a reality relative | `sign · (−1)^{m+N_rep} conj(T^{rep}_{l,−m})` |
| in a real block, `m < 0` | `(−1)^m conj(stored_{l,|m|})`, the reduced storage |
| pinned as imaginary | `i` times the stored real field's |
| a vanishing orbit | zero |

This is worth having on its own: it is the spectral counterpart of phase 4's
component accessor, and it closes the asymmetry §15.4 recorded.

### 16.3 Decisions

**The result has no symmetry.** `∇₁` of a symmetric tensor is not symmetric in
the new slot against the old ones, and inferring a symmetry from an operator is
the same problem `Materialise` declined to solve.

**The result keeps the operand's reality.** `∇₁` of a real tensor is real, so
only the reduced set is computed and the rest follows.

**The `σ = 0` block is stored and zero.** A rank-`(q+1)` tensor whose `e₀` slot
vanishes is an ordinary tensor, and storing it as one keeps the result
composable with contraction, further gradients and the transform at a cost of
a third of the storage. Introducing a "tangential tensor" type to avoid that
would buy less than it complicates.

**Divergence and curl are not separate operators.** They are contractions of
the gradient (D&T C.6.2), and contraction is pointwise, so a caller evaluates
the gradient and contracts in the spatial domain with the machinery phase 3
already provides.
