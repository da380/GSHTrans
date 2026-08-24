# Code review: a working list

A file-by-file read of `GSHTrans/src`, done fresh from the code rather than
from the plans. The output is a list of things to examine or change, not the
changes themselves; each is numbered so it can be picked up, argued with or
dropped.

Three things were asked for beyond ordinary code quality: **Doxygen comments**,
**removal of refactoring residue**, and comments that say *what* and *why*
rather than narrating how the code came to be.

---

## Status

Most of the list has been worked through. What follows records the outcome so
the review stays usable as a record rather than reading as a list of things
still to do.

**Done.**

- **[R0.1]** decided as recommended, and applied: the measured facts and the
  constraints stay, the plan-step tags and the development history go. Every
  `[C..]`, `[D..]`, `[E..]`, `[I..]`, `[R..]`, `M..`, `P..`, `step X` and
  plan-section citation has gone from `GSHTrans/src`. References to the theory
  note and to `docs/gshtrans-reference.tex` stay, since those point at the
  mathematics; the theory note is now named as `docs/canonical-components.tex`
  once in each header that defers to it ([RX.4]).
- **[R0.2]** settled and complete. `docs/Doxyfile.in` is configured by CMake,
  `cmake --build build --target docs` builds the documentation, and
  `EXTRACT_ALL` is off with `WARN_IF_UNDOCUMENTED` on. **Every public entity
  in the library is documented**, the run is clean, and `WARN_AS_ERROR` is on
  with a CI leg to keep it that way. Private helpers are documented in the
  source and excluded from the output.
- **[RX.1]**, **[RX.2]** the include hygiene; **[R3.1]**–**[R3.3]** the dead
  tags; **[RX.3]** the phase vocabulary; **[R8.1]** the `GridBase` references;
  **[R12.2]** the 3-j narrative.
- **[RX.5]** the shared detail namespaces. The two `AnyRadial` predicates are
  one `HasRadialSlot` in `MultiIndex.h`, where it is a fact about the alphabet;
  `Expansion/BundleMaps.h` has its own `SpectralBundleDetails`; and `Omega`
  moves from `ContravariantDerivative.h` into `Eth.h`, so nothing defines into
  another header's namespace.
- **[RX.6]** the kernel asymmetry. `ForwardLoopKernel` and `InverseLoopKernel`
  are named private members beside the matrix pair, and the public entry
  points are 67 and 51 lines rather than 220 and 149.
- **[R6.1]**, **[R6.2]** the BLAS hazards, documented at the declarations.
  There is still no configure-time check, and there is no portable way to make
  one from inside the header: an ILP64 BLAS is a constraint on how the library
  is configured.
- **[R7.1]** the unused `<omp.h>`; **[R7.2]** the `Policies.h` comments;
  **[R2.2]** the `IsFastFFTSize` instance, now stated as the grid's `nPhi`.

**Answered rather than acted on.**

- **[R2.3]** `MinusOneToPower`'s default return type *is* taken — at seven call
  sites in the library and in the tests and examples, wherever the sign is
  wanted as an integer. The default stays.
- **[RX.6]**, second half: `RadialResample.h`'s `Resample` is long but already
  decomposed, into a per-line `run` and a per-piece `fit`. Extracting either
  would mean threading a dozen parameters through to gain nothing a reader
  needs. Left as it is.

**Still open, each needing a decision rather than an edit.**

- **[R1.1]** the version number. `1.0.0` is stale and only the value is; the
  mechanism works. Deliberately not chosen here.
- **[R7.3]** whether `Scheme` is renamed or the opening comment amended. The
  opening now says the tags are the exception and why, which may be enough.
- **[R7.4]** `Chunking::Count` taking `int copies` against `std::ptrdiff_t`
  everywhere else.
- **[R8.3]** `SphericalGrid.h` at 1841 lines. Smaller now that the entry points
  are, but still the largest file by a wide margin.
- **[R11.2]** the `Tuning.h` signatures, and **[R16.3]** the duplication
  between the two layered field types.

The clang-format question that used to sit here is settled: the tree conforms
to the checked-in `.clang-format`, and a CI leg keeps it there.

---

## 0. Two decisions to take first, because they affect every file

### [R0.1] What happens to the plan references — needs a decision

Measured across `GSHTrans/src`: **3607 comment lines in 12,853 total**, about
28%, carrying roughly **250 references to planning documents and decision
tags** — `core-plan.md`, `[C18]`, `F9`, `M6`, `section 22.1`, and so on.

Much of that is exactly the narrative to remove. But not all of it is alike,
and the distinction decides how mechanical this can be:

| kind | example | keep? |
|---|---|---|
| a measured fact justifying a design | "the table is 648 MB against 2.1 MB for a field, so copying one by accident is an out-of-memory event" | **keep** — it is the *why* |
| a constraint a reader must not break | "the recursion's output for one `(n, θ)` spans every order at once, so a per-order block cannot be had without keeping the table" | **keep** |
| what the code used to be | "this used to be held by value with defaulted copy" | **cut** |
| which plan step did it | "(core-plan.md F9)", "[C18]", "M6 measured" | **cut, or reduce to one pointer per file** |
| a measurement's history | "the first version measured worse and was taken out" | **cut** |

*My recommendation*, for a decision rather than an assumption: keep the facts
and the constraints, cut the history, and allow **at most one** plan pointer
per header — in the file-level comment, of the form "the design is recorded in
`docs/core-plan.md` §11" — rather than a tag on each paragraph. That keeps the
audit trail discoverable without the code narrating itself.

The plans stay as they are: they are the record, and they are where the
history belongs.

### [R0.2] Doxygen: coverage is one file in thirty-nine

`3j.h` has 58 `@brief` and full `@param`/`@details`. **Every other header has
none.** There is no `Doxyfile`, and no CI leg building the documentation.

So this is not a top-up, it is the whole job, and it needs a house style
settled before it starts or thirty-nine files will drift. Points to settle:

- `@file` blocks on every header, or only where the file's subject needs one?
- `@param`/`@return` everywhere, or only where the name is not self-evident?
  Most of this library's functions have long, explicit names.
- Are the *private* helpers documented? Much of the substance is private.
- Doxygen groups (`@defgroup`) for the five umbrellas, so the generated
  output has the same shape as the library?
- A `Doxyfile` and a CI leg, so that a malformed block fails rather than
  rendering wrongly and unnoticed. FFTWpp already does this and is worth
  copying from.

`3j.h` is the obvious model, since it is the one file already in the target
state.

---

## 1. `Version.h` (47 lines)

**[R1.1]** The version is `1.0.0` and has not moved while the layered fields,
tangential tensors, radial operators, interpolation, tuning and the 3-j
rewrite all landed. The file's own comment says the macros exist so a
dependent project can ask which facilities are present — which they cannot
usefully do at `1.0.0`. Decide a number and bump it; CMake already refuses to
configure if the two disagree, so the mechanism works and only the value is
stale.

**[R1.2]** Needs Doxygen: `@file`, and `@brief` on the four constants and the
three macros.

**[R1.3]** Minor: the comment lists which facilities arrived when. That is the
narrative form; the same point is made by "a dependent project may need to
know which facilities are present".

## 2. `Utility.h` (42 lines)

**[R2.1]** Needs Doxygen. Three small functions, all worth `@brief` and
`@param`; `IsFastFFTSize` also wants `@return`.

**[R2.2]** `IsFastFFTSize`'s comment cites "514 = 2 × 257 is the case that
matters here" without saying where *here* is. It is the grid's `nPhi` at
`lMax = 256`. Either say so or drop the instance.

**[R2.3]** `MinusOneToPower` defaults its return type to `std::ptrdiff_t`,
but every call site in the library asks for a floating-point type. Worth
checking whether the default is ever taken; if not, drop it and let the caller
say.

## 3. `Concepts.h` (99 lines)

**[R3.1] Dead tags — three of them.** `UpperIndexFirst` and `AngleFirst` occur
exactly once each in the whole repository, at their own definitions. Nothing
names them. Delete.

**[R3.2] `ScalarFunctionS2Expansion`** likewise: defined here, used nowhere,
in `src`, `tests` or `examples`. Delete, or find the caller it was written
for.

**[R3.3] `RowMajor` and `WignerStorage` are all but dead.** `RowMajor` appears
only at its definition and in the `WignerStorage` concept; nothing ever passes
it. `Wigner`'s `Storage` parameter therefore has exactly one reachable value,
`ColumnMajor`, and the concept exists to constrain a choice that is not a
choice. Either delete the axis — which removes a template parameter from a
heavily instantiated class — or write down what `RowMajor` would mean and what
would use it.

**[R3.4]** `#include <stdexcept>` is unused — nothing in the file throws.

**[R3.5]** Needs Doxygen throughout: the tag types especially, since their
names alone do not say what they select.

**[R3.6]** The tag types are grouped under a banner reading "Tag classes", but
they select four unrelated things: Wigner storage order, index ranges, angle
ranges and value kinds. Worth splitting the banner so a reader can see which
tags go with which axis.

## 4. `Views.h` (93 lines)

**[R4.1] The only file in the library with no prose at all.** Four class
templates, two banners, no explanation of what a `GSHView` is, what it views,
who owns the storage, or why the const and non-const versions are separate
types rather than one parameterised on constness. This is the file most in
need of the Doxygen pass.

**[R4.2] Missing includes.** `std::next` is used with no `<iterator>`, and
`std::ptrdiff_t` with no `<cstddef>`. It compiles only through `Indexing.h`;
this is the defect class the self-sufficiency test was added for, and the test
does not catch it because the transitive include is always present.

**[R4.3]** The banner comment style — `/*---- ... ----/` — differs from the
`//---- ... ----//` used everywhere else.

**[R4.4]** `GSHView::begin()`/`end()` are non-const, so a `const GSHView` is
not iterable while a `ConstGSHView` is. Probably an oversight rather than a
decision; worth deciding.

**[R4.5]** Worth documenting that indexing is unchecked, since these are the
lowest-level accessors in the library and the omission is deliberate.

## 5. `Indexing.h` (155 lines)

**[R5.1] The most opaque arithmetic in the library, and it has no
explanation.** `OffsetForDegree` is two closed forms — one for `All`, one for
`NonNegative` — each a branch on whether `mMax` exceeds `|n|`, each built from
triangular numbers. Nothing says what layout they encode, and a reader cannot
recover it from the formulas without deriving it. This is the single highest
value comment in the library to write.

**[R5.2]** `GSHIndices() = default` leaves `_lMax`, `_mMax` and `_n`
uninitialised. Three classes declare a `GSHIndices` member, and all three
initialise it in their constructors, so the default appears to exist only to
permit the member declaration. Either it is unnecessary — delete it — or the
members should be given initialisers, which costs nothing and removes a way to
read rubbish.

**[R5.3]** The bounds checks are `assert`-only. That matches the library's
stated policy, but the policy is stated in a planning document rather than
here, and this is the layer where a reader most needs to know.

**[R5.4]** Needs Doxygen throughout; `Indices()` in particular, which returns a
joined view of `(l, m)` pairs and whose type is unguessable.

## 6. `Blas.h` (96 lines)

**[R6.1] A real portability hazard, undocumented.** The Fortran symbols are
declared with `const int*` dimensions. An ILP64 BLAS — MKL's 64-bit-integer
build, which is a common choice for large problems — passes 64-bit integers,
and these declarations would then be silently wrong rather than failing to
link. Worth at minimum a comment, and possibly a configure-time check.

**[R6.2]** Related: the symbols are spelled with a trailing underscore, which
is the common convention but not universal. If a platform is ever added where
it differs, this is where it breaks.

**[R6.3]** The prose here is good and mostly *why* — the column-major
reversal, and why the transpose flag is the whole difference between the two
directions. It converts to Doxygen almost as it stands.

## 7. `Policies.h` (413 lines)

**[R7.1] `#include <omp.h>` is unused.** No `omp_` call and no pragma appears
in the file. Thirteen headers include `Policies.h`, so this drags OpenMP's
header into nearly the whole library for nothing — and it is the one header
clang cannot parse when it comes from GCC, which is the practical cost.

**[R7.2] The file is the clearest case for [R0.1].** Its comments are long and
carry most of the library's measurement history — `[C9]`, `[C10]`, `[C12]`,
`[C17]`, `P2`, `P8`, "section 10", the 2.3× and the 2.2× and the anchors that
produced the chunk formula. The *rules* in them are essential; the record of
which run produced which number is not.

**[R7.3]** `Scheme` is shaped unlike its four neighbours: a class with no
instances, whose named constructors return nested tag types. That is
deliberate and documented, but it means `Scheme` is not a policy *value* in
the way `Chunking` is, and the file's own opening comment describes all of
them as values. Worth either renaming it or amending the opening.

**[R7.4]** `Chunking::Count` takes `int copies`, while everything around it
uses `std::ptrdiff_t`. Minor inconsistency in a hot signature.

---

## Cross-cutting

### [RX.1] Twelve headers compile only through transitive includes

Each of these uses a standard facility whose header it does not include:

| header | missing |
|---|---|
| `Views.h` | `<cstddef>`, `<iterator>` |
| `Indexing.h` | `<cstddef>`, `<cmath>` |
| `SphericalGrid.h` | `<cstddef>`, `<complex>`, `<iterator>` |
| `WignerMatrices.h` | `<algorithm>`, `<cmath>` |
| `GaussLegendreGrid.h` | `<complex>` |
| `SpinField/SpinField.h` | `<complex>`, `<utility>` |
| `SpinField/SpinFieldView.h` | `<complex>`, `<utility>` |
| `Expansion/SpinExpansion.h` | `<cmath>` |
| `Expansion/Interpolate.h` | `<utility>` |
| `Tensor/TensorExpr.h` | `<stdexcept>` |
| `Layered/RadialDerivatives.h` | `<utility>` |
| `Layered/RadialResample.h` | `<utility>` |

This is the defect class the self-sufficiency test was added for, and **the
test does not catch it**: compiling a header alone still pulls the missing
facility in through one of the library's own headers. Worth fixing, and worth
asking whether the test can be strengthened — an include-what-you-use pass in
CI would catch the whole class at once.

### [RX.2] Nine unused includes

`Policies.h` `<omp.h>`; `SphericalGrid.h` `<numeric>`; `Utility.h`
`<initializer_list>`; `Expansion/IntrinsicDerivative.h` `<array>`;
`SpinField/SpinFieldNodes.h` `<string>`; `SpinField/SpinFieldOverloads.h`
`<functional>`; `Layered/RadialOperator.h`, `Layered/RadialResample.h` and
`Layered/RadialSplineDerivative.h` `<string>`.

### [RX.3] The code speaks in development phases

**27 references to "phase 1", "phase 2", "phase 4" across 8 headers** — "a
phase-1 node", "phase 4's second buffer", "what phase 2 left for phase 4".

These are the field-algebra plan's chapter numbers. A reader of the library has
no idea what phase 4 is, and the phases have no meaning once the work is done.
Each one names something real that has a proper name: *a phase-1 node* is **a
spin-weighted node**; *phase 4's second buffer* is **the real buffer a reality
reduction needs**; *what phase 2 left for phase 4* is **the reality reduction**.

This is the clearest single instance of the refactoring residue, and it is
mechanical to fix.

### [RX.4] "The theory note" is never named

13 references in 9 headers to "theory note section 6", "theory note
eq:reality". The document is `docs/canonical-components.tex`, and a reader has
no way to know that. Name it once per file, as `3j.h` names its references.

Unlike the plan references of [R0.1], **these should stay**: they point at the
authority on the mathematics, which the code genuinely defers to, and that is
a *why* rather than a history.

### [RX.5] Two headers share a detail namespace with unrelated contents

`Tensor/BundleMaps.h` and `Expansion/BundleMaps.h` both open
`GSHTrans::BundleDetails`. They are the same namespace, so the two files'
private helpers sit together, and a future name collision between them is a
confusing error rather than an impossible one. The two also carry
near-duplicate predicates — `AnyRadial<Alphas...>()` and
`AnyRadialIn<Indices>()` — which compute the same thing over different
spellings of a multi-index.

`EthDetails` is likewise opened by both `Expansion/Eth.h` and
`Expansion/ContravariantDerivative.h`. That one may well be deliberate reuse;
worth confirming and saying so.

Separately, the naming convention is `<Thing>Details` everywhere except
`SpinFieldOps` and `IndexRules`. Those two may be deliberate — they are not
really "details" — but the exception should be a decision.

### [RX.6] The loop kernel is inline; the matrix kernel is not

`SphericalGrid.h` has `ForwardMatrixKernel` and `InverseMatrixKernel` as named
private members, and the *loop* kernel written inline inside the public
`ForwardTransformation` and `InverseTransformation` — which are **220 and 149
lines**. The two kernels are peers and should read as peers. Extracting
`ForwardLoopKernel`/`InverseLoopKernel` would make the public entry points
short enough to see whole, and would make the two paths comparable in the
source as [C12] makes them comparable in measurement.

`Layered/RadialResample.h`'s `Resample` is 116 lines and is the other function
worth looking at.

---

## 8. `SphericalGrid.h` (1841 lines)

**[R8.1] Refactoring residue, three sites.** A banner still reads *"Methods
needed to inherit from GridBase"*; `CoefficientSizeFor` is documented as
*"named distinctly from GridBase::CoefficientSize, which it would otherwise
hide"*; and the class comment narrates the split. `GridBase` no longer exists.
The banner is simply wrong, and the `CoefficientSizeFor` note now describes a
hiding relationship between two members of *one* class — which is a different
question, and worth revisiting: can the two be one overload set?

**[R8.2]** See [RX.6]: the two longest functions in the library are here.

**[R8.3]** It is the largest file by a factor of two. Once the loop kernels are
extracted, consider whether the Fourier stages and the workspace machinery
want their own header — the file currently holds the public transform, two
kernels, two Fourier stages, the plan cache, the scratch buffers and the
validation.

**[R8.4]** Heavy Doxygen job, and the one where `@details` will earn its place:
the batch descriptors, the chunking rule and the threading contract are all
things a caller must understand and cannot guess.

## 9. `GaussLegendreGrid.h` (153 lines)

**[R9.1]** In good shape, being newly written. Needs Doxygen and one plan
reference removed (`[C26]`).

**[R9.2]** `Nodes()` builds a `GaussQuad::Quadrature1D` only to call
`Transform` on it and then copy the results out. Whether the mapping is worth
doing directly on the pair `GaussQuadrature` returns is a small simplification
to consider.

## 10. `Wigner.h` (689) and `WignerMatrices.h` (330)

**[R10.1]** `Wigner`'s `Storage` template parameter has one reachable value —
see [R3.3]. Removing it simplifies a heavily instantiated class.

**[R10.2]** `WignerMatrices.h` is missing `<algorithm>` and `<cmath>`
([RX.1]).

**[R10.3]** Both need Doxygen; the recursion in particular, where the seed
row, the boundary terms and the orthonormalisation each need `@details`.

**[R10.4]** These two hold the same values in two layouts and are checked
against each other by test. Worth a cross-reference in each file's `@file`
block so a reader finds the other.

## 11. `Tuning.h` (467 lines)

**[R11.1]** Newest file after `GaussLegendreGrid.h`; comments are already
why-shaped but carry `[C12]`, `[C17]`, `[C18]`, `[C19]`, `[C20]`, `[C22]` —
the densest plan-tag file in the library after `Policies.h`.

**[R11.2]** `TuneKernelLoopOnly` takes eight parameters, seven of them
forwarded unchanged. A small struct of grid parameters would make both it and
`TuneKernel` readable, and `TuneKernel`'s own nine-parameter signature is at
the edge of what a caller can use correctly.

**[R11.3]** `TunedChunking` and `TunedKernel` are near-parallel result structs
that do not share a shape. Worth deciding whether they should.

## 12. `3j.h` (806 lines)

**[R12.1] The Doxygen model for the rest of the library**, and already in the
target state.

**[R12.2] But its `@file` block narrates the development.** It explains at
length what the two superseded schemes were, where each failed, and that
Schulten–Gordon replaced them. Under the brief, that goes: what stays is what
the algorithm *is* and why it is stable, and one sentence that the layout
Woodhouse's `wig2.f` used is offered as a convention.

**[R12.3]** The `@file` block is attached to a `namespace GSHTrans {` opening
rather than the top of the file, and names the file `Wigner3j.hpp`, which is
not its name.

## 13. `SpinField/` (962 lines over five headers)

**[R13.1]** `SpinWeighted.h` is the best-commented file in the library and
needs the least work beyond Doxygen conversion.

**[R13.2]** `SpinField.h` opens by explaining that the class is *named*
`SpinField` rather than `CanonicalComponentField`. The substantive half — that
"canonical component" means something else for rank ≥ 2, so the name would be
wrong — is worth keeping; the half about a rename is not.

**[R13.3]** `SpinField.h` and `SpinFieldView.h` both miss `<complex>` and
`<utility>` ([RX.1]).

**[R13.4]** `SpinFieldOverloads.h` has an unused `<functional>`;
`SpinFieldNodes.h` an unused `<string>`.

**[R13.5]** `SpinFieldOps` breaks the `<Thing>Details` naming convention
([RX.5]).

## 14. `Tensor/` (2114 lines over five headers)

**[R14.1]** Comments are strong and mostly *why*. The work here is Doxygen,
the phase vocabulary ([RX.3], concentrated in `Orbits.h` and `TensorField.h`),
and the `BundleDetails` duplication ([RX.5]).

**[R14.2]** `TensorExpr.h` is missing `<stdexcept>`.

**[R14.3]** `TensorExpr.h` opens `TensorDetails` twice, at lines 73 and 567.
Probably fine, but worth one namespace block or a note saying why two.

**[R14.4]** `MultiIndex.h`'s `SymmetryDetails` and `TensorExpr.h`'s
`TensorDetails` both hold symmetry machinery. Worth checking the split is
where a reader would expect it.

## 15. `Expansion/` (1642 lines over seven headers)

**[R15.1]** `Interpolate.h` is missing `<utility>`; `SpinExpansion.h` is
missing `<cmath>`; `IntrinsicDerivative.h` has an unused `<array>`.

**[R15.2]** `SpinExpansion.h` defines `SpinExpansionBase` *and* `SpinExpansion`
*and* view types. The base's role — and why the split exists — is not stated,
and this is the only place in the library where such a base survives.

**[R15.3]** `Interpolate.h` carries `[I1]` to `[I9]` throughout ([R0.1]).

**[R15.4]** `Eth.h` and `ContravariantDerivative.h` both open `EthDetails`
([RX.5]).

## 16. `Layered/` (2322 lines over nine headers)

**[R16.1]** `RadialOperator.h`, `RadialResample.h` and
`RadialSplineDerivative.h` each carry an unused `<string>`;
`RadialDerivatives.h` and `RadialResample.h` are missing `<utility>`.

**[R16.2]** `Resample` is 116 lines ([RX.6]).

**[R16.3]** `LayeredSpinField.h` exposes fourteen accessors and
`LayeredTensorField.h` six, with `FieldSize`, `MaxDegree`, `NumberOfRadii` and
`RadiusIndices` duplicated between them. Whether a shared base or a concept
would remove that is worth asking — the `LayeredStack` concept already exists
and may be the place.

**[R16.4]** The directory is the largest after `SphericalGrid.h` and the least
uniformly documented; `RadialMajor.h` and `RadialOperator.h` are well
explained, the two layered field types much less so.

---

## Suggested order

The list is deliberately not sequenced by file, because the cheap mechanical
passes cut across all of them and are worth doing first — they are the ones
that can be checked by the compiler rather than by reading.

1. **[R0.1] and [R0.2] decided.** Everything else depends on both.
2. **The include hygiene**, [RX.1] and [RX.2]. Mechanical, compiler-checkable,
   and it is worth asking whether include-what-you-use in CI replaces the
   self-sufficiency test — which does not catch this class.
3. **The dead code**: [R3.1], [R3.2], [R3.3]. Deletions, each independently
   verifiable by the suite.
4. **The residue**: [RX.3] the phase vocabulary, [R8.1] the `GridBase`
   references, [R12.2] the 3-j narrative. Comment-only, so the suite cannot
   check them — read them.
5. **The structural items**: [RX.6] the kernel asymmetry, [R8.3] the file
   size, [R11.2] the signatures, [R16.3] the layered duplication. Each is a
   judgement call and each deserves its own decision.
6. **Doxygen**, last, because writing it before the above means writing some of
   it twice.

## What this review is not

A first pass. Every file was opened and its structure and commentary read;
the larger ones — `SphericalGrid.h`, `TensorField.h`, `TensorExpr.h`,
`LayeredTensorField.h` — were surveyed rather than read line by line, so their
entries are lighter than their size warrants and should not be taken as
"little to do here". The Status section above records what has since been
acted on.

No correctness defect was found. That is worth saying plainly: the findings
are hygiene, documentation, dead code and residue, not bugs.
