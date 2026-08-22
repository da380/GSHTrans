# GaussQuad 1.1.1 — what GSHTrans needs to know

A hand-over note on the state of `da380/GaussQuad` as of 2026-08-21, written
for whoever picks up the GSHTrans side. GSHTrans depends on GaussQuad for the
Gauss-Legendre nodes that are the colatitudes of every grid, and for the
Gauss-Lobatto-Legendre rules the spectral-element codes use on the radial side.
Everything quoted below was measured rather than inferred.

**Short version: nothing in GSHTrans has to change, and the thing GaussQuad
was blocking — GSHTrans shipping a CMake package of its own — is now
unblocked.**

---

## 1. Where the code is

| ref | what it is |
|---|---|
| `v1.0.1` | the library as it was before the refactor, tagged as a record |
| `main` | the refactored library, version `1.1.0` |
| `develop` | `1.1.1`: six further commits, pending release |

Pin `v1.0.1` if anything ever needs the old node values back.

---

## 2. Compatibility with GSHTrans

GSHTrans uses six entry points, and **all six are unchanged in signature and
semantics**:

- `LegendrePolynomial<Real>::GaussQuadrature(n)`
- `Quadrature1D<Real>::Transform(f, df)`
- `Quadrature1D<Real>::Points()`, `Weights()`, `X(i)`, `W(i)`

Every addition described below is additive. The install prefix was built and a
consumer compiled against it calling all six, resolving GaussQuad only through
`find_package`.

Two things to be aware of rather than to act on:

- **Node values have moved in the last few ulps.** The `n ≤ 100` /`n > 100`
  algorithm split is gone and one algorithm now serves every degree, so the
  nodes are not bit-identical to `v1.0.1`. David checked the dependent codes
  and none is sensitive to node values at that level. That, together with
  `v1.0.1` remaining available to pin, is why the version is `1.1.x` rather
  than `2.0.0`.

- **`Transform` keeps its Jacobian convention**, and is deliberately *not*
  deprecated, because GSHTrans calls it and that would fill its build with
  warnings. The convention is easy to get wrong and now documented in place:
  the points are mapped first and `df` is evaluated at the **mapped** points,
  so `df` must be a function of the new variable. For affine maps `df` is
  constant and the distinction does not arise — which is the GSHTrans case.

---

## 3. The CMake package — the part that was blocking GSHTrans

`find_package(GaussQuad)` now works:

```cmake
find_package(GaussQuad REQUIRED)
target_link_libraries(your_target PRIVATE GaussQuad::GaussQuad)
```

Previously there were no `install()` rules at all, and `${INCLUDE_INSTALL_DIR}`
was never defined anywhere, so the `INSTALL_INTERFACE` was empty regardless.
That is what stopped GSHTrans from being given install and export rules of its
own — a package cannot export a target linking a dependency that has no
package.

Also worth knowing on the GSHTrans side:

- **`GaussQuad::GaussQuad` alias** exists, so `add_subdirectory` and
  `find_package` can be used interchangeably.
- **`CMAKE_CXX_STANDARD` no longer leaks.** The standard is set with
  `target_compile_features(... cxx_std_23)` on the target, so adding GaussQuad
  as a subdirectory no longer changes the parent project's standard.
- **Developer targets default off when not top-level**, so
  `add_subdirectory(GaussQuad)` builds no tests or examples. Force with
  `-DGAUSSQUAD_BUILD_TESTS=OFF -DGAUSSQUAD_BUILD_EXAMPLES=OFF` to be
  explicit about it.
- **CI** covers g++-13 and g++-14 in Debug and Release, plus a job that
  installs to a prefix and builds a consumer through `find_package`, so a
  broken export set fails in GaussQuad rather than in GSHTrans.

**Next piece of work: give GSHTrans its own install and export rules.** That
is now purely a GSHTrans-side task.

---

## 4. No Eigen anywhere

The Eigen dependency is gone entirely — not just from Gauss–Legendre, which is
all that had been expected to be possible. `NumericConcepts` is the only
dependency, and the installed `GaussQuadTargets.cmake` links that and nothing
else.

The `n > 100` Newton path and the `n ≤ 100` Golub–Welsch path were both
replaced by one implicit-QL eigensolver on the Jacobi matrix that accumulates
only the first eigenvector row — `O(n²)` time and `O(n)` storage, against
`O(n³)` and `O(n²)`. The two `llt().solve()` calls in Radau/Lobatto became a
15-line Thomas solve.

If GSHTrans wants Eigen it must ask for it itself; it will no longer arrive
transitively.

---

## 5. What GSHTrans might actually want to use

Nothing here is required. Listed because it is directly relevant to what
GSHTrans does with quadrature.

### Faster construction, and accuracy that no longer depends on `n`

| `n` | build, before | after | `\|Σw − 2\|` before | after |
|---|---|---|---|---|
| 1025 | 152 ms | 23 ms | `3.4e-13` | `3.0e-15` |
| 4097 | 2620 ms | 305 ms | `3.3e-12` | `2.7e-15` |

The weight-sum error is now flat in `n` rather than growing with it. This is
structural, not a tuning: the accumulated eigenvector row is only ever acted on
by plane rotations, which preserve its norm. `long double` came out better than
before as well.

For the Gauss–Legendre grids this means the colatitudes of a high-degree grid
cost roughly a tenth of what they did.

### Exact endpoints on the Gauss–Lobatto–Legendre rules

Previously the Radau and Lobatto endpoints were *not* exactly `±1` — at
`n = 20`, `x_max` came back as `1.0000000000000036`, i.e. outside the interval.
They are now assigned exactly.

**This is the one that matters for the radial SEM bases.** A shared element
boundary is now the same number from both sides, and anything evaluating
`sqrt(1-x²)` or `log(1-x)` at a node no longer risks a NaN from a node a few
ulps past the endpoint.

Symmetry is also exact now when `alpha == beta`: `x_i = −x_{n+1-i}` and
`w_i = w_{n+1-i}` hold bitwise, so downstream code may rely on it.

### Rules placed on an interval directly

```cpp
auto q = GaussQuad::GaussLegendreQuadrature1D<double>(n, a, b);
auto r = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(n, a, b);
```

Equivalent to building the rule and mapping it, but with no Jacobian to get
right. For the radial element rules this is likely to be a straight
simplification of whatever currently calls `Transform`. There is also a
non-mutating `Transformed(f, df)` and a `MappedTo(a, b)`.

### An `O(n)` Gauss–Legendre algorithm, if construction cost ever matters

```cpp
auto q = GaussQuad::GaussLegendreQuadrature1D<double>(
    n, GaussQuad::Method::GlaserLiuRokhlin);
```

0.15 ms against 19 ms at `n = 1025`, and 2.3 ms against 4.2 s at `n = 16385`.
The nodes agree with the default to about `1e-15`, but the weights drift as
`O(nε)` where the default's stay flat — so at large `n` it is the faster rule,
not the better one. `GolubWelsch` remains the default and there is deliberately
no automatic crossover. Probably not worth it for GSHTrans unless grid setup
ever shows up in a profile.

### Other rules now available

Gauss–Laguerre, Gauss–Radau–Laguerre (node fixed at the origin) and
Gauss–Hermite, alongside the Jacobi, Legendre and Chebyshev families, each with
Gauss, Radau and Lobatto variants.

---

## 6. Behaviour changes to be aware of

- **Bad arguments now throw `std::invalid_argument`** where they used to
  `assert`. Since everything downstream is compiled with `NDEBUG`, a bad degree
  or a non-integrable weight was previously silently undefined. Anything that
  *builds* a rule validates; the inner-loop evaluation functions keep their
  assertions.

- **Every Gauss–Chebyshev rule used to return silent `NaN`** — two independent
  removable singularities at `alpha + beta = -1`. Both are fixed. If GSHTrans
  has ever avoided the Chebyshev rules, this is why, and it no longer applies.

- **`ChebyshevPolynomial::operator()` and `Derivative` never compiled.** Fixed.

- **The 2-point Lobatto rule** (the trapezoid rule) was rejected by an
  over-strict `assert(n > 2)`. It now works.

- **`Integrate` is now `const`**, so a `const Quadrature1D` is usable.

- **The `Integrable` concept was relaxed to convertibility.** It required
  `f*w` and `f+f` to be exactly the value type, which rejected every
  expression-template integrand — an Eigen vector returns a proxy from those
  operators. If GSHTrans ever wanted to integrate a vector-valued function and
  could not, that was why.

- **`Integrate(f)` is `∫ w(x) f(x) dx`**, not `∫ f(x) dx`. It is the plain
  integral only for the Legendre rules, which are the ones GSHTrans uses. This
  was previously undocumented and is a trap for Chebyshev and Jacobi users.

- **`Method::GlaserLiuRokhlin` throws** rather than returning an infinite
  weight when `n` exceeds what the precision can represent — in `float`, around
  `n = 10⁴`. `GolubWelsch` has no such limit. Only relevant if the `O(n)`
  algorithm is used at all.

---

## 7. Testing and documentation

The suite is 72 tests over `float`, `double` and `long double`, up from 12.
They are identity-based, so they need no reference implementation: weight sums,
orthogonality at the exact highest degree each rule reaches, strict
monotonicity, positivity, finiteness, exact endpoints, exact symmetry, closed
forms, and cross-agreement between two independent Gauss–Legendre algorithms.

They were validated by being run against the pre-refactor code, where they
catch every bug listed in §6.

- `README.md` — API and build instructions.
- `docs/algorithms.md` — the mathematics, the measured accuracy and cost, and
  the reasoning behind the implementation choices.
- Doxygen reference pages: configure with `-DGAUSSQUAD_BUILD_DOCS=ON` and build
  the `docs` target.

---

## 8. Known open items

None of these block GSHTrans.

1. **CI covers only g++-13 and g++-14**, because no clang is installed on the
   GaussQuad development machine, and an unverified matrix leg would fail on
   first push.
   Adding `clang++-18` is a one-line change once it can be tested. The
   `Nodes()` view is already guarded on `__cpp_lib_ranges_zip`, so a standard
   library without `std::views::zip` degrades rather than breaks.
2. **`Interpolation` is no longer used at all.** It was fetched at
   `GIT_TAG main` by the test target, which was the last unpinned dependency
   and the route by which Eigen arrived. The tests now carry their own random
   polynomial and fetch only GoogleTest.
