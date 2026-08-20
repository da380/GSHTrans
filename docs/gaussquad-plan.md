# GaussQuad: notes and suggested work

A hand-over document. Everything below was measured or read from the source at
`da380/GaussQuad@main` as fetched on 2026-08-20; nothing is inferred from the
interface alone. Recommendations are ordered by value, not by effort.

## What it is, and what depends on it

462 lines, header-only: `JacobiPolynomial` and `LegendrePolynomial`
(evaluation, derivative, roots, and Gauss, Gauss–Radau and Gauss–Lobatto
quadrature), plus a `Quadrature1D` holder with `X`, `W`, `Points`, `Weights`,
`Integrate` and `Transform`.

GSHTrans uses six entry points and no more:
`LegendrePolynomial<Real>{}.GaussQuadrature(n)`, then `Quadrature1D`'s
`Transform`, `Points`, `Weights`, `X` and `W`. The Gauss–Legendre nodes are the
colatitudes of every grid in the library, so their accuracy is the accuracy of
every transform.

**It matters at both ends of the target applications.** The angular side needs
Gauss–Legendre; the radial side of the spectral-element codes needs
Gauss–Lobatto–Legendre, which this library also provides. Work here serves both.

## What was found

### There are two algorithms, switched at n = 100, and they do not agree

`GaussQuadrature(n)` takes the Golub–Welsch path for `n ≤ 100` — build the
Jacobi matrix, call `Eigen::SelfAdjointEigenSolver` — and a Newton path for
`n > 100`, finding the roots of `P_n` with an asymptotic initial guess and a
Maehly deflation term, then forming the weights from `P'_n`.

Measured, either side of the switch:

| `n` | path | `Σw − 2` |
|---|---|---|
| 95 | matrix | `3.6×10⁻¹⁵` |
| 100 | matrix | `8.9×10⁻¹⁶` |
| **101** | **Newton** | **`1.2×10⁻¹³`** |
| 200 | Newton | `−1.1×10⁻¹³` |

Two orders of magnitude, at a step change in `n` of one. Nothing in the test
suite covers the crossover, so it has never been visible.

### The accuracy then degrades with degree

| `n` | 65 | 257 | 1025 | 2049 |
|---|---|---|---|---|
| `Σw − 2` | `1.1×10⁻¹⁵` | `5.7×10⁻¹³` | `3.5×10⁻¹³` | `4.0×10⁻¹²` |
| build time | 0.4 ms | 8.2 ms | 145 ms | 596 ms |

The cause is not subtle: the Newton path evaluates `P_n` and `P'_n` by
recurrence at every iteration, and that loses precision at roughly `O(nε)`.
`4×10⁻¹²` at `n = 2049` is about what `nε` predicts. Modern asymptotic methods
never evaluate the polynomial at all and hold near machine precision at any
degree.

### The cost is quadratic

`O(n²)`, from the ratios above. That is invisible at the sizes GSHTrans uses
today — 8 ms against a 648 MB Wigner table at `lMax = 256` — but it is the
*dominant* construction cost for a grid that generates its Wigner values rather
than storing them, and it is 0.6 s at `lMax = 2048`.

### A latent defect in the root finder, which is not the accuracy problem

In `Zeros`:

```cpp
auto sum = 0;                                   // deduces int
for (auto i = 0; i < k; i++)
  sum += static_cast<Real>(1) / (r - zeros[i]); // truncated on every +=
```

The Maehly deflation sum is accumulated in **integer** arithmetic. It should
be `Real sum = 0;`.

**It is worth reporting honestly that fixing this changes nothing measurable.**
I patched a local copy and re-ran the table above: `1.17×10⁻¹³` before,
`1.17×10⁻¹³` after. The reason is that the deflation enters as `der − sum·fun`,
and at convergence `fun → 0`, so the term vanishes exactly where it would have
mattered. What the defect costs is robustness rather than accuracy — deflation
exists to stop Newton falling back onto an already-found root, and truncated to
an integer it largely does not do that. It should still be fixed; it is one
word.

### Eigen is a dependency of this library alone

`OrthogonalPolynomial.h` includes `Eigen/Cholesky`, `Eigen/Core` and
`Eigen/Eigenvalues`, and the CMake fetches Eigen and links `Eigen3::Eigen`.
GSHTrans includes no Eigen header anywhere; it inherits the dependency entirely
through GaussQuad. Eigen is used for the `n ≤ 100` Gauss path and for
Gauss–Radau and Gauss–Lobatto, which have no Newton alternative and are always
solved as matrix problems.

### There are no install rules

No `install()` and no `export()` in `CMakeLists.txt`, so `find_package(GaussQuad)`
cannot work. This blocks more than it looks: a downstream `INTERFACE` target
can only be exported if everything it links is exported or imported, so
**GSHTrans cannot ship a CMake package while this is true**, and its
consumption is limited to `add_subdirectory` and `FetchContent`.

## Suggested work

**1. Add install and export rules.** Smallest change here, and it unblocks
downstream packaging. Twenty lines: `install(TARGETS … EXPORT …)`, an
`install(DIRECTORY GaussQuad …)`, an exported target set under a `GaussQuad::`
namespace, and a generated config file. Worth doing even if nothing else on
this list happens.

**2. `Real sum = 0`.** One word, and the deflation starts working.

**3. Replace the quadrature algorithms with an asymptotic method.** This is the
substantial item, and it is three improvements at once:

- `O(n²)` becomes `O(n)`, with `O(1)` per node;
- the accuracy stops depending on `n`, holding near `10⁻¹⁶` at any degree;
- **the Eigen dependency disappears**, because the method uses no linear
  algebra — which would leave the whole set header-only and small.

Bogaert (2014) is the reference: explicit asymptotic expansions for the nodes
and weights of Gauss–Legendre, accurate to near machine precision from small
`n` upward, and the basis of most modern implementations. Public-domain
implementations exist and are a few hundred lines. Keeping the `Quadrature1D`
interface exactly means nothing downstream changes.

The honest caveat: this covers Gauss–Legendre. Gauss–Radau and Gauss–Lobatto,
and Jacobi weights generally, do not have such tidy expansions, so the matrix
method would have to stay for those — and with it Eigen, unless they are given
a Newton treatment of their own. Since Lobatto is used at small orders in the
radial element bases, where the matrix method is both fast and accurate, the
pragmatic split is: asymptotic for Gauss–Legendre at any `n`, matrix for the
rest at the small `n` they are actually used at, and a documented note saying
so.

**4. Test what has never been tested.** The failures above were all found by two
identities that need no reference implementation:

- `Σᵢ wᵢ = 2` on `[-1, 1]`;
- `Σᵢ wᵢ P_l(xᵢ)² = 2/(2l+1)` at the highest degree the rule integrates
  exactly.

Run them at `n = 10, 95, 100, 101, 105, 1000, 10000`, with the crossover
deliberately straddled, and for Gauss, Radau and Lobatto alike. Had these been
present, the two-orders-of-magnitude step at `n = 101` would have been caught
when it was introduced.

**5. Consider whether the two paths should exist at all.** A single algorithm
that is accurate and fast at every `n` removes a discontinuity, a branch, and
an entire dependency. If the asymptotic route is taken for Gauss–Legendre, the
`n ≤ 100` special case has no remaining purpose.
