// 04 -- Lazy evaluation
//
// Operators build expression nodes rather than fields. Nothing is computed
// until a value is asked for, so an intermediate is never materialised unless
// you ask for it.

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <iostream>

// Asked through a concept, since a bare requires-expression in a plain
// function makes GCC report the failure rather than answer it. See example 03.
template <typename F, typename E>
concept TimesAssignable = requires(F f, E e) { f *= e; };

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  auto grid = Grid(16, 2);

  auto u = SpinField<1, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::sin(theta) * std::cos(phi), 0.3};
  });
  auto v = SpinField<1, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::cos(theta), std::sin(phi)};
  });

  // An expression is a value of its own type, cheap to hold and to copy. This
  // allocates nothing and computes nothing.
  auto e = conj(u) * v + 2.0 * abs2(u);
  static_assert(decltype(e)::UpperIndex == 0);

  // It is evaluated point by point, on demand.
  std::cout << "e at one point " << (e[3, 4]) << "\n";

  // Materialise breaks the chain when a result is wanted more than once --
  // otherwise the tree is walked again at every read.
  auto stored = Materialise(e);
  std::cout << "materialised   " << (stored[3, 4]) << "\n\n";

  // Assignment evaluates into existing storage. The destination may appear on
  // the right: every node in this layer is pointwise *and index-preserving*,
  // reading its operands only at the point it is writing, so evaluating in
  // place needs no temporary. That is a theorem about the node set, not an
  // accident, and it is why no re-indexing view is allowed to sneak in here.
  auto w = SpinField<1, Grid>(grid, [](auto, auto) { return Complex{1.0, 0.0}; });
  const auto before = w[3, 4];
  w = w + conj(u) * v * w;
  std::cout << "in place: " << before << " -> " << (w[3, 4]) << "\n";

  // Compound assignment, constrained exactly as the binary form is: the
  // right-hand side must land at the destination's upper index, so a scalar
  // field may multiply in place and a spin-1 field may not.
  w *= abs2(v);
  static_assert(!TimesAssignable<decltype(w), decltype(v)>);
  static_assert(TimesAssignable<decltype(w), decltype(abs2(v))>);

  // Map applies an arbitrary callable pointwise, and lives at upper index
  // zero because that is the only place a nonlinear function of a component
  // is again a component. pow, exp and log are all this.
  auto damped = Map(abs2(u), [](auto x) { return std::exp(-x); });
  static_assert(decltype(damped)::UpperIndex == 0);
  std::cout << "exp(-|u|^2)    " << (damped[3, 4]) << "\n";
}
