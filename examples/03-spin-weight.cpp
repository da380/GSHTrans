// 03 -- Spin weight
//
// A field on the sphere carries an upper index N -- a spin weight -- saying
// how it responds to a rotation of the local frame. The library tracks N in
// the type, so the index arithmetic of an expression is checked where it is
// written.

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <iostream>

// The library's constraints are requires-clauses on the operators, so an
// unlawful combination is an overload-resolution failure rather than an error
// inside a node -- which means it can be *asked about*. The question has to be
// put through a concept: inside a plain function GCC reports "no match for
// operator+" eagerly instead of answering false.
template <typename L, typename R>
concept Addable = requires(L l, R r) { l + r; };

template <typename L, typename R>
concept Divisible = requires(L l, R r) { l / r; };

template <typename A>
concept HasRealPart = requires(A a) { GSHTrans::real(a); };

template <typename A>
concept Integrable = requires(A a) { GSHTrans::Integrate(a); };

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  auto grid = Grid(16, 2);

  auto u = SpinField<2, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::sin(theta) * std::cos(phi), std::sin(2 * phi)};
  });
  auto v = SpinField<2, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::cos(theta), 0.5 * std::sin(phi)};
  });
  auto s = SpinField<0, Grid>(grid, [](auto theta, auto) {
    return Complex{1 + 0.25 * std::cos(theta)};
  });

  // The rules, read off the types. Conjugation *reverses* the upper index: a
  // quantity carrying exp(-i N psi) has a conjugate carrying exp(+i N psi).
  static_assert(decltype(u)::UpperIndex == 2);
  static_assert(decltype(conj(u))::UpperIndex == -2);
  static_assert(decltype(u * v)::UpperIndex == 4);
  static_assert(decltype(u * s)::UpperIndex == 2);
  static_assert(decltype(abs2(u))::UpperIndex == 0);

  // So conj(u) * v lands at zero, and only there can it be integrated: the
  // integral over the sphere of a field of nonzero upper index vanishes
  // identically, so the library will not let you ask for it.
  std::cout << "<u, v>  = " << Integrate(conj(u) * v) << "\n"
            << "||u||   = " << std::sqrt(Integrate(abs2(u))) << "\n\n";

  // What the type system refuses. Each of these is an overload-resolution
  // failure at the point of writing, not a runtime check.
  using Spin2 = decltype(u);
  using Scalar = decltype(s);

  static_assert(!Addable<Spin2, Scalar>);  // unequal upper indices
  static_assert(Addable<Spin2, Spin2>);    // equal ones are fine
  static_assert(!HasRealPart<Spin2>);      // covariant only at N = 0
  static_assert(HasRealPart<Scalar>);
  static_assert(!Divisible<Scalar, Spin2>);  // a divisor must be at N = 0
  static_assert(Divisible<Spin2, Scalar>);
  static_assert(!Integrable<Spin2>);  // vanishes identically
  static_assert(Integrable<Scalar>);

  // And the constraint that cannot be written at all: a real-valued field at
  // nonzero upper index. Real-valuedness is not preserved by the frame
  // rotation, so it is not a property any component of any tensor can have.
  static_assert(SpinWeighted<SpinField<0, Grid, RealValued>>);

  std::cout << "u + s does not compile; conj(u) * v does.\n";
}
