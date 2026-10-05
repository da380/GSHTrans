// 03 -- Spin weight
//
// A field on the sphere carries an upper index N -- a spin weight -- saying
// how it responds to a rotation of the local frame. The library tracks N in
// the type, so the index arithmetic of an expression is checked where it is
// written.
//
// The mathematics, briefly. At each point take the complex tangent vectors
// e_+- = -+(theta_hat +- i phi_hat) / sqrt(2) and e_0 = r_hat. Rotating the
// tangent frame through psi about r_hat sends e_+- to exp(-+ i psi) e_+-, and
// a component of a tensor in this basis picks up a phase exp(-i N psi),
// where N is the sum of its indices. That N is the field's upper index (spin
// weight): a scalar has N = 0, the e_+- components of a vector have N = +-1,
// and so on (example 08). A field of upper index N expands in the
// generalised spherical harmonics Y^N_lm, which is what the transform of
// example 06 computes. See docs/gshtrans-reference.tex, sections
// "Conventions" and "What the type system enforces".
//
// What this shows
//   The index rules of the field algebra read off the types, and the
//   operations the type system refuses.
//
// Assumes
//   Example 02 (a SpinField, Integrate).
//
// Introduced
//   Complex-valued SpinField at N != 0, conj, abs2, the UpperIndex member,
//   the SpinWeighted concept, testing for an operation with a concept.
//
// Output
//   An inner product and a norm of spin-2 fields; everything else in the
//   file is checked at compile time.

#include <GSHTrans/GSHTrans.hpp>
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

  // nMax = 2 so that the grid can store fields at upper index -2..2.
  auto grid = Grid(16, 2);

  // Fields at N != 0 are complex-valued; ComplexValued is the default. The
  // functions here are arbitrary, chosen only to have something to compute.
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
  // A product adds upper indices, since the phases multiply, and abs2 (|u|^2,
  // that is conj(u) * u) lands at zero. Every node -- field, view or
  // expression -- has a static UpperIndex.
  static_assert(decltype(u)::UpperIndex == 2);
  static_assert(decltype(conj(u))::UpperIndex == -2);
  static_assert(decltype(u * v)::UpperIndex == 4);
  static_assert(decltype(u * s)::UpperIndex == 2);
  static_assert(decltype(abs2(u))::UpperIndex == 0);

  // u * v sits at N = 4, beyond this grid's nMax. That is fine for an
  // expression, which stores nothing; materialising or transforming it would
  // need a grid carrying N = 4, and a field built on one that does not
  // throws std::invalid_argument.

  // So conj(u) * v lands at zero, and only there can it be integrated: the
  // integral over the sphere of a field of nonzero upper index vanishes
  // identically, so the library will not let you ask for it. Integrate of
  // conj(u) * v is the L2 inner product <u, v>.
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
  // SpinField<2, Grid, RealValued> is rejected by a static_assert inside the
  // class -- a hard error, so it cannot be probed like the cases above -- and
  // the SpinWeighted concept, which every node satisfies, requires the same.
  // Only the lawful case can be asserted here.
  static_assert(SpinWeighted<SpinField<0, Grid, RealValued>>);

  std::cout << "u + s does not compile; conj(u) * v does.\n";
}
