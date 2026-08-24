// 13 -- Raising and lowering
//
// eth and eth-bar change the upper index by one. They are the reason the
// library has two representations rather than one: everything in the field
// algebra is local in (theta, phi), and these are local in (l, m), so "the
// gradient of a product" is necessarily evaluated in both.

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <iostream>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  using Int = std::ptrdiff_t;

  constexpr auto lMax = Int{12};
  auto grid = Grid(lMax, 2);

  // A scalar field, expanded.
  auto f = SpinField<0, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::cos(theta) + std::sin(theta) * std::cos(phi), 0.0};
  });
  auto e = Expand(f, lMax);

  // Raising takes it to upper index +1, lowering to -1. Those two are the
  // canonical components of the surface gradient.
  auto raised = Raise(e);
  auto lowered = Lower(e);
  static_assert(decltype(raised)::UpperIndex == 1);
  static_assert(decltype(lowered)::UpperIndex == -1);

  // The degree ranges look after themselves. Raising from N loses degree |N|,
  // and the factor there is sqrt(0), so nothing is silently dropped.
  std::cout << "f       degrees " << e.MinDegree() << " .. " << e.MaxDegree()
            << "\n"
            << "eth f   degrees " << raised.MinDegree() << " .. "
            << raised.MaxDegree() << "\n\n";

  // Applying both gives the surface Laplacian, whose eigenvalue on Y_{lm} is
  // -l(l+1). This is the identity that pins the relative sign of the two
  // operators and both of their magnitudes.
  auto laplacian = Lower(Raise(e));
  auto worst = Real{0};
  for (auto l : e.Degrees()) {
    for (auto m : e.Orders(l)) {
      const auto expected = -static_cast<Real>(l * (l + 1)) * Complex{e[l, m]};
      worst = std::max(worst, std::abs(Complex{laplacian[l, m]} - expected));
    }
  }
  std::cout << "eth-bar eth f against -l(l+1) f: " << worst << "\n";

  // The other identity that survives the convention: [eth, eth-bar] = -2N.
  // The *overall* sign of the pair is a convention this library cannot settle
  // from the inside -- both identities are unchanged if you flip both -- and
  // it is implemented as the theory note states it.
  auto g = SpinExpansion<1, Grid>(grid, lMax);
  for (auto l : g.Degrees()) {
    for (auto m : g.Orders(l)) {
      g[l, m] = Complex{std::cos(0.3 * l + m), std::sin(0.2 * l - m)};
    }
  }
  worst = 0;
  for (auto l : g.Degrees()) {
    for (auto m : g.Orders(l)) {
      const auto commutator =
          Complex{Raise(Lower(g))[l, m]} - Complex{Lower(Raise(g))[l, m]};
      worst = std::max(worst, std::abs(commutator + 2.0 * Complex{g[l, m]}));
    }
  }
  std::cout << "[eth, eth-bar] against -2N:      " << worst << "\n\n";

  // Back to the spatial domain: the gradient's components are ordinary spin
  // fields, at upper index +1 and -1.
  auto gradientPlus = Evaluate(raised);
  static_assert(decltype(gradientPlus)::UpperIndex == 1);
  std::cout << "the raised field is a spin-1 field of " << gradientPlus.Size()
            << " samples\n";
}
