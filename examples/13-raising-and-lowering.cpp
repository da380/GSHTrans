// 13 -- Raising and lowering
//
// What this shows. The spin-raising and spin-lowering operators eth and
// eth-bar, which take a field of upper index N to one of upper index N + 1 and
// N - 1. On the basis functions they are
//
//     eth     Y^N_{lm} = -sqrt((l - N)(l + N + 1)) Y^{N+1}_{lm}
//     eth-bar Y^N_{lm} = +sqrt((l + N)(l - N + 1)) Y^{N-1}_{lm}
//
// so in the spectral domain each is a multiplication of the coefficient at
// (l, m) by a factor depending on l and N, with no coupling between different
// (l, m). That is why they live on expansions (example 12) and not in the
// pointwise field algebra: everything in the field algebra is local in
// (theta, phi), these are local in (l, m), and "the gradient of a product"
// therefore needs both representations -- the product in one, the derivative
// in the other. docs/gshtrans-reference.tex, "Raising and lowering", has the
// details and the sign convention.
//
// Read first. 03 (upper indices), 06 (the transform), 12 (expansions).
//
// Introduced. Raise(e) and Lower(e), which return SpinExpansions at N + 1
// and N - 1. The two identities that pin their normalisation: eth-bar eth is
// the surface Laplacian on a scalar, and [eth, eth-bar] = -2N.
//
// The output. The degree ranges before and after raising; both identities
// holding to rounding; and the size of the raised field once evaluated back
// to the grid.
//
// What this does not show. Raise and Lower act on one spin-weighted function.
// Applied component by component to a tensor they do not give the tensor's
// gradient; example 15 explains why and what to use instead.

#include <GSHTrans/GSHTrans.hpp>
#include <cmath>
#include <complex>
#include <iostream>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  using Int = std::ptrdiff_t;

  // Upper indices up to |N| = 2 are needed below: the scalar at 0, its raised
  // and lowered forms at +-1, and the spin-1 field g and its raised form at 2.
  constexpr auto lMax = Int{12};
  auto grid = Grid(lMax, 2);

  // A smooth scalar (N = 0) field, expanded. cos(theta) and
  // sin(theta) cos(phi) are degree-1 harmonics, so the expansion has content
  // at l = 1 only.
  auto f = SpinField<0, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::cos(theta) + std::sin(theta) * std::cos(phi), 0.0};
  });
  auto e = Expand(f, lMax);

  // Raise is eth, taking upper index 0 to +1; Lower is eth-bar, taking it to
  // -1. For a scalar these two are, up to a factor of -1/sqrt(2) and
  // +1/sqrt(2) respectively, the +1 and -1 canonical components of the surface
  // gradient (example 15 makes the factor explicit). The new upper index is
  // part of the result's type.
  auto raised = Raise(e);
  auto lowered = Lower(e);
  static_assert(decltype(raised)::UpperIndex == 1);
  static_assert(decltype(lowered)::UpperIndex == -1);

  // The degree ranges look after themselves. A field at upper index 1 starts
  // at l = 1, so raising from N = 0 drops degree 0 -- and the raising factor
  // there is sqrt((0 - 0)(0 + 0 + 1)) = 0, so nothing is lost. Lowering goes
  // the other way: the result starts at a lower degree, and the coefficient
  // that has no source is zero.
  std::cout << "f       degrees " << e.MinDegree() << " .. " << e.MaxDegree()
            << "\n"
            << "eth f   degrees " << raised.MinDegree() << " .. "
            << raised.MaxDegree() << "\n\n";

  // Applying eth then eth-bar to a scalar gives the surface Laplacian, whose
  // eigenvalue on Y_{lm} is -l(l+1): the product of the two factors above at
  // N = 0 is -sqrt(l(l+1)) * sqrt(l(l+1)). This identity fixes the relative
  // sign of the two operators and both of their magnitudes. Lower(Raise(e))
  // reads inside out: Raise first, then Lower.
  auto laplacian = Lower(Raise(e));
  auto worst = Real{0};
  for (auto l : e.Degrees()) {
    for (auto m : e.Orders(l)) {
      const auto expected = -static_cast<Real>(l * (l + 1)) * Complex{e[l, m]};
      worst = std::max(worst, std::abs(Complex{laplacian[l, m]} - expected));
    }
  }
  std::cout << "eth-bar eth f against -l(l+1) f: " << worst << "\n";

  // The second identity: the commutator [eth, eth-bar] = eth eth-bar -
  // eth-bar eth is multiplication by -2N. Checked on a spin-1 expansion
  // filled with arbitrary coefficients, so N = 1 and the expected result is
  // -2 g.
  //
  // Neither identity sees the *overall* sign of the pair -- both survive
  // flipping eth and eth-bar together. That sign is fixed jointly with the
  // sign in e_{+-} = -+(theta-hat +- i phi-hat)/sqrt(2), by requiring the
  // surface gradient read back into (theta-hat, phi-hat) to be the gradient;
  // tests/TestConventions.cpp checks it, and the reference states it under
  // "What remains a convention".
  auto g = SpinExpansion<1, Grid>(grid, lMax);
  for (auto l : g.Degrees()) {
    for (auto m : g.Orders(l)) {
      g[l, m] = Complex{std::cos(0.3 * l + m), std::sin(0.2 * l - m)};
    }
  }
  worst = 0;
  for (auto l : g.Degrees()) {
    for (auto m : g.Orders(l)) {
      // Each Raise or Lower builds a whole new expansion, so this inner loop
      // is quadratic in the number of coefficients. That is fine for a check
      // at lMax = 12; in real code compute each composite once.
      const auto commutator =
          Complex{Raise(Lower(g))[l, m]} - Complex{Lower(Raise(g))[l, m]};
      worst = std::max(worst, std::abs(commutator + 2.0 * Complex{g[l, m]}));
    }
  }
  std::cout << "[eth, eth-bar] against -2N:      " << worst << "\n\n";

  // Back to the spatial domain. Evaluate (the inverse transform) turns the
  // raised expansion into an ordinary spin-1 field, on which the whole field
  // algebra of examples 02-04 is available.
  auto gradientPlus = Evaluate(raised);
  static_assert(decltype(gradientPlus)::UpperIndex == 1);
  std::cout << "the raised field is a spin-1 field of " << gradientPlus.Size()
            << " samples\n";
}
