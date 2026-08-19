// The spin-field algebra: lazy, index-checked expressions over a grid.

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <iostream>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  constexpr auto band = 8;
  constexpr auto nMax = 2;

  // Enough quadrature headroom for a product of two band-limited fields.
  auto grid = Grid::ForBand(band, nMax, 2.0);

  auto u = SpinField<2, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::sin(theta) * std::cos(phi), std::sin(2 * phi)};
  });
  auto v = SpinField<2, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::cos(theta), 0.5 * std::sin(phi)};
  });
  auto w = SpinField<0, Grid, RealValued>(
      grid, [](auto theta, auto) { return 1.0 + 0.25 * std::cos(theta); });

  // conj(u) carries upper index -2, so the product lands at zero and can be
  // integrated. Nothing is evaluated until Integrate walks the expression.
  const auto pairing = Integrate(conj(u) * v);
  const auto energy = Integrate(abs2(u));

  std::cout << "<u, v> = " << pairing << "\n";
  std::cout << "||u||  = " << std::sqrt(energy) << "\n";

  // Scalar fields multiply anything without changing its upper index.
  auto scaled = Materialise(u * w);
  static_assert(decltype(scaled)::UpperIndex == 2);

  // Assignment evaluates in place; the destination may appear on the right.
  // The right-hand side must still land at upper index 2: u * conj(u) is at
  // zero, so multiplying it in leaves the index alone. Writing
  // conj(scaled) * (u * conj(u)) instead would land at -2 and not compile.
  scaled += scaled * (u * conj(u));

  std::cout << "||u w||  = " << std::sqrt(Integrate(abs2(scaled))) << "\n";

  FFTWpp::CleanUp();
}
