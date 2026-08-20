// 12 -- The spectral side
//
// An expansion is the spectral counterpart of a field: coefficients indexed by
// degree and order, with the upper index in the type. It is where the raising
// and lowering operators live, because they are local in (l, m) and not in
// position.

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

  // Coefficients are addressed by (l, m). The degrees start at |N|.
  auto e = SpinExpansion<1, Grid>(grid, lMax);
  std::cout << "degrees " << e.MinDegree() << " .. " << e.MaxDegree() << ", "
            << e.Size() << " coefficients\n";

  // A single harmonic: everything zero but one coefficient.
  e[3, 2] = Complex{1.0, 0.0};
  auto field = Evaluate(e);
  static_assert(decltype(field)::UpperIndex == 1);

  // Transforming back recovers it, since a single harmonic is band-limited by
  // construction.
  auto back = Expand(field, lMax);
  std::cout << "Y^1_{3,2} round trip: " << (back[3, 2]) << "\n";

  // Orthonormality, read off the coefficients: every other one is zero.
  auto worst = Real{0};
  for (auto l : back.Degrees()) {
    for (auto m : back.Orders(l)) {
      if (l == 3 && m == 2) continue;
      worst = std::max(worst, std::abs(back[l, m]));
    }
  }
  std::cout << "largest other coefficient: " << worst << "\n\n";

  // A tensor has an expansion too: one buffer, per-component blocks, and the
  // components sharing an upper index are transformed as one batch.
  using Tensor = TensorField<2, Symmetric<2>, RealTensor, Grid>;
  auto t = Tensor(grid);
  t.Component<0, 0>()[2, 2] = 1.0;
  t.Component<-1, -1>()[2, 2] = Complex{0.5, -0.25};

  auto expansion = Expand(t, lMax);
  std::cout << "tensor expansion: " << expansion.Size() << " coefficients for "
            << Tensor::StoredComponents << " stored components\n";

  // A pinned component is a real field, so its block uses the reduced m >= 0
  // storage -- the same saving in the spectral domain as in the spatial one.
  auto pinned = expansion.Component<0, 0>();
  auto general = expansion.Component<-1, -1>();
  std::cout << "  the pinned block  " << pinned.Size() << " coefficients\n"
            << "  a complex block   " << general.Size() << "\n";

  FFTWpp::CleanUp();
}
