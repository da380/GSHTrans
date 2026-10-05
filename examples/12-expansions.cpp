// 12 -- The spectral side
//
// What this shows. An expansion is the spectral counterpart of a field: the
// coefficients f^N_{lm} of
//
//     f = sum_{l >= |N|} sum_{|m| <= l} f^N_{lm} Y^N_{lm}(theta, phi),
//
// indexed by degree l and order m, with the upper index N in the type just as
// it is on a SpinField. The harmonics Y^N_{lm} are orthonormal over the sphere
// (docs/gshtrans-reference.tex, "Generalized spherical harmonics"). The
// expansion is where the raising and lowering operators of example 13 live,
// because they act at fixed (l, m) and not at fixed position.
//
// Read first. 02 (fields), 06 (Expand and Evaluate, and why an arbitrary field
// is projected rather than reproduced), 08 and 10 (tensor fields and their
// reality reduction).
//
// Introduced. SpinExpansion<N, Grid>: MinDegree, MaxDegree, Degrees, Orders,
// Size, and the coefficient accessor e[l, m]. The tensor expansion that
// Expand(tensor, lMax) returns, and its Component<...>() blocks.
//
// The output. The degree range and size of a spin-1 expansion; a single
// harmonic surviving a round trip through the spatial domain with every
// other coefficient at rounding level; and the coefficient counts of a real
// symmetric rank-2 tensor's expansion, showing the reduced storage of its
// pinned components.

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

  // A grid resolving degree 12 and upper indices up to |N| = 2, which is
  // enough for the spin-1 expansion and for the rank-2 tensor below.
  constexpr auto lMax = Int{12};
  auto grid = Grid(lMax, 2);

  // A zero expansion of upper index 1, up to degree lMax. Coefficients are
  // addressed by (l, m) with -l <= m <= l, and the degrees start at |N| = 1:
  // d^l_{Nm} vanishes identically below l = |N|, so there is no harmonic there
  // and nothing to store. Size() counts the (l, m) pairs that are stored.
  auto e = SpinExpansion<1, Grid>(grid, lMax);
  std::cout << "degrees " << e.MinDegree() << " .. " << e.MaxDegree() << ", "
            << e.Size() << " coefficients\n";

  // A single harmonic, Y^1_{3,2}: one coefficient set to one, the rest zero.
  // Evaluate is the inverse transform; it returns a SpinField whose upper
  // index is the expansion's, checked here at compile time.
  e[3, 2] = Complex{1.0, 0.0};
  auto field = Evaluate(e);
  static_assert(decltype(field)::UpperIndex == 1);

  // Expand is the forward transform. A field built from coefficients up to
  // lMax is band-limited by construction, so expanding it again recovers them
  // to rounding -- the round trip that example 06 shows can fail for a field
  // given as an arbitrary function.
  auto back = Expand(field, lMax);
  std::cout << "Y^1_{3,2} round trip: " << (back[3, 2]) << "\n";

  // Orthonormality, read off the coefficients: Y^1_{3,2} projects onto no
  // other harmonic, so every other coefficient is zero to rounding. Degrees()
  // and Orders(l) iterate exactly the stored (l, m).
  auto worst = Real{0};
  for (auto l : back.Degrees()) {
    for (auto m : back.Orders(l)) {
      if (l == 3 && m == 2) continue;
      worst = std::max(worst, std::abs(back[l, m]));
    }
  }
  std::cout << "largest other coefficient: " << worst << "\n\n";

  // A tensor field has an expansion too, mirroring the TensorField layout of
  // examples 08-10: one buffer, holding one coefficient block per *stored*
  // component, each block sized by that component's own upper index. The
  // stored components sharing an upper index are transformed as one batch
  // (example 07), since a batch must share its upper index.
  using Tensor = TensorField<2, Symmetric<2>, RealTensor, Grid>;
  auto t = Tensor(grid);
  // Two components set at one grid point -- a spatial field is indexed
  // [iTheta, iPhi], not [l, m].
  t.Component<0, 0>()[2, 2] = 1.0;
  t.Component<-1, -1>()[2, 2] = Complex{0.5, -0.25};

  auto expansion = Expand(t, lMax);
  std::cout << "tensor expansion: " << expansion.Size() << " coefficients for "
            << Tensor::StoredComponents << " stored components\n";

  // Component<...>() on an expansion is a spin expansion over that block.
  // (0, 0) is a pinned component of a real tensor (example 10): a real field
  // at N = 0, so its block uses the reduced m >= 0 storage, the negative
  // orders following from f_{l,-m} = (-1)^m conj(f_{lm}). (-1, -1) sits at
  // N = -2 and is a general complex field, stored for every m and only from
  // l = 2 upward. The two sizes printed show the difference.
  auto pinned = expansion.Component<0, 0>();
  auto general = expansion.Component<-1, -1>();
  std::cout << "  the pinned block  " << pinned.Size() << " coefficients\n"
            << "  a complex block   " << general.Size() << "\n";
}
