// 15 -- The surface gradient
//
// What this shows. The derivative a tensor field wants: the contravariant
// derivative of Phinney & Burridge and Dahlen & Tromp (C.151)-(C.153), which
// takes the coefficients of a rank-q tensor T to those of the rank-(q+1)
// tensor grad_1 T on the unit sphere. With
// Omega^{+-N}_l = sqrt((l +- N)(l -+ N + 1)/2), its components are
//
//     (grad_1 T)^{sigma a_1...a_q}_{lm} = Omega^{-sigma N}_l T^{a_1...a_q}_{lm}
//                                  - sum_i T^{a_1...(a_i + sigma)...a_q}_{lm}
//
// for sigma = +-1, where N is the upper index of T^{a_1...a_q} (so sigma = -1
// takes Omega^{+N} and sigma = +1 takes Omega^{-N}), and a shifted index
// outside {-1, 0, +1} contributes nothing.
// The sigma = 0 component is zero: the surface gradient has no radial part.
//
// It is *not* eth applied to each component. On a scalar the two agree up to
// sqrt(2) (the normalisation of e_{+-}), but from rank one upward the sum
// above -- the connection terms -- appears, because the canonical basis
// vectors e_{-1}, e_0, e_{+1} themselves vary over the sphere and
// differentiating a tensor differentiates its basis too.
// docs/gshtrans-reference.tex, "The contravariant derivative, and how it
// differs", has the derivation.
//
// Read first. 08-10 (tensor fields, canonical components, the metric trace,
// the reality reduction), 12 (expansions) and 13 (eth and eth-bar).
//
// Introduced. TensorExpansion built directly; SurfaceGradient(expansion);
// Coefficient<...>(l, m), which reads any component of a tensor expansion,
// stored or derived.
//
// The output. The radial component of a gradient (zero); grad^+ of a scalar
// against -eth/sqrt(2); the metric trace of the second gradient against the
// surface Laplacian; the storage of the gradient of a real vector; and a
// derived spectral component checked against the reality relation.

#include <GSHTrans/GSHTrans.hpp>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  using Int = std::ptrdiff_t;

  constexpr auto lMax = Int{10};
  auto grid = Grid(lMax, 2);

  // A scalar field, as a rank-0 complex tensor expansion with arbitrary
  // coefficients. It has one component, addressed with an empty multi-index:
  // Component<>(), and later Coefficient<>(l, m).
  auto f = TensorExpansion<0, NoSymmetry<0>, ComplexTensor, Grid>(grid, lMax);
  for (auto l = Int{0}; l <= lMax; l++) {
    for (auto m = -l; m <= l; m++) {
      f.Component<>()[l, m] = Complex{std::cos(0.4 * l + m), 0.2 * m};
    }
  }

  // The gradient prepends a slot: rank 0 becomes rank 1, and the new leading
  // index sigma takes values -1, 0, +1. A component's upper index is the sum
  // of its slot indices (example 08), so (grad f)^{+1} sits at N = +1 and
  // (grad f)^{-1} at N = -1. The result has no symmetry and keeps the
  // operand's reality.
  auto grad = SurfaceGradient(f);
  static_assert(decltype(grad)::Rank == 1);
  static_assert(decltype(grad.Component<-1>())::UpperIndex == -1);
  static_assert(decltype(grad.Component<1>())::UpperIndex == 1);

  std::cout << std::scientific << std::setprecision(3);

  // There is no radial component. The surface gradient has no e_0 part at
  // all (D&T C.145), so that block is stored and zero -- which keeps the
  // result an ordinary rank-1 tensor, composable with contraction, further
  // gradients and the transform. Example 16 fills that slot with d/dr to
  // make the full three-dimensional gradient.
  std::cout << "radial component of grad f: "
            << std::abs(grad.Coefficient<0>(4, 2)) << "\n";

  // On a scalar the operator is eth, up to the sqrt(2) that is the
  // normalisation of e_+- = -+(theta-hat +- i phi-hat)/sqrt(2): from the
  // factors in example 13, grad^+ f = -(eth f)/sqrt(2) and
  // grad^- f = +(eth-bar f)/sqrt(2). Raise acts on a SpinExpansion, so the
  // coefficients are copied into one first.
  auto asSpin = SpinExpansion<0, Grid>(grid, lMax);
  for (auto l : asSpin.Degrees()) {
    for (auto m : asSpin.Orders(l)) asSpin[l, m] = f.Coefficient<>(l, m);
  }
  const auto rootTwo = std::numbers::sqrt2_v<Real>;
  const auto viaEth = -Complex{Raise(asSpin)[4, 2]} / rootTwo;
  std::cout << "grad^+ against -eth/sqrt(2): "
            << std::abs(grad.Coefficient<1>(4, 2) - viaEth) << "\n\n";

  // Applying it twice and taking the metric trace gives the surface
  // Laplacian, eigenvalue -l(l+1). The trace contracts the two slots with
  // g_{ab} = (-1)^a delta_{a+b,0}, so it is
  //
  //     -T^{-1,+1} + T^{0,0} - T^{+1,-1}.
  //
  // This runs through the connection terms: d^- of the +1 component
  // subtracts the 0 component, which eth knows nothing about, so the identity
  // would fail if this were eth component by component.
  auto second = SurfaceGradient(grad);
  static_assert(decltype(second)::Rank == 2);

  auto worst = Real{0};
  for (auto l = Int{0}; l <= lMax; l++) {
    for (auto m = -l; m <= l; m++) {
      const auto trace = -second.Coefficient<-1, 1>(l, m) +
                         second.Coefficient<0, 0>(l, m) -
                         second.Coefficient<1, -1>(l, m);
      const auto expected =
          -static_cast<Real>(l * (l + 1)) * f.Coefficient<>(l, m);
      worst = std::max(worst, std::abs(trace - expected));
    }
  }
  std::cout << "trace of the second gradient against -l(l+1) f: " << worst
            << "\n\n";

  // The gradient of a real vector field is a real rank-2 tensor, and only the
  // reduced set is stored: the reality condition supplies the rest
  // (example 10). Here the vector is filled on the grid -- index [iTheta,
  // iPhi] -- through its stored components: v^{-1} is complex, and v^0 is
  // pinned real; v^{+1} = -conj(v^{-1}) is derived and not written.
  using RealVector = TensorField<1, NoSymmetry<1>, RealTensor, Grid>;
  auto v = RealVector(grid);
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      v.Component<-1>()[iTheta, iPhi] =
          Complex{std::cos(0.3 * iTheta), std::sin(0.2 * iPhi)};
      v.Component<0>()[iTheta, iPhi] = 1.0 + 0.5 * std::cos(0.4 * iTheta);
    }
  }

  const auto& vector = v;
  auto gradV = SurfaceGradient(Expand(vector, lMax));
  static_assert(decltype(gradV)::Rank == 2);
  static_assert(std::same_as<decltype(gradV)::Reality, RealTensor>);

  std::cout << std::defaultfloat;
  std::cout << "grad of a real vector: rank " << decltype(gradV)::Rank << ", "
            << decltype(gradV)::StoredComponents << " stored components of "
            << decltype(gradV)::Components << ", " << gradV.Size()
            << " coefficients\n";

  // Every component is readable through Coefficient whether or not it is
  // stored, which is what the derived half of a real tensor needs in the
  // spectral domain. Deriving one there is not the pointwise relation it is
  // on the sphere: for a real tensor it is T^{-N}_{l,-m} = (-1)^m
  // conj(T^N_{lm}), which *reverses the order*, so the derived coefficient at
  // +m comes from the stored one at -m (docs/gshtrans-reference.tex,
  // "Reality"). At m = 2 the sign is +1, so the two printed values are
  // conjugates. Component<...>() gives only stored blocks, for this reason.
  const auto stored = gradV.Coefficient<-1, -1>(4, -2);
  const auto derived = gradV.Coefficient<1, 1>(4, 2);
  std::cout << "  stored  (grad v)^{-1,-1} at (l,m) = (4,-2) " << stored << "\n"
            << "  derived (grad v)^{+1,+1} at (l,m) = (4, 2) " << derived
            << "\n"
            << "  the relation says the second is conj of the first: "
            << (std::abs(derived - std::conj(stored)) < 1.0e-12 ? "yes" : "no")
            << "\n";
}
