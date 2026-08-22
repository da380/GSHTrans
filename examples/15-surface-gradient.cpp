// 15 -- The surface gradient
//
// The derivative a tensor field wants, in the formalism of Phinney & Burridge
// and Dahlen & Tromp. It is *not* eth applied to each component: for a scalar
// the two agree up to sqrt(2), but from rank one upward the operator also
// carries connection terms, because the canonical basis vectors themselves
// vary over the sphere.

#include <GSHTrans/All>
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

  // A scalar field, as a rank-0 tensor. Its expansion is one block.
  auto f = TensorExpansion<0, NoSymmetry<0>, ComplexTensor, Grid>(grid, lMax);
  for (auto l = Int{0}; l <= lMax; l++) {
    for (auto m = -l; m <= l; m++) {
      f.Component<>()[l, m] = Complex{std::cos(0.4 * l + m), 0.2 * m};
    }
  }

  // The gradient prepends a slot: rank 0 becomes rank 1, and the new index
  // takes values -1, 0, +1 with the component's upper index the signed sum.
  auto grad = SurfaceGradient(f);
  static_assert(decltype(grad)::Rank == 1);
  static_assert(decltype(grad.Component<-1>())::UpperIndex == -1);
  static_assert(decltype(grad.Component<1>())::UpperIndex == 1);

  std::cout << std::scientific << std::setprecision(3);

  // There is no radial component. The surface gradient has no e_0 part at
  // all, so that block is stored and zero -- which keeps the result an
  // ordinary tensor, composable with everything else.
  std::cout << "radial component of grad f: "
            << std::abs(grad.Coefficient<0>(4, 2)) << "\n";

  // On a scalar the operator is eth, up to the sqrt(2) that is the
  // normalisation of e_+- = -+(theta-hat +- i phi-hat)/sqrt(2).
  auto asSpin = SpinExpansion<0, Grid>(grid, lMax);
  for (auto l : asSpin.Degrees()) {
    for (auto m : asSpin.Orders(l)) asSpin[l, m] = f.Coefficient<>(l, m);
  }
  const auto rootTwo = std::numbers::sqrt2_v<Real>;
  const auto viaEth = -Complex{Raise(asSpin)[4, 2]} / rootTwo;
  std::cout << "grad^+ against -eth/sqrt(2): "
            << std::abs(grad.Coefficient<1>(4, 2) - viaEth) << "\n\n";

  // Applying it twice and taking the metric trace gives the surface
  // Laplacian, eigenvalue -l(l+1). This runs through the connection terms:
  // d^- of the +1 component subtracts the 0 component, which eth knows
  // nothing about, so the identity would fail if this were eth component by
  // component.
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
  // reduced set is stored: the reality condition supplies the rest.
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
  std::cout << "grad of a real vector: rank " << decltype(gradV)::Rank
            << ", " << decltype(gradV)::StoredComponents
            << " stored components of "
            << decltype(gradV)::Components << ", "
            << gradV.Size() << " coefficients\n";

  // Every component is readable whether or not it is stored, which is what
  // the derived half of a real tensor needs in the spectral domain -- and
  // deriving one there is not the pointwise relation it is on the sphere. It
  // is T^{-N}_{l,-m} = (-1)^m conj(T^N_{lm}), which *reverses the order*, so
  // the derived coefficient at +m comes from the stored one at -m.
  const auto stored = gradV.Coefficient<-1, -1>(4, -2);
  const auto derived = gradV.Coefficient<1, 1>(4, 2);
  std::cout << "  stored  (grad v)^{-1,-1} at (l,m) = (4,-2) " << stored << "\n"
            << "  derived (grad v)^{+1,+1} at (l,m) = (4, 2) " << derived
            << "\n"
            << "  the relation says the second is conj of the first: "
            << (std::abs(derived - std::conj(stored)) < 1.0e-12 ? "yes" : "no")
            << "\n";
}
