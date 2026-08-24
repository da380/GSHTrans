// 17 -- Tangential tensors, and the derivative that is closed on them
//
// Most geophysical surface objects have no radial slot: surface strain and
// stress, the metric of the sphere, the spin-2 fields of geodesy and the CMB.
// Their indices run over {-1, +1} rather than {-1, 0, +1}, so a rank-p one has
// 2^p components rather than 3^p.
//
// That is not a general tensor with some components left at zero. It is
// another object in another bundle, and the two are joined by an embedding
// rather than by a coincidence of storage. What makes the type worth having is
// not the saving but the derivative: on tangential tensors the intrinsic
// covariant derivative is closed, and the ambient one is not.

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  using Int = std::ptrdiff_t;

  constexpr auto lMax = Int{8};
  auto grid = Grid(lMax, 3);

  std::cout << std::scientific << std::setprecision(3);

  //------------------------------------------------------------------------//
  // The bundle is smaller, and the type says which one it is
  //------------------------------------------------------------------------//

  using General2 = Rank2TensorField<Grid, ComplexTensor>;
  using Tangential2 = TangentialRank2Field<Grid, ComplexTensor>;

  std::cout << "components of a rank-2 tensor\n"
            << "  general    " << General2::Components << '\n'
            << "  tangential " << Tangential2::Components << "\n\n";

  // A radial index is not a component of a tangential tensor -- not a
  // component that is zero, but one that does not exist. Asking for it does
  // not compile, and `Represents` is what a traversal asks first.
  static_assert(Tangential2::Represents<-1, 1>);
  static_assert(!Tangential2::Represents<0, 1>);

  // Reality reduces further, and here it reduces further than it does in the
  // general bundle: negation has no fixed point without a zero letter, so
  // nothing is pinned to a single real number and the reality reduction's
  // second buffer is
  // empty.
  using RealTangential2 = TangentialRank2Field<Grid>;
  std::cout << "a real rank-2 tensor, reals per point\n"
            << "  general    " << Rank2TensorField<Grid>::RealsPerPoint << '\n'
            << "  tangential " << RealTangential2::RealsPerPoint << '\n'
            << "  tangential, symmetric "
            << TangentialSymmetricField<Grid>::RealsPerPoint
            << "   (a real symmetric 2x2 matrix)\n\n";

  //------------------------------------------------------------------------//
  // Crossing between the bundles is explicit
  //------------------------------------------------------------------------//

  auto t = Tangential2(grid);
  const auto fill = [&](auto&& u, Real tag) {
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        u[iTheta, iPhi] =
            Complex{tag + std::cos(0.3 * iTheta), std::sin(0.2 * iPhi)};
      }
    }
  };
  fill(t.Component<-1, -1>(), 1.0);
  fill(t.Component<-1, 1>(), 2.0);
  fill(t.Component<1, -1>(), 3.0);
  fill(t.Component<1, 1>(), 4.0);

  const auto& tangential = t;

  // Embed widens the alphabet. It is a lazy node: nothing is allocated, and
  // the components with a radial slot are ones the embedded tensor does not
  // have rather than zeros it stores.
  const auto embedded = Embed(tangential);
  static_assert(std::same_as<decltype(embedded)::SlotSet, AllSlots>);
  static_assert(!decltype(embedded)::Represents<0, 1>);

  // Tangential projects back, and the round trip is the identity.
  const auto back = Tangential(Embed(tangential));
  std::cout << "Tangential(Embed(t)) - t at one point: "
            << std::abs(back.Component<-1, 1>()[2, 2] -
                        tangential.Component<-1, 1>()[2, 2])
            << "\n\n";

  // The embedding is also how a product crosses bundles, which otherwise does
  // not compile: a tangential factor and a general one are not tensors over a
  // common alphabet.
  auto v = VectorField<Grid, ComplexTensor>(grid);
  fill(v.Component<0>(), 5.0);
  const auto& general = v;
  const auto product = TensorProduct(Embed(tangential), general);
  static_assert(decltype(product)::Rank == 3);

  //------------------------------------------------------------------------//
  // Two derivatives, and only one of them is closed
  //------------------------------------------------------------------------//

  // The ambient surface gradient differentiates the tensor *and* its basis,
  // and the basis leaves the tangent plane. So it takes a general operand and
  // returns a general result: a tangential tensor is embedded first, visibly.
  //
  // The intrinsic derivative is the Levi-Civita connection of the induced
  // metric, and it is closed. Gauss says it is the tangential block of the
  // ambient one -- exactly, with D&T's Omega, which is why it is not `eth`'s
  // normalisation up to a root two but the block itself.
  auto u =
      TensorExpansion<1, NoSymmetry<1>, ComplexTensor, Grid, TangentialSlots>(
          grid, lMax);
  for (Int i = 0; i < u.Size(); i++) {
    u.Data()[i] = Complex{std::cos(0.17 * i), std::sin(0.37 * i)};
  }
  const auto& line = u;

  const auto intrinsic = IntrinsicDerivative(line);
  const auto ambient = SurfaceGradient(Embed(line));

  static_assert(std::same_as<decltype(intrinsic)::SlotSet, TangentialSlots>);
  static_assert(std::same_as<decltype(ambient)::SlotSet, AllSlots>);

  auto worst = Real{0};
  auto curvature = Real{0};
  for (auto l = Int{0}; l <= lMax; l++) {
    for (auto m = -l; m <= l; m++) {
      // The tangential block of the ambient gradient is the intrinsic one.
      worst = std::max(worst, std::abs(intrinsic.Coefficient<1, -1>(l, m) -
                                       ambient.Coefficient<1, -1>(l, m)));
      // And the block with a radial slot is the extrinsic curvature: minus
      // the field with that slot replaced.
      curvature = std::max(curvature, std::abs(ambient.Coefficient<1, 0>(l, m) +
                                               line.Coefficient<1>(l, m)));
    }
  }
  std::cout << "the split of the ambient gradient\n"
            << "  tangential block against IntrinsicDerivative " << worst
            << '\n'
            << "  radial block against -T                      " << curvature
            << '\n';

  // Closed means it applies again without leaving the bundle, which the
  // ambient one cannot do without another embedding.
  const auto second = IntrinsicDerivative(IntrinsicDerivative(line));
  std::cout << "\nIntrinsicDerivative applied twice: rank "
            << decltype(second)::Rank << ", still tangential\n";

  return 0;
}
