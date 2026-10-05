// 17 -- Tangential tensors, and the derivative that is closed on them
//
// What this shows. Most geophysical surface objects have no radial slot:
// surface strain and stress, the metric of the sphere, the spin-2 fields of
// geodesy and the CMB. Their indices run over {-1, +1} rather than
// {-1, 0, +1}, so a rank-p one has 2^p components rather than 3^p. The slot
// alphabet is a template parameter of the tensor types: AllSlots (the
// default) or TangentialSlots.
//
// A tangential tensor is not a general tensor with some components left at
// zero. It is another object in another bundle, and the two are joined by an
// explicit embedding (Embed) and projection (Tangential) rather than by a
// coincidence of storage. What makes the type worth having is not only the
// saving but the derivative. There are two:
//
//   - SurfaceGradient (example 15) is the *ambient* derivative. It
//     differentiates the basis too, and the basis leaves the tangent plane,
//     so its result has radial slots even when its operand has none.
//   - IntrinsicDerivative is the Levi-Civita connection of the sphere's own
//     metric. On a tangential tensor every connection term of the ambient
//     formula would land on a radial slot, which does not exist, so what is
//     left is a pure Omega multiplication -- eth and eth-bar component by
//     component, up to sqrt(2) -- and the result is tangential again.
//
// The two are related by the Gauss formula: the tangential block of the
// ambient gradient is the intrinsic one, and the block with one radial slot
// is the extrinsic curvature, which on the unit sphere is exactly minus the
// operand. docs/gshtrans-reference.tex, "The intrinsic derivative, where eth
// is exactly right", gives the derivation.
//
// Read first. 08-10 (tensor fields and the reality reduction), 13 (eth),
// 15 (the surface gradient).
//
// Introduced. TangentialRank2Field, TangentialSymmetricField and the
// TangentialSlots parameter; Represents<...>; Embed and Tangential, on fields
// and on expansions; IntrinsicDerivative.
//
// The output. Component and storage counts in the two bundles; an
// Embed/Tangential round trip; and the split of the ambient gradient into the
// intrinsic derivative and the curvature term, both to rounding.

#include <GSHTrans/GSHTrans.hpp>
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
  // not compile. `Represents<...>` is the compile-time question "does this
  // tensor have that component?", which generic code asks before touching
  // one.
  static_assert(Tangential2::Represents<-1, 1>);
  static_assert(!Tangential2::Represents<0, 1>);

  // The reality reduction of example 10 applies here too, and more simply.
  // For a real tensor T^{-a...} = (-1)^N conj(T^{a...}), which pairs each
  // multi-index with its negation. In the general bundle the all-zero index
  // is its own negation, so it is pinned to a single real field; without a
  // zero letter, and without symmetry, there is no such fixed point, so the
  // components pair off exactly and nothing is pinned: four complex
  // components give 4 reals a point. With symmetry, (-1, +1) and (+1, -1)
  // are one component that is its own negation, so it is pinned real, and a
  // real symmetric tangential tensor comes to 3 reals a point.
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

  // A complex tangential rank-2 field, every component filled with
  // arbitrary values on the grid.
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

  // Embed widens the alphabet. On fields it is a lazy node (example 04):
  // nothing is allocated, and the components with a radial slot are ones the
  // embedded tensor does not represent rather than zeros it stores. Every
  // traversal -- materialising, a product, a contraction -- reads an
  // unrepresented component as zero and skips it.
  const auto embedded = Embed(tangential);
  static_assert(std::same_as<decltype(embedded)::SlotSet, AllSlots>);
  static_assert(!decltype(embedded)::Represents<0, 1>);

  // Tangential projects back by dropping every component with a radial slot,
  // so Tangential(Embed(t)) is t. Also lazy; the printed difference is read
  // at one grid point.
  const auto back = Tangential(Embed(tangential));
  std::cout << "Tangential(Embed(t)) - t at one point: "
            << std::abs(back.Component<-1, 1>()[2, 2] -
                        tangential.Component<-1, 1>()[2, 2])
            << "\n\n";

  // The embedding is also how a product crosses bundles, which otherwise does
  // not compile: a tangential factor and a general one are not tensors over a
  // common alphabet. TensorProduct of a rank-2 and a rank-1 is rank 3.
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
  // On expansions Embed copies into a general expansion rather than building
  // a lazy node; the spectral side has no lazy nodes.
  //
  // The intrinsic derivative is the Levi-Civita connection of the induced
  // metric, and it is closed. By the Gauss formula it is exactly the
  // tangential block of the ambient one, so it carries D&T's Omega factors
  // and not eth's, which differ by sqrt(2).
  //
  // The operand is a tangential vector expansion -- the fifth template
  // argument selects the alphabet -- filled with arbitrary coefficients
  // straight through its buffer.
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
      // The tangential block of the ambient gradient is the intrinsic one;
      // (+1, -1) is checked here as a representative.
      worst = std::max(worst, std::abs(intrinsic.Coefficient<1, -1>(l, m) -
                                       ambient.Coefficient<1, -1>(l, m)));
      // And the block with a radial slot is the extrinsic curvature: minus
      // the field with that slot replaced by the derivative index, here
      // (grad T)^{+1, 0} = -T^{+1}.
      curvature = std::max(curvature, std::abs(ambient.Coefficient<1, 0>(l, m) +
                                               line.Coefficient<1>(l, m)));
    }
  }
  std::cout << "the split of the ambient gradient\n"
            << "  tangential block against IntrinsicDerivative " << worst
            << '\n'
            << "  radial block against -T                      " << curvature
            << '\n';

  // Closed means it applies again without leaving the bundle: the second
  // derivative is a rank-2 tangential tensor. The ambient operator cannot do
  // that; its result already lives in the general bundle.
  const auto second = IntrinsicDerivative(IntrinsicDerivative(line));
  std::cout << "\nIntrinsicDerivative applied twice: rank "
            << decltype(second)::Rank << ", still tangential\n";

  return 0;
}
