// 11 -- Elasticity
//
// The case rank 4 exists for: a stress field from an elastic tensor and a
// strain, which is a double contraction of a tensor product. Every
// intermediate's rank and upper index is checked by the compiler.
//
// The mathematics, briefly. Hooke's law, sigma^{ij} = c^{ijkl} e_{kl}. In
// canonical components every index is contravariant, so lowering the
// strain's indices is a contraction with the metric g_ab = (-1)^a
// delta_{a+b,0}: sigma^{ij} = sum c^{ijkl} g_{km} g_{ln} e^{mn}. That is the
// tensor product c (x) e, of rank 6, contracted twice. See
// docs/gshtrans-reference.tex, section "Canonical components".
//
// What this shows
//   Real rank-4 and rank-2 fields with their symmetries, the tensor product,
//   contraction of chosen slots, and materialising the result with a stated
//   symmetry and reality.
//
// Assumes
//   Examples 08-10.
//
// Introduced
//   TensorProduct, Contract<J, K>, Materialise<Symmetry, Reality>.
//
// Output
//   Storage counts for the elastic tensor and the strain, one component of
//   the stress at a point, the integral of its trace, and the storage of the
//   materialised stress.

#include <GSHTrans/GSHTrans.hpp>
#include <cmath>
#include <complex>
#include <iostream>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  // Rank 4 needs upper indices to 4. The rank-6 product below reaches
  // upper index 6, but it is never stored, and an expression may carry upper
  // indices the grid does not.
  auto grid = Grid(16, 4);

  // The elastic symmetry c_{ijkl} = c_{jikl} = c_{ijlk} = c_{klij} is given
  // by three generators; the orbit machinery works out that it leaves 21
  // independent components, which is the familiar count. As a RealTensor
  // those 21 real degrees of freedom are stored as 8 complex components and
  // 5 pinned real ones, 13 in all. The strain is a real symmetric rank-2
  // tensor: 6 reals, as 2 complex and 2 real components.
  using Elastic = TensorField<4, ElasticSymmetry, RealTensor, Grid>;
  using Strain = TensorField<2, Symmetric<2>, RealTensor, Grid>;

  std::cout << "elastic tensor: " << Elastic::StoredComponents
            << " stored components, " << Elastic::RealsPerPoint
            << " reals a point\n"
            << "strain:         " << Strain::StoredComponents << " stored, "
            << Strain::RealsPerPoint << " reals a point\n\n";

  auto c = Elastic(grid);
  auto e = Strain(grid);

  // Fill the stored components. A pinned component is a real field, so the
  // write is of a real number; the others are complex. Of those written
  // below, c^{0000}, c^{-+00}, e^{00} and e^{-+} are pinned and e^{--} is
  // complex. Every other component stays zero or is derived from these.
  const auto fill = [&](auto&& u, Real tag) {
    using Node = std::remove_cvref_t<decltype(u)>;
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        if constexpr (std::same_as<typename Node::Value, RealValued>) {
          u[iTheta, iPhi] = tag * (1 + 0.1 * std::cos(0.3 * iTheta));
        } else {
          u[iTheta, iPhi] = Complex{tag * std::cos(0.3 * iTheta),
                                    0.1 * tag * std::sin(0.2 * iPhi)};
        }
      }
    }
  };
  fill(c.Component<0, 0, 0, 0>(), 2.0);
  fill(c.Component<-1, 1, 0, 0>(), 0.5);
  fill(e.Component<0, 0>(), 1.0);
  fill(e.Component<-1, 1>(), 0.25);
  fill(e.Component<-1, -1>(), 0.1);

  // Read-only from here on, through const references (example 09); the
  // expressions below refer to these named tensors, which must outlive them.
  const auto& elastic = c;
  const auto& strain = e;

  // c^{ijkl} e^{mn}, then contract (k, m) and what is left of (l, n). The
  // metric comes with each contraction, so this is the correct
  // index-raised-and-lowered product and not a naive sum over components.
  // TensorProduct concatenates the slots, (c (x) e)^{ijklmn} = c^{ijkl}
  // e^{mn}, so its upper indices add.
  auto product = TensorProduct(elastic, strain);
  static_assert(decltype(product)::Rank == 6);

  // Contract<J, K> contracts slots J and K, counted from zero, against the
  // metric. Slots 2 and 4 of the product are k and m; what remains is
  // i j l n, in that order.
  auto once = Contract<2, 4>(product);
  static_assert(decltype(once)::Rank == 4);

  // Slots 2 and 3 of that are l and n, leaving sigma^{ij}.
  auto stress = Contract<2, 3>(once);
  static_assert(decltype(stress)::Rank == 2);
  static_assert(decltype(stress.Component<1, 1>())::UpperIndex == 2);

  std::cout << "stress^{00} at a point " << (stress.Component<0, 0>()[3, 3])
            << "\n";

  // And the trace of the stress, contracted against the metric and integrated
  // over the sphere -- a scalar built from a rank-6 intermediate without any
  // of it being stored.
  std::cout << "integral of tr(stress) " << Integrate(Trace(stress)) << "\n\n";

  // Nothing has been stored yet: the point printed and the integral each
  // walked the lazy chain on demand. Materialise evaluates it into a field.
  // The stress is symmetric because c^{ijkl} = c^{jikl}, and real because c
  // and e are, so it may be stored as a real symmetric tensor: six reals a
  // point rather than the eighteen of a general complex rank-2 tensor. As in
  // example 09, the stated symmetry and reality are the caller's assertion
  // and are not checked.
  auto stored = Materialise<Symmetric<2>, RealTensor>(stress);
  std::cout << "materialised stress: " << decltype(stored)::StoredComponents
            << " stored components, " << decltype(stored)::RealsPerPoint
            << " reals a point\n";
}
