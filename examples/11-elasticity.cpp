// 11 -- Elasticity
//
// The case rank 4 exists for: a stress field from an elastic tensor and a
// strain, which is a double contraction of a tensor product. Every
// intermediate's rank and upper index is checked by the compiler.

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <iostream>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  // Rank 4 needs upper indices to 4.
  auto grid = Grid(16, 4);

  // The elastic symmetry c_{ijkl} = c_{jikl} = c_{ijlk} = c_{klij} is given
  // by three generators; the orbit machinery works out that it leaves 21
  // independent components, which is the familiar count.
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
  // write is of a real number; the others are complex.
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

  const auto& elastic = c;
  const auto& strain = e;

  // c^{ijkl} e^{mn}, then contract (k, m) and what is left of (l, n). The
  // metric comes with each contraction, so this is the correct
  // index-raised-and-lowered product and not a naive sum over components.
  auto product = TensorProduct(elastic, strain);
  static_assert(decltype(product)::Rank == 6);

  auto once = Contract<2, 4>(product);
  static_assert(decltype(once)::Rank == 4);

  auto stress = Contract<2, 3>(once);
  static_assert(decltype(stress)::Rank == 2);
  static_assert(decltype(stress.Component<1, 1>())::UpperIndex == 2);

  std::cout << "stress^{00} at a point " << (stress.Component<0, 0>()[3, 3])
            << "\n";

  // And the trace of the stress, contracted against the metric and integrated
  // over the sphere -- a scalar built from a rank-6 intermediate without any
  // of it being stored.
  std::cout << "integral of tr(stress) " << Integrate(Trace(stress)) << "\n\n";

  // Nothing has been computed yet beyond the two points printed: the whole
  // chain is lazy. Materialise evaluates it, and asking for a real symmetric
  // result stores six reals a point rather than eighteen.
  auto stored = Materialise<Symmetric<2>, RealTensor>(stress);
  std::cout << "materialised stress: " << decltype(stored)::StoredComponents
            << " stored components, " << decltype(stored)::RealsPerPoint
            << " reals a point\n";

  FFTWpp::CleanUp();
}
