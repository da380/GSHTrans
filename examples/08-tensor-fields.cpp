// 08 -- Tensor fields
//
// A tensor field owns one buffer and hands out its canonical components as
// ordinary spin fields. A component is labelled by a multi-index, and its
// upper index is the signed sum of that multi-index -- which for rank two and
// above are different things.

#include <GSHTrans/All>
#include <complex>
#include <iostream>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  // A rank-2 tensor has components at upper index +-2, so the grid must carry
  // them. Building one on a narrower grid throws at construction.
  auto grid = Grid(16, 2);

  using Tensor = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto t = Tensor(grid);

  std::cout << "components " << Tensor::Components << ", stored "
            << Tensor::StoredComponents << "\n\n";

  // The upper index of each component, which is the count the theory note
  // tabulates: 1, 2, 3, 2, 1 across N = -2 .. 2.
  static_assert(decltype(t.Component<-1, -1>())::UpperIndex == -2);
  static_assert(decltype(t.Component<-1, 0>())::UpperIndex == -1);
  static_assert(decltype(t.Component<-1, 1>())::UpperIndex == 0);
  static_assert(decltype(t.Component<0, 0>())::UpperIndex == 0);
  static_assert(decltype(t.Component<1, -1>())::UpperIndex == 0);
  static_assert(decltype(t.Component<1, 1>())::UpperIndex == 2);

  // Three *distinct* components share N = 0. Their node types are identical,
  // because a node is labelled by its upper index and nothing else -- which is
  // exactly why a collection labelled only by N does not determine a tensor,
  // and why components are addressed by multi-index here.
  t.Component<-1, 1>()[0, 0] = Complex{1.0, 0.0};
  t.Component<0, 0>()[0, 0] = Complex{2.0, 0.0};
  t.Component<1, -1>()[0, 0] = Complex{3.0, 0.0};
  std::cout << "the three N = 0 components: " << (t.Component<-1, 1>()[0, 0])
            << " " << (t.Component<0, 0>()[0, 0]) << " "
            << (t.Component<1, -1>()[0, 0]) << "\n\n";

  // Symmetry is storage rather than bookkeeping: a symmetric tensor stores six
  // of its nine components, and the transposed one *is* the same field.
  using Symmetric2 = TensorField<2, Symmetric<2>, ComplexTensor, Grid>;
  auto s = Symmetric2(grid);
  std::cout << "symmetric: stored " << Symmetric2::StoredComponents << " of "
            << Symmetric2::Components << "\n";
  s.Component<0, 1>()[1, 1] = Complex{7.0, -2.0};
  std::cout << "wrote T^{01}, read T^{10}: " << (s.Component<1, 0>()[1, 1])
            << "\n\n";

  // An antisymmetric tensor goes further: its diagonal vanishes identically,
  // and asking for a component that is identically zero is a compile error
  // rather than a silently zero field. The trait is there so a compile-time
  // traversal can skip them.
  using Skew = TensorField<2, Antisymmetric<2>, ComplexTensor, Grid>;
  static_assert(Skew::Vanishes<0, 0>);
  static_assert(!Skew::Vanishes<0, 1>);
  std::cout << "antisymmetric: stored " << Skew::StoredComponents
            << ", and the diagonal does not exist\n";

  // The elastic symmetry of a rank-4 tensor, c_{ijkl} = c_{jikl} = c_{ijlk} =
  // c_{klij}, comes out at the twenty-one independent components everyone
  // already knows -- without the machinery being told anything about
  // elasticity.
  using Elastic = TensorField<4, ElasticSymmetry, ComplexTensor, Grid>;
  std::cout << "elastic rank 4: stored " << Elastic::StoredComponents
            << " of " << Elastic::Components << "\n";

  FFTWpp::CleanUp();
}
