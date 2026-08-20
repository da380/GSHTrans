// 09 -- Tensor algebra
//
// An operation on tensors is an operation on the component accessor, so the
// algebra is lazy and needs no second expression system: what comes back from
// a component is a spin field, and its index arithmetic is already checked.

#include <GSHTrans/All>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  using Int = std::ptrdiff_t;

  auto grid = Grid(16, 2);
  using Tensor = TensorField<2, NoSymmetry<2>, ComplexTensor, Grid>;
  auto t = Tensor(grid);

  const auto fill = [&](auto&& u, Real tag) {
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        u[iTheta, iPhi] = Complex{tag + std::cos(0.3 * iTheta),
                                  std::sin(0.2 * iPhi) - tag};
      }
    }
  };
  fill(t.Component<-1, 1>(), 1.0);
  fill(t.Component<0, 0>(), 2.0);
  fill(t.Component<1, -1>(), 3.0);
  fill(t.Component<0, 1>(), 4.0);
  fill(t.Component<1, 0>(), 5.0);

  const auto& tensor = t;

  // Transposition relabels slots and copies nothing.
  auto transposed = Transpose(tensor);
  std::cout << "T^{01}          " << (tensor.Component<0, 1>()[2, 2]) << "\n"
            << "(T^T)^{01}      " << (transposed.Component<0, 1>()[2, 2])
            << "\n\n";

  // The trace contracts against the metric g_{ab} = (-1)^a delta_{a+b,0}, so
  // it is -T^{-+} + T^{00} - T^{+-}. Every term sits at upper index zero
  // because the contracted pair contributes a + (-a) = 0 -- which is why the
  // sum is admissible at all, and why a contraction that paired its slots
  // wrongly would fail to compile rather than give a wrong answer.
  auto trace = Trace(tensor);
  static_assert(decltype(trace)::UpperIndex == 0);
  std::cout << "trace at a point " << (trace[2, 2]) << "\n"
            << "integral of the trace " << Integrate(trace) << "\n\n";

  // Symmetric and antisymmetric parts, as projections onto a symmetry group.
  auto sym = Symmetrise<Symmetric<2>>(tensor);
  auto skew = Symmetrise<Antisymmetric<2>>(tensor);
  const auto a = tensor.Component<0, 1>()[2, 2];
  const auto b = tensor.Component<1, 0>()[2, 2];
  std::cout << "symmetric part  " << (sym.Component<0, 1>()[2, 2])
            << "  expected " << 0.5 * (a + b) << "\n"
            << "skew part       " << (skew.Component<0, 1>()[2, 2])
            << "  expected " << 0.5 * (a - b) << "\n\n";

  // Materialise is where a lazy tensor stops being lazy. The symmetry is the
  // caller's to state: the product of two symmetric tensors is not symmetric,
  // and inferring one from an expression tree is a research problem. Asking
  // for it stores only the components that symmetry keeps.
  auto stored = Materialise<Symmetric<2>>(sym);
  std::cout << "materialised symmetric: " << decltype(stored)::StoredComponents
            << " components rather than " << Tensor::StoredComponents << "\n";

  // Higher ranks permute as asked -- on a grid that carries their upper
  // indices, which for rank 4 means up to 4.
  auto wider = Grid(16, 4);
  using Rank4 = TensorField<4, NoSymmetry<4>, ComplexTensor, Grid>;
  auto r = Rank4(wider);
  const auto& quartic = r;
  auto rotated = Permute<std::array<Int, 4>{2, 3, 0, 1}>(quartic);
  static_assert(decltype(rotated)::Rank == 4);
  std::cout << "rank-4 slot permutation compiles and relabels\n";

  FFTWpp::CleanUp();
}
