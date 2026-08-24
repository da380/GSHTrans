// 10 -- Real tensors
//
// A real tensor's components are not independent: T^{-alpha} = (-1)^N
// conj(T^alpha). Storing one component per orbit of that relation costs
// exactly the tensor's real degrees of freedom, and the rest are derived.

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <iostream>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  auto grid = Grid(16, 2);

  // The vector case, which the theory note spells out: u^0 is a real field
  // and u^- = -conj(u^+), so three reals a point rather than six.
  using Vector = TensorField<1, NoSymmetry<1>, RealTensor, Grid>;
  std::cout << "real vector: " << Vector::StoredComponents
            << " stored components, " << Vector::RealsPerPoint
            << " reals a point\n";

  auto v = Vector(grid);
  const auto fill = [&](auto&& u, Real tag) {
    using Node = std::remove_cvref_t<decltype(u)>;
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        if constexpr (std::same_as<typename Node::Value, RealValued>) {
          u[iTheta, iPhi] = tag + std::cos(0.3 * iTheta);
        } else {
          u[iTheta, iPhi] =
              Complex{tag + std::cos(0.3 * iTheta), std::sin(0.2 * iPhi)};
        }
      }
    }
  };
  fill(v.Component<-1>(), 1.0);
  fill(v.Component<0>(), 2.0);

  const auto& u = v;

  // The radial component is a real *field*, not a complex one that happens to
  // be real. It is at upper index zero, which is the only place the library
  // permits a real-valued field -- and the only place a self-paired component
  // can sit, since permutation preserves the slot sum and negation reverses
  // it.
  static_assert(std::same_as<decltype(u.Component<0>())::Value, RealValued>);
  static_assert(
      std::same_as<decltype(u.Component<-1>())::Value, ComplexValued>);

  // The derived component. conj reverses the upper index, which is why the
  // relation works out: a component stored at -1 derives one at +1.
  static_assert(decltype(u.Component<-1>())::UpperIndex == -1);
  static_assert(decltype(u.Component<1>())::UpperIndex == 1);
  std::cout << "u^-  " << (u.Component<-1>()[2, 2]) << "\n"
            << "u^+  " << (u.Component<1>()[2, 2]) << "   = -conj(u^-)\n\n";

  // A derived component is read-only: writing it would mean conjugating on
  // the way in, which a view cannot do. Write the one it derives from.
  static_assert(Vector::Writable<-1>);
  static_assert(!Vector::Writable<1>);

  // The counts, which are the independent real degrees of freedom at every
  // rank and symmetry, with nothing special-cased.
  std::cout
      << "reals per point\n"
      << "  rank 2                "
      << (TensorField<2, NoSymmetry<2>, RealTensor, Grid>::RealsPerPoint)
      << "   (3^2)\n"
      << "  rank 2 symmetric      "
      << (TensorField<2, Symmetric<2>, RealTensor, Grid>::RealsPerPoint)
      << "   (a real symmetric 3x3 matrix)\n"
      << "  rank 2 antisymmetric  "
      << (TensorField<2, Antisymmetric<2>, RealTensor, Grid>::RealsPerPoint)
      << "\n"
      << "  rank 3 symmetric      "
      << (TensorField<3, Symmetric<3>, RealTensor, Grid>::RealsPerPoint) << "\n"
      << "  rank 4                "
      << (TensorField<4, NoSymmetry<4>, RealTensor, Grid>::RealsPerPoint)
      << "  (3^4)\n"
      << "  rank 4 elastic        "
      << (TensorField<4, ElasticSymmetry, RealTensor, Grid>::RealsPerPoint)
      << "   (the elastic constants)\n";
}
