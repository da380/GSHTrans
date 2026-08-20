// A tensor field on the sphere: one buffer, components as spin fields.

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  constexpr auto lMax = std::ptrdiff_t{16};

  // A rank-2 tensor has components at upper index +-2, so the grid has to
  // carry them. Asking for a narrower grid is an error at construction.
  auto grid = Grid(lMax, 2, FFTWpp::Measure);

  // Symmetric, so six of the nine components are stored and the other three
  // are the same fields seen through transposed indices.
  using Tensor = TensorField<2, Symmetric<2>, ComplexTensor, Grid>;
  auto t = Tensor(grid);

  std::cout << "stored components: " << Tensor::StoredComponents << " of "
            << Tensor::Components << "\n";

  // Components are addressed by multi-index and come back as ordinary
  // spin-weighted fields, at the upper index the multi-index implies.
  static_assert(decltype(t.Component<1, 1>())::UpperIndex == 2);
  static_assert(decltype(t.Component<-1, 1>())::UpperIndex == 0);

  auto fill = [&](auto&& u, Real phase) {
    for (auto iTheta : grid.CoLatitudeIndices()) {
      for (auto iPhi : grid.LongitudeIndices()) {
        const auto [theta, phi] =
            std::pair(grid.CoLatitudes()[iTheta], grid.Longitudes()[iPhi]);
        u[iTheta, iPhi] = Complex{std::cos(theta + phase) * std::cos(phi),
                                  std::sin(theta) * std::sin(phi + phase)};
      }
    }
  };
  fill(t.Component<-1, -1>(), 0.0);
  fill(t.Component<-1, 0>(), 0.3);
  fill(t.Component<-1, 1>(), 0.6);
  fill(t.Component<0, 0>(), 0.9);
  fill(t.Component<0, 1>(), 1.2);
  fill(t.Component<1, 1>(), 1.5);

  // Symmetry is storage, not bookkeeping: the transposed component is the
  // same field, so writing one wrote both.
  std::cout << "T^{01} == T^{10}: "
            << ((t.Component<0, 1>()[3, 4] == t.Component<1, 0>()[3, 4]) ? "yes"
                                                                        : "no")
            << "\n";

  // The trace, contracted with the metric of the theory note: with
  // g_{ab} = (-1)^a delta_{a+b,0} it is -T^{-+} + T^{00} - T^{+-}. Every term
  // is at upper index zero, so the sum is too, and only then can it be
  // integrated over the sphere. This is a lazy expression: nothing is
  // evaluated until Integrate walks it.
  const auto& tensor = t;
  auto trace = -tensor.Component<-1, 1>() + tensor.Component<0, 0>() -
               tensor.Component<1, -1>();
  static_assert(decltype(trace)::UpperIndex == 0);

  std::cout << "mean trace  = " << Integrate(trace) / (4 * std::numbers::pi)
            << "\n";
  std::cout << "trace norm  = " << std::sqrt(Integrate(abs2(trace))) << "\n";

  // The whole tensor through the spectral domain and back. Components sharing
  // an upper index are transformed as one batch, which is what the buffer's
  // component order is arranged for, and the threading is asked for
  // explicitly -- the library never creates threads on its own.
  auto coefficients = std::vector<Complex>(t.CoefficientSize(lMax));
  std::cout << "coefficients: " << coefficients.size() << " for "
            << Tensor::StoredComponents << " components\n";

  // Transforming an arbitrary field is a *projection*, not a round trip: a
  // component at upper index N has no content below degree |N|, since no
  // generalized spherical harmonic there exists, so whatever the fill above
  // put at low degree is discarded. Projecting once makes the field
  // band-limited; after that the round trip is an identity.
  t.ForwardTransformation(lMax, coefficients, Execution::Parallel(4));
  t.InverseTransformation(lMax, coefficients, Execution::Parallel(4));

  const auto before = t.Component<0, 1>()[3, 4];
  t.ForwardTransformation(lMax, coefficients, Execution::Parallel(4));
  t.InverseTransformation(lMax, coefficients, Execution::Parallel(4));
  const auto after = t.Component<0, 1>()[3, 4];

  std::cout << "round trip drift, once band-limited: "
            << std::abs(after - before) << "\n";

  FFTWpp::CleanUp();
}
