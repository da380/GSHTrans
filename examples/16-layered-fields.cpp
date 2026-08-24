// 16 -- Three-dimensional fields, and the gradient
//
// A field on a ball is a stack of angular fields, one per radius. Two
// dimensions stay the primitive: a slice of the stack *is* an ordinary spin
// field, so everything in examples 01 to 15 applies to it unchanged, and
// nothing here is a second version of anything there.
//
// What the radial axis adds is two things. It is the batch axis, so a whole
// stack transforms in one call rather than in a loop. And it is where the
// library stops: the radial derivative belongs to the application's
// discretisation and arrives as a callable.

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <span>
#include <vector>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;
  using Int = std::ptrdiff_t;

  constexpr auto lMax = Int{16};
  auto grid = Grid(lMax, 2);

  // The radial half: nodes, and optionally weights. Deliberately thin -- the
  // element basis, the differentiation matrices and any factorisation stay in
  // the application, because none of them is a spherical-harmonic question.
  const auto nR = Int{40};
  auto radii = std::vector<Real>(nR);
  auto weights = std::vector<Real>(nR);
  const auto h = (1.0 - 0.4) / (nR - 1);
  for (Int i = 0; i < nR; i++) {
    radii[i] = 0.4 + h * i;
    weights[i] = (i == 0 || i == nR - 1) ? h / 2 : h;
  }
  auto radial = RadialGrid<Real>(radii, weights);

  std::cout << std::scientific << std::setprecision(3);

  //------------------------------------------------------------------------//
  // A stack, and its slices
  //------------------------------------------------------------------------//

  auto stack = LayeredSpinField<0, Grid>(radial, grid);
  std::cout << "radii " << stack.NumberOfRadii() << ", angular points "
            << stack.FieldSize() << ", total " << stack.Size() << "\n";

  // Slice(i) is a view over the stack's own storage, so writing through it
  // writes the stack. It is an ordinary phase-1 node: the index algebra and
  // the lazy evaluation lift with no extra machinery.
  for (auto i : stack.RadiusIndices()) {
    const auto r = radial.Radius(i);
    auto slice = stack.Slice(i);
    for (Int j = 0; j < slice.Size(); j++) slice.Data()[j] = r * r;
  }

  // Integrating over radius with the grid's own weights. The r^2 of the
  // volume element is *not* inserted here: fold it into the weights if that
  // is what you want, because guessing would be wrong half the time.
  const auto overRadius = IntegrateRadially(stack);
  std::cout << "int_{0.4}^{1} r^2 dr = " << overRadius[0].real() << "  (exact "
            << (1.0 - 0.4 * 0.4 * 0.4) / 3 << ")\n\n";

  //------------------------------------------------------------------------//
  // The radial axis is the batch axis
  //------------------------------------------------------------------------//

  // One call, not a loop over radii. The Wigner block for this upper index is
  // streamed once for the whole stack rather than once per radius, and the
  // result is bit-for-bit what the loop would have given: batching is a
  // scheduling change and not a numerical one.
  auto expansion = Expand(stack, lMax, Execution::Parallel());
  std::cout << "coefficients per radius: " << expansion.CoefficientSize()
            << ", stacked: " << expansion.Size() << "\n\n";

  //------------------------------------------------------------------------//
  // The seam
  //------------------------------------------------------------------------//

  // A radial operator maps one radial line to another. This library never
  // supplies one: a finite-difference derivative is banded, a spectral-element
  // one is block-diagonal, and a solve may carry a factorisation you want to
  // reuse. Here, a three-point centred difference.
  const auto ddr = [h](std::span<const Complex> in, std::span<Complex> out) {
    const auto n = static_cast<Int>(in.size());
    for (Int i = 0; i < n; i++) {
      if (i == 0) {
        out[i] = (in[1] - in[0]) / h;
      } else if (i == n - 1) {
        out[i] = (in[n - 1] - in[n - 2]) / h;
      } else {
        out[i] = (in[i + 1] - in[i - 1]) / (2 * h);
      }
    }
  };

  // It applies on either side of the transform, because it acts along r alone
  // and so commutes with the angular one. The spectral side is usually where
  // you want it: a coefficient buffer is complex whatever the field's reality,
  // so the operator never has to be generic over Real and Complex.
  auto slope = ApplyRadially(expansion, ddr, Execution::Parallel());

  // The stack holds r^2 at every angular point, so the ratio of the two
  // (l, m) = (0, 0) coefficients is (d/dr r^2) / r^2 = 2/r, whatever the
  // normalisation of Y_00 happens to be.
  {
    const auto r = radial.Radius(20);
    const auto ratio = (slope[20, 0, 0] / expansion[20, 0, 0]).real();
    std::cout << "d/dr of r^2 at r = " << r << ": " << ratio * r * r
              << "  (exact " << 2 * r << ")\n\n";
  }

  //------------------------------------------------------------------------//
  // The gradient
  //------------------------------------------------------------------------//

  // D&T's gradient in the canonical basis is grad = e_0 d_r + r^{-1} grad_1,
  // and the two halves are on opposite sides of the seam. grad_1 is example
  // 15's surface gradient, unchanged; d_r is yours.
  //
  //     sigma = +-1 :  r^{-1} (grad_1 f)^sigma
  //     sigma =   0 :  df/dr
  //
  // A scalar is a rank-0 tensor and its gradient is a rank-1 one, so the
  // answer has the three components a vector field on the ball has.
  constexpr auto l = Int{6};
  constexpr auto m = Int{3};
  const auto a = Real{3};

  auto f = LayeredScalarExpansion<Grid, ComplexTensor>(radial, grid, lMax);
  for (auto i : radial.RadiusIndices()) {
    f.ComponentStack<>()[i, l, m] = std::pow(radial.Radius(i), a);
  }

  // An exact d/dr for this particular field, so that what is printed below is
  // the algebra and not a difference formula's truncation error.
  const auto exact = [&](Real power) {
    return
        [power, &radial](std::span<const Complex> in, std::span<Complex> out) {
          for (std::size_t i = 0; i < in.size(); i++) {
            out[i] = power * in[i] / radial.Radius(static_cast<Int>(i));
          }
        };
  };

  auto grad = Gradient(f, exact(a));
  static_assert(decltype(grad)::Rank == 1);

  const auto r = radial.Radius(20);
  const auto omega = std::sqrt(l * (l + 1.0) / 2);
  std::cout << "at r = " << r << ", for f = r^3 Y_{" << l << "," << m << "}:\n";
  std::cout << "  (grad f)^0  = " << (grad.Coefficient<0>(20, l, m)).real()
            << "   expect " << a * std::pow(r, a - 1) << "\n";
  std::cout << "  (grad f)^+1 = " << (grad.Coefficient<1>(20, l, m)).real()
            << "   expect " << omega * std::pow(r, a - 1) << "\n\n";

  //------------------------------------------------------------------------//
  // Why the connection terms matter
  //------------------------------------------------------------------------//

  // Take the gradient again. This time the operand is a vector whose e_0
  // component is not zero, so grad_1 has to differentiate the basis as well:
  // the terms that move a slot between e_0 and e_{+-} are what make the answer
  // right. Contracting with the metric g_{ab} = (-1)^a delta_{a+b,0} must give
  // the Laplacian, and there is nowhere for an error in those terms to hide.
  auto second = Gradient(grad, exact(a - 1));
  static_assert(decltype(second)::Rank == 2);

  const auto trace = second.Coefficient<0, 0>(20, l, m) -
                     second.Coefficient<1, -1>(20, l, m) -
                     second.Coefficient<-1, 1>(20, l, m);
  std::cout << "  lap f = " << trace.real() << "   expect "
            << (a * (a + 1) - l * (l + 1.0)) * std::pow(r, a - 2) << "\n";

  // The gradient carries an explicit r^{-1}, so it is not defined at r = 0.
  // That is a singularity of the basis, not of the field, and the library says
  // so rather than returning an infinity.
}
