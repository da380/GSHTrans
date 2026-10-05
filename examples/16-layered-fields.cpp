// 16 -- Three-dimensional fields, and the gradient
//
// What this shows. A field on a ball or spherical shell is a stack of angular
// fields, one per radius, stored radius-major: nR slices laid end to end, each
// an ordinary angular field. Two dimensions stay the primitive: a slice of the
// stack *is* a spin field (a view over the stack's storage), so everything in
// examples 01 to 15 applies to it unchanged, and nothing here is a second
// version of anything there.
//
// What the radial axis adds is two things. It is the batch axis (example 07),
// so a whole stack transforms in one call rather than in a loop. And it is
// the seam between the library and the radial discretisation: a radial
// operator, such as d/dr, is any callable mapping one radial line to another,
// and the library applies it along the radial axis without needing to know
// how it was built. Example 18 shows the ready-made operators and the rules
// for writing one; here they are written by hand.
//
// Read first. 07 (batches), 12 (expansions), 15 (the surface gradient).
//
// Introduced. RadialGrid; LayeredSpinField and its Slice(i);
// IntegrateRadially; Expand on a stack; ApplyRadially; LayeredScalarExpansion
// and Gradient, the full three-dimensional gradient
// grad = e_0 d/dr + r^{-1} grad_1.
//
// The output. The stack's sizes; a trapezoidal radial integral against its
// exact value; the stacked coefficient count; a finite-difference d/dr; and
// the components of grad f and the trace of grad grad f for f = r^3 Y_{6,3},
// each against its closed form.

#include <GSHTrans/GSHTrans.hpp>
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

  // The radial half: nodes in increasing order, and optionally quadrature
  // weights. Here 40 equally spaced radii on [0.4, 1] with trapezoidal
  // weights. A RadialGrid is a value-semantic handle like the angular grid:
  // copies share one implementation, and two stacks are on the same radial
  // grid when their handles agree. It can also record an element partition
  // (example 19); the differentiation itself stays outside it.
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

  // FieldSize() is the number of angular samples in one slice and Size() the
  // whole stack, nR * FieldSize(). Slice(i) is a view over the stack's own
  // storage at radius index i, so writing through it writes the stack. It is
  // an ordinary spin-weighted node, so the index algebra and the lazy
  // evaluation of examples 03-05 apply to it with no extra machinery. Here
  // every point of slice i is set to r_i^2.
  for (auto i : stack.RadiusIndices()) {
    const auto r = radial.Radius(i);
    auto slice = stack.Slice(i);
    for (Int j = 0; j < slice.Size(); j++) slice.Data()[j] = r * r;
  }

  // Integrating over radius with the grid's own weights leaves one value per
  // angular point. The r^2 of the volume element is *not* inserted: fold it
  // into the weights to integrate over the ball, and leave it out to
  // integrate a radial profile. Here the integrand is already r^2, so the
  // result is int r^2 dr; the difference from the exact value printed is the
  // trapezoidal rule's error, not the library's.
  const auto overRadius = IntegrateRadially(stack);
  std::cout << "int_{0.4}^{1} r^2 dr = " << overRadius[0].real() << "  (exact "
            << (1.0 - 0.4 * 0.4 * 0.4) / 3 << ")\n\n";

  //------------------------------------------------------------------------//
  // The radial axis is the batch axis
  //------------------------------------------------------------------------//

  // One call, not a loop over radii. Every slice shares grid, degree and
  // upper index, which is what a batch requires, so the Wigner values for
  // this upper index are read once for the whole stack rather than once per
  // radius. Batching widens the inner loop without reordering any sum, so the
  // coefficients are the ones a loop over radii would give
  // (tests/TestLayered.cpp checks they are equal). Threading is opt-in, as
  // everywhere: Execution::Parallel() asks for it.
  //
  // The result is a layered expansion: one block of CoefficientSize()
  // coefficients per radius, radius-major, indexed [i, l, m].
  auto expansion = Expand(stack, lMax, Execution::Parallel());
  std::cout << "coefficients per radius: " << expansion.CoefficientSize()
            << ", stacked: " << expansion.Size() << "\n\n";

  //------------------------------------------------------------------------//
  // The seam
  //------------------------------------------------------------------------//

  // A radial operator maps one radial line -- the nR values at one angular
  // point, or of one (l, m) coefficient -- to another. It is a callable
  // rather than a matrix because discretisations differ: a finite-difference
  // derivative is banded, a spectral-element one is block-diagonal, and a
  // solve may carry a factorisation you want to reuse. Ready-made ones exist
  // (example 18); this one is written out to show the shape: a three-point
  // centred difference, one-sided at the two ends.
  //
  // It is called once per line, possibly from several threads at once, so
  // it must not modify shared state. The in and out spans are always
  // distinct, so it may read in after writing out.
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

  // ApplyRadially gathers each line, applies the operator and scatters the
  // result into a new stack of the same shape. It works on either side of the
  // transform, because the operator acts along r alone and so commutes with
  // the angular transform. The spectral side is usually where you want it: a
  // coefficient buffer is complex whatever the field's reality, so the
  // operator only ever sees std::complex data, as this one assumes.
  auto slope = ApplyRadially(expansion, ddr, Execution::Parallel());

  // The stack holds r^2 at every angular point, so only (l, m) = (0, 0) is
  // non-zero and the ratio of the two (0, 0) coefficients at a radius is
  // (d/dr r^2) / r^2 = 2/r, whatever the normalisation of Y_00 happens to be.
  // A centred difference is exact on a quadratic, so the match is to
  // rounding at an interior radius.
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
  // 15's surface gradient, unchanged; d_r is the radial operator passed in.
  // For a scalar f,
  //
  //     sigma = +-1 :  r^{-1} (grad_1 f)^sigma
  //     sigma =   0 :  df/dr
  //
  // A scalar is a rank-0 tensor and its gradient is a rank-1 one, so the
  // answer has the three components a vector field on the ball has. The
  // whole computation is spectral: the input is a layered tensor expansion,
  // here a scalar one, filled with f = r^3 at (l, m) = (6, 3) for every
  // radius. ComponentStack<>() is its single component, indexed [i, l, m].
  constexpr auto l = Int{6};
  constexpr auto m = Int{3};
  const auto a = Real{3};

  auto f = LayeredScalarExpansion<Grid, ComplexTensor>(radial, grid, lMax);
  for (auto i : radial.RadiusIndices()) {
    f.ComponentStack<>()[i, l, m] = std::pow(radial.Radius(i), a);
  }

  // An exact d/dr for this particular field, so that what is printed below is
  // the algebra and not a difference formula's truncation error: on a line
  // holding c r^p, d/dr is p/r times the line. exact(p) returns that
  // operator for a given power.
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

  // The expected values: (grad f)^0 = d/dr r^3 = 3 r^2, and
  // (grad f)^{+1} = r^{-1} Omega^0_l r^3 = Omega r^2 with
  // Omega = sqrt(l(l+1)/2), the factor of example 15. Coefficient<sigma>(i,
  // l, m) reads a component at radius index i.

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
  // Every component of grad f goes as r^2, so the radial operator this time
  // is exact(a - 1). The expected value is the Laplacian of r^a Y_lm,
  // [a(a+1) - l(l+1)] r^{a-2}.
  auto second = Gradient(grad, exact(a - 1));
  static_assert(decltype(second)::Rank == 2);

  const auto trace = second.Coefficient<0, 0>(20, l, m) -
                     second.Coefficient<1, -1>(20, l, m) -
                     second.Coefficient<-1, 1>(20, l, m);
  std::cout << "  lap f = " << trace.real() << "   expect "
            << (a * (a + 1) - l * (l + 1.0)) * std::pow(r, a - 2) << "\n";

  // The gradient carries an explicit r^{-1}, so it is not defined at r = 0.
  // That is a singularity of the basis, not of the field: Gradient throws
  // std::invalid_argument if the radial grid's first radius is not positive,
  // rather than returning an infinity. This grid starts at 0.4, clear of it.
}
