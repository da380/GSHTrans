// 18 -- The radial seam, and the operators that come ready made
//
// Example 16 showed the seam: the library applies whatever callable it is
// given along the radial axis and owns no discretisation. That is still the
// arrangement, and it is the point rather than a gap -- a finite-difference
// derivative is banded, a spectral-element one is block-diagonal, and a caller
// with a factorisation wants to apply it rather than hand over a matrix.
//
// What this example adds is that a few operators now come with the library, so
// nothing has to be written before a gradient will run. They are also worth
// reading: each is a worked example of the two obligations the seam states,
// for a caller whose discretisation is none of them.
//
// The library is not trying to own the radial half. Production codes use a
// finite-element basis, finite differences, or a radial spectral basis --
// Chebyshev, say, which would sit behind this same seam as a transform, a
// multiply and a transform back. These are conveniences.

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

  constexpr auto lMax = Int{8};
  auto grid = Grid(lMax, 2);

  // Deliberately uneven: every operator here takes the spacing as it finds it,
  // and a rule that quietly assumed a uniform grid would pass on one.
  const auto radii = std::vector<Real>{0.40, 0.55, 0.70, 1.00, 1.30, 1.45};
  auto radial = RadialGrid<Real>(radii);

  std::cout << std::scientific << std::setprecision(3);

  const auto line = [&](const auto& op, const std::vector<Real>& in) {
    auto out = std::vector<Real>(in.size());
    op(std::span<const Real>(in), std::span<Real>(out));
    return out;
  };

  //------------------------------------------------------------------------//
  // Three operators, and what distinguishes them
  //------------------------------------------------------------------------//

  // Finite differences: banded, exact for polynomials up to the order asked
  // for, centred where there is room and one-sided at the ends -- which is
  // what makes it usable on a grid that *has* ends, those being exactly the
  // radii a boundary condition is applied at. The one to reach for by
  // default.
  const auto fd = FiniteDifferenceDerivative<Real>(radial, 4);

  // The differentiation matrix of the nodes: exact to the highest degree any
  // operator on them could be, and the right thing on a few nodes -- the
  // Gauss-Lobatto points of one element, say. Wrong on many: the cost is
  // quadratic and a global polynomial through many points diverges.
  const auto lagrange = LagrangeDerivative<Real>(radial);

  auto cubic = std::vector<Real>{};
  for (auto r : radii) cubic.push_back(r * r * r);

  const auto byDifference = line(fd, cubic);
  const auto byMatrix = line(lagrange, cubic);

  std::cout << "d/dr of r^3, against the exact 3r^2\n"
            << "     r      finite difference     matrix\n";
  for (std::size_t i = 0; i < radii.size(); i++) {
    const auto exact = 3 * radii[i] * radii[i];
    std::cout << "  " << radii[i] << "      "
              << std::abs(byDifference[i] - exact) << "      "
              << std::abs(byMatrix[i] - exact) << '\n';
  }

#ifdef GSHTRANS_HAVE_INTERPOLATION
  // The spline: global like the matrix but linear in the number of radii, and
  // unlike a global polynomial it does not fall apart as the nodes multiply.
  // Built on Interpolation's factorised system, so it exists only when that
  // optional dependency does.
  //
  // Its end conditions matter. Natural forces the second derivative to zero at
  // the ends, which is wrong for a curve whose is not; not-a-knot stays fourth
  // order right up to them.
  const auto natural = SplineDerivative<Real>(radial);
  const auto notAKnot = SplineDerivative<Real>(
      radial, BoundaryCondition::NotAKnot, BoundaryCondition::NotAKnot);

  const auto bySplineN = line(natural, cubic);
  const auto bySplineK = line(notAKnot, cubic);
  std::cout << "\nthe same, by spline\n"
            << "     r         natural         not-a-knot\n";
  for (std::size_t i = 0; i < radii.size(); i++) {
    const auto exact = 3 * radii[i] * radii[i];
    std::cout << "  " << radii[i] << "      "
              << std::abs(bySplineN[i] - exact) << "      "
              << std::abs(bySplineK[i] - exact) << '\n';
  }
  std::cout << "  (natural is wrong at the ends by construction; the cubic's\n"
               "   second derivative is not zero there and natural says it is)\n";
#endif

  //------------------------------------------------------------------------//
  // Writing your own, which is the general case
  //------------------------------------------------------------------------//

  // A radial operator is any callable taking one line to another. Two
  // obligations come with that, and both follow from how the library calls it:
  // once per line, from inside a parallel region, through a const reference.
  //
  //   - it must be safe to call concurrently, so any scratch is thread_local
  //     and never a mutable member;
  //   - it is called SliceSize() times per application -- tens of thousands --
  //     so whatever depends only on the nodes is computed once, at
  //     construction, and an allocation inside the call is a defect.
  //
  // Here is one that obeys both, in a dozen lines: the exact derivative of
  // c * r^a, which is what a test wants when it needs no truncation error.
  struct PowerDerivative {
    Real power;
    std::vector<Real> r;

    void operator()(std::span<const Complex> in,
                    std::span<Complex> out) const {
      for (std::size_t i = 0; i < in.size(); i++) out[i] = power * in[i] / r[i];
    }
  };

  //------------------------------------------------------------------------//
  // And they compose with everything else
  //------------------------------------------------------------------------//

  // The seam applies any of them along the radial axis of a stack, in either
  // domain: a field's line is angular samples and an expansion's is
  // coefficients, and the radial axis does not distinguish them.
  auto f = LayeredSpinField<0, Grid, ComplexValued>(radial, grid);
  for (Int i = 0; i < f.NumberOfRadii() * f.SliceSize(); i++) {
    f.Data()[i] = Complex{std::cos(0.05 * i), std::sin(0.11 * i)};
  }
  const auto df = ApplyRadially(f, fd, Execution::Parallel());
  std::cout << "\nApplyRadially over " << f.SliceSize()
            << " lines, threaded: " << df.NumberOfRadii() << " radii out\n";

  // The full gradient takes one and supplies the rest. Note what production
  // code often does instead: expand, apply the radial operator to the
  // coefficients, evaluate back. The pieces are separate on purpose, and
  // Gradient is the assembled convenience over them.
  constexpr auto l = Int{3};
  const auto a = Real{2};
  auto scalar = LayeredScalarExpansion<Grid, ComplexTensor>(radial, grid, lMax);
  for (auto i : radial.RadiusIndices()) {
    scalar.ComponentStack<>()[i, l, 1] = std::pow(radial.Radius(i), a);
  }

  const auto gradient = Gradient(scalar, fd);
  const auto second = Gradient(gradient, fd);
  // `a` would print in scientific notation under the format set above, and a
  // power of two is easier to read as a two.
  std::cout << "\nthe Laplacian of r^2 Y_" << l << "^1, "
            << "against [a(a+1) - l(l+1)] r^(a-2)\n";
  for (auto i : radial.RadiusIndices()) {
    const auto trace = second.Coefficient<0, 0>(i, l, 1) -
                       second.Coefficient<1, -1>(i, l, 1) -
                       second.Coefficient<-1, 1>(i, l, 1);
    const auto expected =
        (a * (a + 1) - l * (l + 1.0)) * std::pow(radial.Radius(i), a - 2);
    std::cout << "  r = " << radial.Radius(i) << "   error "
              << std::abs(trace.real() - expected) << '\n';
  }
  std::cout << "  (r^2 is degree two, so a five-point rule is exact on it and\n"
               "   the error is rounding rather than truncation)\n";

  return 0;
}
