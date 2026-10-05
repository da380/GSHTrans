// 18 -- The radial seam, and the operators that come ready made
//
// What this shows. Example 16 introduced the seam: the library applies
// whatever callable it is given along the radial axis, and does not own the
// radial discretisation. That is deliberate -- a finite-difference derivative
// is banded, a spectral-element one is block-diagonal, and a caller with a
// factorisation wants to apply it rather than hand over a matrix.
//
// A few radial derivatives are supplied (GSHTrans/Layered/
// RadialDerivatives.hpp and RadialSplineDerivative.hpp), so a gradient can be
// taken without writing one first. They are conveniences, not the library's
// claim on the radial half: production codes use a finite-element basis,
// finite differences or a radial spectral basis -- Chebyshev, say, which
// would sit behind this same seam as a transform, a multiply and a transform
// back. They are also worth reading as worked examples of the two rules any
// radial operator must follow.
//
// The example then shows the second layout a radial code needs. A layered
// field is stored radius-major, [r][(l, m)], which is what the angular
// transform wants. A radial solve wants [(l, m)][r], where each radial line
// is a contiguous vector; ExpandToLines, ApplyToLines and EvaluateLines work
// in that layout directly.
//
// Read first. 07 (threading policy), 16 (layered fields, ApplyRadially,
// Gradient).
//
// Introduced. FiniteDifferenceDerivative, LagrangeDerivative and, when the
// Interpolation dependency is present, SplineDerivative; the RadialOperator
// contract; ExpandToLines, ApplyToLines, EvaluateLines.
//
// The output. Errors in d/dr of r^3 for each operator on an uneven grid;
// the Laplacian of r^2 Y_3^1 from two applications of Gradient; a
// hand-written operator applied through ApplyRadially, exact to rounding; and
// the lines route agreeing exactly with the gather route.

#include <GSHTrans/GSHTrans.hpp>
#include <algorithm>
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

  // Six radii, deliberately unevenly spaced: every operator here takes the
  // spacing as it finds it, and a rule that quietly assumed a uniform grid
  // would pass a test on a uniform one. No weights are given, since nothing
  // here integrates.
  const auto radii = std::vector<Real>{0.40, 0.55, 0.70, 1.00, 1.30, 1.45};
  auto radial = RadialGrid<Real>(radii);

  std::cout << std::scientific << std::setprecision(3);

  // Apply an operator to one line by hand, which is all a radial operator
  // is: a const call taking an input span and an output span of nR values.
  const auto line = [&](const auto& op, const std::vector<Real>& in) {
    auto out = std::vector<Real>(in.size());
    op(std::span<const Real>(in), std::span<Real>(out));
    return out;
  };

  //------------------------------------------------------------------------//
  // Three operators, and what distinguishes them
  //------------------------------------------------------------------------//

  // Finite differences with Fornberg's weights: banded, exact for
  // polynomials of degree up to the order asked for (here 4, a five-point
  // stencil), centred where there is room and one-sided at the ends -- which
  // is what makes it usable on a grid that *has* ends, those being exactly
  // the radii a boundary condition is applied at. Its cost per line grows
  // with nR only linearly. The one to reach for by default.
  const auto fd = FiniteDifferenceDerivative<Real>(radial, 4);

  // The differentiation matrix of the nodes: d/dr of the polynomial of
  // degree nR - 1 through all of them. Exact to the highest degree any
  // operator on these nodes could be, and the right thing on a few nodes --
  // the Gauss-Lobatto points of one element, say. Wrong on many: the cost is
  // nR^2 per line, and a global polynomial through many equally spaced
  // points diverges.
  //
  // Both are exact on a cubic, so the errors printed are rounding.
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
  // The cubic spline: global like the matrix but linear in the number of
  // radii, and unlike a global polynomial it does not fall apart as the nodes
  // multiply. Built on Interpolation's factorised spline system, so it exists
  // only when that optional dependency does (GSHTRANS_HAVE_INTERPOLATION).
  //
  // Its end conditions matter. Natural, the default, forces the second
  // derivative to zero at the ends, which is wrong for a curve whose is not;
  // not-a-knot stays fourth order right up to them, and reproduces a cubic
  // exactly.
  const auto natural = SplineDerivative<Real>(radial);
  const auto notAKnot = SplineDerivative<Real>(
      radial, BoundaryCondition::NotAKnot, BoundaryCondition::NotAKnot);

  const auto bySplineN = line(natural, cubic);
  const auto bySplineK = line(notAKnot, cubic);
  std::cout << "\nthe same, by spline\n"
            << "     r         natural         not-a-knot\n";
  for (std::size_t i = 0; i < radii.size(); i++) {
    const auto exact = 3 * radii[i] * radii[i];
    std::cout << "  " << radii[i] << "      " << std::abs(bySplineN[i] - exact)
              << "      " << std::abs(bySplineK[i] - exact) << '\n';
  }
  std::cout
      << "  (natural is wrong at the ends by construction; the cubic's\n"
         "   second derivative is not zero there and natural says it is)\n";
#endif

  //------------------------------------------------------------------------//
  // Writing your own, which is the general case
  //------------------------------------------------------------------------//

  // A radial operator is any callable taking one line to another, as
  // op(std::span<const T> in, std::span<T> out), both of length nR and never
  // the same span. Two obligations come with that, and both follow from how
  // the library calls it: once per line, from inside a parallel region,
  // through a const reference.
  //
  //   - it must be safe to call concurrently, so any scratch is thread_local
  //     or local to the call, and never a mutable member;
  //   - it is called SliceSize() times per application -- once per angular
  //     point, or once per (l, m) coefficient, tens of thousands at
  //     production degree -- so whatever depends only on the nodes is
  //     computed once, at construction, and an allocation inside the call is
  //     a defect.
  //
  // Here is one that obeys both: the exact derivative of c * r^a, which is
  // what a test wants when it needs no truncation error. It is example 16's
  // lambda written as a type, with the radii copied in at construction. It
  // is applied below, once there is a field to apply it to.
  struct PowerDerivative {
    Real power;
    std::vector<Real> r;

    void operator()(std::span<const Complex> in, std::span<Complex> out) const {
      for (std::size_t i = 0; i < in.size(); i++) out[i] = power * in[i] / r[i];
    }
  };

  //------------------------------------------------------------------------//
  // And they compose with everything else
  //------------------------------------------------------------------------//

  // ApplyRadially applies any of them along the radial axis of a stack, in
  // either domain: a field's line is the samples at one angular point and an
  // expansion's is one coefficient across radius, and the radial axis does
  // not distinguish them. The field here is complex, with arbitrary values
  // written straight into its buffer.
  auto f = LayeredSpinField<0, Grid, ComplexValued>(radial, grid);
  for (Int i = 0; i < f.NumberOfRadii() * f.SliceSize(); i++) {
    f.Data()[i] = Complex{std::cos(0.05 * i), std::sin(0.11 * i)};
  }
  const auto df = ApplyRadially(f, fd, Execution::Parallel());
  std::cout << "\nApplyRadially over " << f.SliceSize()
            << " lines, threaded: " << df.NumberOfRadii() << " radii out\n";

  // The full gradient takes a radial operator and supplies the angular part
  // (example 16). Production code often works with the pieces instead:
  // expand, apply the radial operator to the coefficients, evaluate back.
  // The pieces are public on purpose, and Gradient is the assembled
  // convenience over them.
  constexpr auto l = Int{3};
  const auto a = Real{2};
  auto scalar = LayeredScalarExpansion<Grid, ComplexTensor>(radial, grid, lMax);
  for (auto i : radial.RadiusIndices()) {
    scalar.ComponentStack<>()[i, l, 1] = std::pow(radial.Radius(i), a);
  }

  const auto gradient = Gradient(scalar, fd);
  const auto second = Gradient(gradient, fd);
  // The metric trace of the second gradient is the Laplacian, as in
  // example 16. The power 2 is written into the label rather than printed
  // from `a`, which would come out in scientific notation under the format
  // set above.
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

  // The hand-written operator goes through the same seam as the ready-made
  // ones. Applied to r^a Y_l^1 it gives a r^(a-1) Y_l^1, exactly.
  auto power = LayeredSpinExpansion<0, Grid>(radial, grid, lMax);
  for (auto i : radial.RadiusIndices()) {
    power[i, l, 1] = std::pow(radial.Radius(i), a);
  }
  const auto exact = ApplyRadially(power, PowerDerivative{a, radii});
  auto worstPower = Real{0};
  for (auto i : radial.RadiusIndices()) {
    const auto expected = a * std::pow(radial.Radius(i), a - 1);
    worstPower =
        std::max(worstPower, std::abs(Complex{exact[i, l, 1]} - expected));
  }
  std::cout << "\nPowerDerivative through ApplyRadially, worst error "
            << worstPower << '\n';

  //------------------------------------------------------------------------//
  // The same, with a radial line contiguous
  //------------------------------------------------------------------------//

  // ApplyRadially gathers each line out of radius-major storage, applies the
  // operator and scatters the answer back, which is right when a line is
  // touched once. A solver touches it many times -- an iteration, a
  // factorisation applied again and again -- and then the layout to be in is
  // the other one, [(l, m)][r], where a line is a contiguous vector. The
  // RadialMajor buffer holds that layout: a shape and the radial grid, with
  // no angular grid or upper index, since nothing angular is meaningful
  // once the data are cut into radial lines.
  //
  // ExpandToLines transforms straight into that layout and EvaluateLines out
  // of it, so the radius-major coefficients never exist: one whole set of
  // coefficients saved, hundreds of megabytes to gigabytes at production
  // sizes. RadialMajor(Expand(f)) gives the same numbers, bit for bit, by way
  // of the copy. What this saves is memory; whether it also saves time
  // depends on the machine.
  auto lines = ExpandToLines(f, lMax, Execution::Parallel());
  std::cout << "\nExpandToLines: " << lines.NumberOfLines() << " lines of "
            << lines.NumberOfRadii() << " radii, each contiguous\n";

  // In place, as often as is wanted: here d/dr twice, giving d^2/dr^2. The
  // operator is handed a scratch line when the two buffers are one, so it
  // need not cope with its output being its input. Then back to a layered
  // field on the angular grid.
  ApplyToLines(lines, lines, fd, Execution::Parallel());
  ApplyToLines(lines, lines, fd, Execution::Parallel());
  const auto curvature =
      EvaluateLines(lines, grid, lMax, Execution::Parallel());

  // Against the gathering route, which must agree exactly: the arithmetic is
  // the same and only the copying moved. Under the same policy, that is. A
  // threaded transform chunks its sums differently from a sequential one and
  // may differ from it in the last bit, as any two orders of summation may,
  // so the promise is between the two routes and not between two policies.
  // The example exits non-zero if they differ at all.
  const auto expanded = Expand(f, lMax, Execution::Parallel());
  const auto gathered = Evaluate(ApplyRadially(ApplyRadially(expanded, fd), fd),
                                 Execution::Parallel());
  auto worst = Real{0};
  for (Int i = 0; i < curvature.Size(); i++) {
    worst = std::max(worst, std::abs(curvature.Data()[i] - gathered.Data()[i]));
  }
  std::cout << "  d^2/dr^2 through the lines against through the gather: "
            << "differ by " << worst << '\n';
  if (worst != 0) return 1;

  return 0;
}
