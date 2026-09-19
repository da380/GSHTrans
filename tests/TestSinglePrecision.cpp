#include <gtest/gtest.h>

#include <GSHTrans/All>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

// Single precision, one test a layer, and a grid covering only the
// non-negative upper indices.
//
// Neither was instantiated anywhere -- not in a test, an example or a
// benchmark -- although both are offered: `float` by every concept in the
// library, and NRange = NonNegative by dedicated branches in the grid and in
// the Wigner tables. They worked when first tried, so nothing here was written
// after a failure. It is here so that they go on working: a template nobody
// instantiates is checked by nobody, and these are the cheapest instantiations
// that reach each layer's arithmetic.

using namespace GSHTrans;

namespace {

using Int = std::ptrdiff_t;
using Real = float;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

constexpr Int lMax = 16;
constexpr Real tolerance = 2e-5f;

}  // namespace

TEST(SinglePrecision, ASpinFieldRoundTripsAndItsAlgebraHolds) {
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  auto e = SpinExpansion<2, Grid>(grid, lMax);
  for (auto l : e.Degrees()) {
    for (auto m : e.Orders(l)) {
      e[l, m] = Complex{std::sin(0.3f * static_cast<Real>(l + m)),
                        std::cos(0.2f * static_cast<Real>(l - m))};
    }
  }
  const auto field = Evaluate(e);
  const auto back = Expand(field, lMax);
  for (auto l : e.Degrees()) {
    for (auto m : e.Orders(l)) {
      EXPECT_NEAR(std::abs((back[l, m]) - (e[l, m])), 0, tolerance);
    }
  }

  // |f|^2 is a real scalar, 2 f - f is f, and an integer scales a field in
  // single precision as it does in double.
  auto modulus = abs2(field);
  static_assert(decltype(modulus)::UpperIndex == 0);
  SpinField<2, Grid> same = 2 * field - field;
  EXPECT_NEAR(std::abs((same[3, 5]) - (field[3, 5])), 0, tolerance);
  EXPECT_GE((modulus[3, 5]), 0.0f);
}

TEST(SinglePrecision, ATensorContractsAndRoundTrips) {
  auto grid = Grid(lMax, 2, FFTWpp::Estimate);
  using Strain = TensorField<2, Symmetric<2>, RealTensor, Grid>;
  auto strain = Strain(grid);
  auto j = Int{0};
  for (auto& x : strain.Data()) {
    x = Complex{std::sin(0.01f * static_cast<Real>(j)),
                std::cos(0.013f * static_cast<Real>(j))};
    j++;
  }
  for (auto& x : strain.RealData())
    x = std::sin(0.02f * static_cast<Real>(j++));

  // The trace of a real symmetric tensor is real.
  const auto& named = strain;
  auto trace = Trace(named);
  EXPECT_NEAR((Complex{trace[2, 3]}).imag(), 0, tolerance);

  // And its expansion evaluates back to it, component by component.
  auto bandLimited = Evaluate(Expand(named, lMax));
  auto again = Evaluate(Expand(std::as_const(bandLimited), lMax));
  for (std::size_t k = 0; k < again.Data().size(); k++) {
    EXPECT_NEAR(std::abs(again.Data()[k] - bandLimited.Data()[k]), 0,
                10 * tolerance);
  }
}

TEST(SinglePrecision, ALayeredFieldTransformsAndDifferentiates) {
  auto grid = Grid(lMax, 0, FFTWpp::Estimate);
  auto radii = std::vector<Real>{};
  for (auto i = 0; i < 9; i++)
    radii.push_back(0.5f + 0.0625f * static_cast<Real>(i));
  const auto radial = RadialGrid<Real>(radii);

  // f(r, theta, phi) = r^2 Y(theta, phi): its radial derivative is 2 r Y.
  auto e = LayeredSpinExpansion<0, Grid>(radial, grid, lMax);
  for (auto i : radial.RadiusIndices()) {
    const auto r = radial.Radius(i);
    e[i, 2, 1] = Complex{r * r, 0};
    e[i, 3, -2] = Complex{0, r * r};
  }
  const auto ddr = FiniteDifferenceDerivative<Real>(radial, 2);
  const auto derivative = ApplyRadially(e, ddr);
  for (auto i : radial.RadiusIndices()) {
    const auto r = radial.Radius(i);
    EXPECT_NEAR(std::abs((derivative[i, 2, 1]) - Complex{2 * r, 0}), 0, 1e-4f);
    EXPECT_NEAR(std::abs((derivative[i, 3, -2]) - Complex{0, 2 * r}), 0, 1e-4f);
  }

  const auto field = Evaluate(e);
  const auto back = Expand(field, lMax);
  for (auto i : radial.RadiusIndices()) {
    EXPECT_NEAR(std::abs((back[i, 2, 1]) - (e[i, 2, 1])), 0, tolerance);
  }
}

//--------------------------------------------------------------------------//
//                  A grid of non-negative upper indices only                //
//--------------------------------------------------------------------------//

TEST(NonNegativeUpperIndices, RoundTripThroughEveryKernel) {
  using Double = double;
  using ComplexD = std::complex<Double>;
  using HalfGrid = GaussLegendreGrid<Double, All, NonNegative>;
  constexpr Int degree = 12;
  constexpr Int nMax = 2;

  const auto check = [&](const HalfGrid& grid, const char* name) {
    for (auto n = Int{0}; n <= nMax; n++) {
      const auto size =
          static_cast<std::size_t>(grid.CoefficientSize(degree, n));
      auto given = std::vector<ComplexD>(size);
      for (std::size_t k = 0; k < size; k++) {
        given[k] = ComplexD{std::sin(0.37 * static_cast<Double>(k)),
                            std::cos(0.11 * static_cast<Double>(k))};
      }
      auto field = std::vector<ComplexD>(grid.FieldSize());
      auto back = std::vector<ComplexD>(size);
      grid.InverseTransformation(degree, n, given, field);
      grid.ForwardTransformation(degree, n, field, back);
      for (std::size_t k = 0; k < size; k++) {
        EXPECT_NEAR(std::abs(back[k] - given[k]), 0, 1e-12)
            << name << ", upper index " << n << ", coefficient " << k;
      }
    }
  };

  check(HalfGrid(degree, nMax, FFTWpp::Estimate), "stored");
  check(HalfGrid(degree, nMax, FFTWpp::Estimate, Chunking::Automatic(),
                 WignerValues::Generated()),
        "generated");
#ifdef GSHTRANS_HAVE_BLAS
  check(HalfGrid(degree, nMax, FFTWpp::Estimate, Chunking::Automatic(),
                 WignerValues::Stored(), TransformKernel::Matrix()),
        "matrix");
#endif

  // An upper index the grid does not cover is refused and not served wrongly.
  const auto grid = HalfGrid(degree, nMax, FFTWpp::Estimate);
  auto given = std::vector<ComplexD>(
      static_cast<std::size_t>(grid.CoefficientSize(degree, 1)));
  auto field = std::vector<ComplexD>(grid.FieldSize());
  EXPECT_THROW(grid.InverseTransformation(degree, -1, given, field),
               std::invalid_argument);
}
