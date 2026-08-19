
#include <GSHTrans/All>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>

using namespace GSHTrans;

using Int = std::ptrdiff_t;

// One generator for the whole example, seeded explicitly so that a run can be
// repeated. The library no longer supplies random-coefficient generators;
// they were test scaffolding on GridBase and were removed (core-plan.md F12).
std::mt19937_64 gen(20260819u);

Int RandomDegree(Int lMin, Int lMax) {
  return std::uniform_int_distribution<Int>(lMin, lMax)(gen);
}

template <IndexRange NRange>
Int RandomUpperIndex(Int nMax) {
  if constexpr (std::same_as<NRange, All>) {
    return std::uniform_int_distribution<Int>(-nMax, nMax)(gen);
  } else {
    return std::uniform_int_distribution<Int>(0, nMax)(gen);
  }
}

// Random coefficients of a complex-valued field. The (lMax, lMax) coefficient
// is unresolvable while nPhi = 2*lMax, so it is left at zero; that workaround
// goes with core-plan.md step D.
template <typename Grid, typename Range>
void RandomComplexCoefficient(const Grid& grid, Int lMax, Int n, Range& range) {
  using Complex = std::ranges::range_value_t<Range>;
  using Real = RemoveComplex<Complex>;
  auto dist = std::normal_distribution<Real>();
  std::ranges::generate(range,
                        [&dist]() { return Complex{dist(gen), dist(gen)}; });
  if (lMax == grid.MaxDegree()) {
    range[GSHIndices<All>(lMax, lMax, n).Index(lMax, lMax)] = 0;
  }
}

int main() {
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Complex;
  using MRange = All;
  using NRange = All;
  using Grid = GaussLegendreGrid<Real, MRange, NRange>;

  {
    auto lMaxGrid = RandomDegree(4, 256);
    auto lMax = RandomDegree(4, lMaxGrid);
    auto nMax = std::min(lMax, Int(4));

    auto grid = Grid(lMaxGrid, nMax);
    auto n = RandomUpperIndex<NRange>(nMax);

    Int size;
    if constexpr (ComplexFloatingPoint<Scalar>) {
      size = grid.CoefficientSize(lMax, n);
    } else {
      size = grid.RealCoefficientSize(lMax);
    }
    auto flm = FFTWpp::vector<Complex>(size);
    RandomComplexCoefficient(grid, lMax, n, flm);

    auto f = FFTWpp::vector<Scalar>(grid.FieldSize());
    auto glm = FFTWpp::vector<Complex>(size);

    grid.InverseTransformation(lMax, n, flm, f);

    grid.ForwardTransformation(lMax, n, f, glm);

    auto error = std::ranges::max(std::ranges::views::zip_transform(
        [](auto x, auto y) { return std::abs(x - y); }, flm, glm));

    ;

    std::cout << lMaxGrid << " " << lMax << " " << n << " " << error
              << std::endl;
  }

  FFTWpp::CleanUp();
}