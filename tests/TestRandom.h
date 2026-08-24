#ifndef GSH_TRANS_TEST_RANDOM_GUARD_H
#define GSH_TRANS_TEST_RANDOM_GUARD_H

#include <GSHTrans/All>
#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <ranges>
#include <string>

// Seeded random data for the tests.
//
// Everything here is explicit about its seed, because drawing a fresh
// std::random_device per call means a test that fails cannot be rerun on the
// data that failed it. Tests draw one seed from
// TestSeed(), thread it through, and report it on failure; setting
// GSHTRANS_TEST_SEED to that value reproduces the run exactly, and setting it
// to "random" goes back to unseeded exploration.

namespace GSHTransTest {

using Int = std::ptrdiff_t;
using Seed = std::uint64_t;
using Generator = std::mt19937_64;

inline Seed TestSeed() {
  constexpr Seed defaultSeed = 20260819u;
  const auto* environment = std::getenv("GSHTRANS_TEST_SEED");
  if (environment == nullptr) return defaultSeed;
  const auto text = std::string(environment);
  if (text == "random") return std::random_device{}();
  return static_cast<Seed>(std::stoull(text));
}

inline auto MakeGenerator(Seed seed) { return Generator(seed); }

inline Int RandomDegree(Generator& gen, Int lMin, Int lMax) {
  return std::uniform_int_distribution<Int>(lMin, lMax)(gen);
}

template <GSHTrans::IndexRange NRange>
Int RandomUpperIndex(Generator& gen, Int nMax) {
  if constexpr (std::same_as<NRange, GSHTrans::All>) {
    return std::uniform_int_distribution<Int>(-nMax, nMax)(gen);
  } else {
    return std::uniform_int_distribution<Int>(0, nMax)(gen);
  }
}

// Fill `range` with the coefficients of a random complex-valued field of
// degree lMax and upper index n.
template <typename Grid, std::ranges::range Range>
requires requires() {
  requires GSHTrans::ComplexFloatingPoint<std::ranges::range_value_t<Range>>;
  requires std::ranges::output_range<Range, std::ranges::range_value_t<Range>>;
}
void RandomComplexCoefficient(const Grid& grid, Int lMax, Int n, Range& range,
                              Generator& gen) {
  using Complex = std::ranges::range_value_t<Range>;
  using Real = GSHTrans::RemoveComplex<Complex>;
  assert(range.size() == grid.CoefficientSize(lMax, n));

  auto dist = std::normal_distribution<Real>();
  std::ranges::generate(
      range, [&gen, &dist]() { return Complex{dist(gen), dist(gen)}; });
}

// Fill `range` with the reduced (m >= 0) coefficients of a random real-valued
// field of degree lMax. There is no upper index: real-valued fields exist only
// at n = 0.
template <typename Grid, std::ranges::range Range>
requires requires() {
  requires GSHTrans::ComplexFloatingPoint<std::ranges::range_value_t<Range>>;
  requires std::ranges::output_range<Range, std::ranges::range_value_t<Range>>;
}
void RandomRealCoefficient(const Grid& grid, Int lMax, Range& range,
                           Generator& gen) {
  using Complex = std::ranges::range_value_t<Range>;
  using Real = GSHTrans::RemoveComplex<Complex>;
  assert(range.size() == grid.RealCoefficientSize(lMax));

  auto dist = std::normal_distribution<Real>();
  std::ranges::generate(
      range, [&gen, &dist]() { return Complex{dist(gen), dist(gen)}; });

  // The m = 0 coefficients of a real field are real.
  auto indices = GSHTrans::GSHIndices<GSHTrans::NonNegative>(lMax, lMax, 0);
  for (auto l : indices.Degrees()) {
    range[indices.Index(l, 0)].imag(0);
  }
}

}  // namespace GSHTransTest

#endif  // GSH_TRANS_TEST_RANDOM_GUARD_H
