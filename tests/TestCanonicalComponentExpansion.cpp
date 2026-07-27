#include <gtest/gtest.h>

#include <GSHTrans/All>

#include <complex>
#include <cstddef>
#include <ranges>
#include <type_traits>
#include <utility>

namespace {

using namespace GSHTrans;

using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;
using RealExpansion = RealCanonicalComponentExpansion<1, Grid>;
using ComplexExpansion = ComplexCanonicalComponentExpansion<1, Grid>;

void ExpectNear(Complex actual, Complex expected) {
  EXPECT_NEAR(actual.real(), expected.real(), 1.0e-12);
  EXPECT_NEAR(actual.imag(), expected.imag(), 1.0e-12);
}

template <typename Expansion>
void Fill(Expansion& expansion, Complex offset) {
  auto i = std::ptrdiff_t{0};
  for (auto [l, m] : expansion.Indices()) {
    expansion[l, m] =
        offset + Complex{static_cast<Real>(i), -static_cast<Real>(2 * i)};
    ++i;
  }
}

template <typename Expansion, typename Function>
void ExpectValues(const Expansion& expansion, Function expected) {
  for (auto [l, m] : expansion.Indices()) {
    ExpectNear(expansion[l, m], expected(l, m));
  }
}

template <typename Expansion>
void CheckAssignmentAndArithmetic() {
  auto destinationGrid = Grid(3, 1, FFTWpp::Estimate);
  auto sourceGrid = Grid(3, 1, FFTWpp::Estimate);
  auto destination = Expansion(destinationGrid);
  auto source = Expansion(sourceGrid);
  auto other = Expansion(sourceGrid);
  Fill(source, Complex{2.0, -1.0});
  Fill(other, Complex{-0.5, 3.0});

  destination = source;
  EXPECT_EQ(&destination.Grid(), &destinationGrid);
  ExpectValues(destination,
               [&](auto l, auto m) { return source[l, m]; });

  destination += other;
  ExpectValues(destination, [&](auto l, auto m) {
    return source[l, m] + other[l, m];
  });

  destination = source;
  destination -= other;
  ExpectValues(destination, [&](auto l, auto m) {
    return source[l, m] - other[l, m];
  });

  const auto scalar = Complex{1.5, -0.25};
  destination = source;
  destination *= scalar;
  ExpectValues(destination, [&](auto l, auto m) {
    return source[l, m] * scalar;
  });
  destination /= scalar;
  ExpectValues(destination,
               [&](auto l, auto m) { return source[l, m]; });

  destination += Expansion(other);
  ExpectValues(destination, [&](auto l, auto m) {
    return source[l, m] + other[l, m];
  });

  Fill(source, Complex{8.0, 2.0});
  destination = std::move(source);
  EXPECT_EQ(&destination.Grid(), &destinationGrid);
  ExpectValues(destination,
               [&](auto l, auto m) { return source[l, m]; });
}

}  // namespace

TEST(CanonicalComponentExpansion, TraitsAndOrderRanges) {
  static_assert(std::same_as<RealExpansion::Scalar, Complex>);
  static_assert(std::same_as<ComplexExpansion::Scalar, Complex>);
  static_assert(std::same_as<RealExpansion::MRange, NonNegative>);
  static_assert(std::same_as<ComplexExpansion::MRange, All>);
  static_assert(
      !std::is_base_of_v<GSHIndices<NonNegative>, RealExpansion>);
  static_assert(!std::is_base_of_v<GSHIndices<All>, ComplexExpansion>);
  static_assert(!std::is_default_constructible_v<RealExpansion>);
  static_assert(!std::is_default_constructible_v<ComplexExpansion>);
  static_assert(std::same_as<RealScalarExpansion<Grid>,
                             ScalarExpansion<Grid, RealValued>>);
  static_assert(std::same_as<ComplexScalarExpansion<Grid>,
                             ScalarExpansion<Grid, ComplexValued>>);

  auto grid = Grid(3, 1, FFTWpp::Estimate);
  auto realExpansion = RealExpansion(grid);
  auto complexExpansion = ComplexExpansion(grid);

  EXPECT_EQ(realExpansion.MinDegree(), 1);
  EXPECT_EQ(realExpansion.MaxDegree(), 3);
  EXPECT_EQ(realExpansion.MaxOrder(), 3);
  EXPECT_EQ(realExpansion.Size(), GSHIndices<NonNegative>(3, 3, 1).Size());
  EXPECT_EQ(complexExpansion.Size(), GSHIndices<All>(3, 3, 1).Size());

  auto realIndex = std::ptrdiff_t{0};
  for (auto [l, m] : realExpansion.Indices()) {
    EXPECT_GE(m, 0);
    EXPECT_LE(m, l);
    EXPECT_EQ(realExpansion.Index(l, m), realIndex);
    ++realIndex;
  }
  EXPECT_EQ(realIndex, realExpansion.Size());

  auto complexIndex = std::ptrdiff_t{0};
  auto sawNegativeOrder = false;
  for (auto [l, m] : complexExpansion.Indices()) {
    sawNegativeOrder = sawNegativeOrder || m < 0;
    EXPECT_GE(m, -l);
    EXPECT_LE(m, l);
    EXPECT_EQ(complexExpansion.Index(l, m), complexIndex);
    ++complexIndex;
  }
  EXPECT_TRUE(sawNegativeOrder);
  EXPECT_EQ(complexIndex, complexExpansion.Size());

  auto orders = realExpansion.Orders();
  auto order = orders.begin();
  for (auto [l, m] : realExpansion.Indices()) {
    static_cast<void>(l);
    ASSERT_NE(order, orders.end());
    EXPECT_EQ(*order, m);
    ++order;
  }
}

TEST(CanonicalComponentExpansion, ComplexStorageAndDataAccess) {
  using MutableData = decltype(std::declval<RealExpansion&>().Data());
  using ConstData = decltype(std::declval<const RealExpansion&>().Data());
  static_assert(std::ranges::output_range<MutableData, Complex>);
  static_assert(!std::ranges::output_range<ConstData, Complex>);

  auto grid = Grid(3, 1, FFTWpp::Estimate);
  auto realExpansion = RealExpansion(grid);
  auto complexExpansion = ComplexExpansion(grid);
  Fill(realExpansion, Complex{1.0, 4.0});
  Fill(complexExpansion, Complex{-2.0, 5.0});

  auto sawImaginaryValue = false;
  for (auto [l, m] : realExpansion.Indices()) {
    const auto i = realExpansion.Index(l, m);
    ExpectNear(realExpansion.Data()[i], realExpansion[l, m]);
    sawImaginaryValue =
        sawImaginaryValue || realExpansion[l, m].imag() != 0.0;
  }
  EXPECT_TRUE(sawImaginaryValue);
  for (auto [l, m] : complexExpansion.Indices()) {
    const auto i = complexExpansion.Index(l, m);
    ExpectNear(complexExpansion.Data()[i], complexExpansion[l, m]);
  }

  const auto& constExpansion = realExpansion;
  EXPECT_EQ(constExpansion.Data().size(), realExpansion.Size());
  ExpectNear(constExpansion.Data()[0], realExpansion.Data()[0]);
}

TEST(CanonicalComponentExpansion, RealAssignmentAndArithmetic) {
  CheckAssignmentAndArithmetic<RealExpansion>();
}

TEST(CanonicalComponentExpansion, ComplexAssignmentAndArithmetic) {
  CheckAssignmentAndArithmetic<ComplexExpansion>();
}
