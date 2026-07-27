#include <gtest/gtest.h>

#include <GSHTrans/All>

#include <complex>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace {

using namespace GSHTrans;

using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;
using RealField = RealCanonicalComponentField<0, Grid>;
using ComplexField = ComplexCanonicalComponentField<0, Grid>;

void ExpectNear(Real actual, Real expected) {
  EXPECT_NEAR(actual, expected, 1.0e-12);
}

void ExpectNear(Complex actual, Complex expected) {
  EXPECT_NEAR(actual.real(), expected.real(), 1.0e-12);
  EXPECT_NEAR(actual.imag(), expected.imag(), 1.0e-12);
}

template <typename Field, typename Function>
void Fill(Field& field, Function value) {
  for (auto [iTheta, iPhi] : field.PointIndices()) {
    field[iTheta, iPhi] = value(iTheta, iPhi);
  }
}

template <typename Field, typename Function>
void ExpectValues(const Field& field, Function expected) {
  for (auto [iTheta, iPhi] : field.PointIndices()) {
    ExpectNear(field[iTheta, iPhi], expected(iTheta, iPhi));
  }
}

template <typename Field>
void CheckCompoundArithmetic(Field& u, const Field& v,
                             typename Field::Scalar scalar) {
  Field result(u.Grid());

  result = u;
  result += v;
  ExpectValues(result, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] + v[iTheta, iPhi];
  });

  result = u;
  result -= v;
  ExpectValues(result, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] - v[iTheta, iPhi];
  });

  result = u;
  result *= v;
  ExpectValues(result, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] * v[iTheta, iPhi];
  });

  result = u;
  result /= v;
  ExpectValues(result, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] / v[iTheta, iPhi];
  });

  result = u;
  result += scalar;
  ExpectValues(result, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] + scalar;
  });

  result = u;
  result -= scalar;
  ExpectValues(result, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] - scalar;
  });

  result = u;
  result *= scalar;
  ExpectValues(result, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] * scalar;
  });

  result = u;
  result /= scalar;
  ExpectValues(result, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] / scalar;
  });
}

class MoveOnlyUnary {
 public:
  explicit MoveOnlyUnary(Real offset)
      : _offset{std::make_unique<Real>(offset)} {}

  MoveOnlyUnary(const MoveOnlyUnary&) = delete;
  MoveOnlyUnary(MoveOnlyUnary&&) = default;

  Real operator()(Real value) const { return value + *_offset; }

 private:
  std::unique_ptr<Real> _offset;
};

class MoveOnlyScalarCallable {
 public:
  explicit MoveOnlyScalarCallable(Real factor)
      : _factor{std::make_unique<Real>(factor)} {}

  MoveOnlyScalarCallable(const MoveOnlyScalarCallable&) = delete;
  MoveOnlyScalarCallable(MoveOnlyScalarCallable&&) = default;

  Real operator()(Real value, Real scalar) const {
    return *_factor * value + scalar;
  }

 private:
  std::unique_ptr<Real> _factor;
};

class StatefulUnary {
 public:
  explicit StatefulUnary(Real offset) : _offset{offset} {}

  void SetOffset(Real offset) { _offset = offset; }
  Real operator()(Real value) const { return value + _offset; }

 private:
  Real _offset;
};

class StatefulScalarCallable {
 public:
  explicit StatefulScalarCallable(Real factor) : _factor{factor} {}

  void SetFactor(Real factor) { _factor = factor; }
  Real operator()(Real value, Real scalar) const {
    return _factor * value + scalar;
  }

 private:
  Real _factor;
};

class MutableUnary {
 public:
  Real operator()(Real value) { return value; }
};

class MutableScalarCallable {
 public:
  Real operator()(Real value, Real scalar) { return value + scalar; }
};

template <typename Callable>
concept FormsUnaryExpression = requires(RealField& field, Callable&& callable) {
  CanonicalComponentFieldUnary(field,
                               std::forward<Callable>(callable));
};

template <typename Callable>
concept FormsScalarExpression =
    requires(RealField& field, Callable&& callable) {
      CanonicalComponentFieldUnaryWithScalar(
          field, std::forward<Callable>(callable), Real{});
    };

static_assert(FormsUnaryExpression<MoveOnlyUnary>);
static_assert(FormsScalarExpression<MoveOnlyScalarCallable>);
static_assert(!FormsUnaryExpression<MutableUnary>);
static_assert(!FormsScalarExpression<MutableScalarCallable>);

}  // namespace

TEST(CanonicalComponentField, RealArithmeticAndMaterialization) {
  auto grid = Grid(2, 0, FFTWpp::Estimate);
  auto u = RealField(grid);
  auto v = RealField(grid);
  Fill(u, [](auto iTheta, auto iPhi) {
    return 2.0 + static_cast<Real>(iTheta + 2 * iPhi);
  });
  Fill(v, [](auto iTheta, auto iPhi) {
    return 5.0 + static_cast<Real>(2 * iTheta + iPhi);
  });

  auto sum = RealField(u + v);
  auto difference = RealField(u - v);
  auto product = RealField(u * v);
  auto quotient = RealField(u / v);
  auto negative = RealField(-u);
  auto shifted = RealField(u + 3.0);
  auto leftShifted = RealField(3.0 + u);
  auto reduced = RealField(u - 1.0);
  auto scaled = RealField(u * 2.5);
  auto leftScaled = RealField(2.5 * u);
  auto divided = RealField(u / 2.0);
  auto nested = RealField((u + v) - v);

  ExpectValues(sum, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] + v[iTheta, iPhi];
  });
  ExpectValues(difference, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] - v[iTheta, iPhi];
  });
  ExpectValues(product, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] * v[iTheta, iPhi];
  });
  ExpectValues(quotient, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] / v[iTheta, iPhi];
  });
  ExpectValues(negative,
               [&](auto iTheta, auto iPhi) { return -u[iTheta, iPhi]; });
  ExpectValues(shifted, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] + 3.0;
  });
  ExpectValues(leftShifted, [&](auto iTheta, auto iPhi) {
    return 3.0 + u[iTheta, iPhi];
  });
  ExpectValues(reduced, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] - 1.0;
  });
  ExpectValues(scaled, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] * 2.5;
  });
  ExpectValues(leftScaled, [&](auto iTheta, auto iPhi) {
    return 2.5 * u[iTheta, iPhi];
  });
  ExpectValues(divided, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] / 2.0;
  });
  ExpectValues(nested, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi];
  });

  CheckCompoundArithmetic(u, v, 2.0);
}

TEST(CanonicalComponentField, ComplexArithmeticAndUnaryExpressions) {
  auto grid = Grid(2, 0, FFTWpp::Estimate);
  auto u = ComplexField(grid);
  auto v = ComplexField(grid);
  Fill(u, [](auto iTheta, auto iPhi) {
    return Complex{2.0 + static_cast<Real>(iTheta),
                   1.0 + static_cast<Real>(iPhi)};
  });
  Fill(v, [](auto iTheta, auto iPhi) {
    return Complex{4.0 + static_cast<Real>(iPhi),
                   2.0 + static_cast<Real>(iTheta)};
  });

  auto sum = ComplexField(u + v);
  auto difference = ComplexField(u - v);
  auto product = ComplexField(u * v);
  auto quotient = ComplexField(u / v);
  auto conjugate = ComplexField(conj(u));
  auto realPart = RealField(real(u));
  auto imaginaryPart = RealField(imag(u));
  const auto scalar = Complex{1.5, -0.5};
  auto shifted = ComplexField(u - scalar);
  auto scaled = ComplexField(u / scalar);
  auto leftShifted = ComplexField(scalar + u);
  auto leftScaled = ComplexField(scalar * u);

  ExpectValues(sum, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] + v[iTheta, iPhi];
  });
  ExpectValues(difference, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] - v[iTheta, iPhi];
  });
  ExpectValues(product, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] * v[iTheta, iPhi];
  });
  ExpectValues(quotient, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] / v[iTheta, iPhi];
  });
  ExpectValues(conjugate, [&](auto iTheta, auto iPhi) {
    return std::conj(u[iTheta, iPhi]);
  });
  ExpectValues(realPart, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi].real();
  });
  ExpectValues(imaginaryPart, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi].imag();
  });
  ExpectValues(shifted, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] - scalar;
  });
  ExpectValues(scaled, [&](auto iTheta, auto iPhi) {
    return u[iTheta, iPhi] / scalar;
  });
  ExpectValues(leftShifted, [&](auto iTheta, auto iPhi) {
    return scalar + u[iTheta, iPhi];
  });
  ExpectValues(leftScaled, [&](auto iTheta, auto iPhi) {
    return scalar * u[iTheta, iPhi];
  });

  CheckCompoundArithmetic(u, v, scalar);
}

TEST(CanonicalComponentField, AssignmentPreservesDestinationGrid) {
  static_assert(!std::is_default_constructible_v<RealField>);
  static_assert(!std::is_default_constructible_v<ComplexField>);

  auto destinationGrid = Grid(2, 0, FFTWpp::Estimate);
  auto sourceGrid = Grid(2, 0, FFTWpp::Estimate);
  auto destination = RealField(destinationGrid);
  auto source = RealField(sourceGrid);
  Fill(source, [](auto iTheta, auto iPhi) {
    return 7.0 + static_cast<Real>(3 * iTheta + iPhi);
  });

  destination = source;
  EXPECT_EQ(&destination.Grid(), &destinationGrid);
  ExpectValues(destination, [&](auto iTheta, auto iPhi) {
    return source[iTheta, iPhi];
  });

  Fill(source, [](auto iTheta, auto iPhi) {
    return 11.0 + static_cast<Real>(iTheta + 4 * iPhi);
  });
  const auto expectedAfterMove = RealField(source);
  destination = std::move(source);
  EXPECT_EQ(&destination.Grid(), &destinationGrid);
  ExpectValues(destination, [&](auto iTheta, auto iPhi) {
    return expectedAfterMove[iTheta, iPhi];
  });

  destination = destination + destination;
  EXPECT_EQ(&destination.Grid(), &destinationGrid);
  ExpectValues(destination, [&](auto iTheta, auto iPhi) {
    return 2.0 * expectedAfterMove[iTheta, iPhi];
  });
}

TEST(CanonicalComponentField, AssignmentRejectsMismatchedGridSizes) {
  auto smallGrid = Grid(1, 0, FFTWpp::Estimate);
  auto largeGrid = Grid(2, 0, FFTWpp::Estimate);
  auto small = RealField(smallGrid);
  auto large = RealField(largeGrid);

  EXPECT_THROW(small = large, std::invalid_argument);
  EXPECT_THROW(small = std::move(large), std::invalid_argument);
  EXPECT_THROW(small = large + 1.0, std::invalid_argument);
  EXPECT_THROW(large = small + 1.0, std::invalid_argument);
}

TEST(CanonicalComponentField, CallableExpressionsOwnForwardedState) {
  auto grid = Grid(2, 0, FFTWpp::Estimate);
  auto field = RealField(grid);
  Fill(field, [](auto iTheta, auto iPhi) {
    return 1.0 + static_cast<Real>(iTheta + iPhi);
  });

  auto unaryExpression =
      CanonicalComponentFieldUnary(field, MoveOnlyUnary{4.0});
  auto scalarExpression = CanonicalComponentFieldUnaryWithScalar(
      field, MoveOnlyScalarCallable{3.0}, 2.0);

  auto unaryResult = RealField(unaryExpression);
  auto scalarResult = RealField(scalarExpression);
  ExpectValues(unaryResult, [&](auto iTheta, auto iPhi) {
    return field[iTheta, iPhi] + 4.0;
  });
  ExpectValues(scalarResult, [&](auto iTheta, auto iPhi) {
    return 3.0 * field[iTheta, iPhi] + 2.0;
  });
}

TEST(CanonicalComponentField, LvalueCallablesAreCopiedIntoExpressions) {
  auto grid = Grid(2, 0, FFTWpp::Estimate);
  auto field = RealField(grid);
  Fill(field, [](auto iTheta, auto iPhi) {
    return 1.0 + static_cast<Real>(iTheta + iPhi);
  });

  auto unaryCallable = StatefulUnary{4.0};
  auto scalarCallable = StatefulScalarCallable{3.0};
  auto unaryExpression =
      CanonicalComponentFieldUnary(field, unaryCallable);
  auto scalarExpression = CanonicalComponentFieldUnaryWithScalar(
      field, scalarCallable, 2.0);

  unaryCallable.SetOffset(40.0);
  scalarCallable.SetFactor(30.0);

  auto unaryResult = RealField(unaryExpression);
  auto scalarResult = RealField(scalarExpression);
  ExpectValues(unaryResult, [&](auto iTheta, auto iPhi) {
    return field[iTheta, iPhi] + 4.0;
  });
  ExpectValues(scalarResult, [&](auto iTheta, auto iPhi) {
    return 3.0 * field[iTheta, iPhi] + 2.0;
  });
}
