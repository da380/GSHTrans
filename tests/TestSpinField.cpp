#include <gtest/gtest.h>

#include <GSHTrans/All>
#include <complex>
#include <cstddef>
#include <cmath>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

// Test family 1: the compile-time algebra.
//
// Almost everything here is a static_assert, so the runtime bodies are empty
// and the value is that the file compiles. It grows through the phase-1 steps:
// this instalment covers the node concept and the traits underneath it, before
// any node exists to use them.

using namespace GSHTrans;
using Int = std::ptrdiff_t;

//--------------------------------------------------------------------------//
//                 A stub node, declared but never defined                   //
//--------------------------------------------------------------------------//

// Concept satisfaction needs declarations, not definitions, so nothing here is
// implemented. The point is to exercise SpinWeighted independently of the
// terminal, which arrives in step 2.
template <Int N, typename ValueTag>
struct StubNode {
  static constexpr Int UpperIndex = N;
  using Value = ValueTag;
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = ScalarFor<Real, Value>;
  using GridType = GaussLegendreGrid<Real, All, All>;

  const GridType& Grid() const;
  Scalar operator[](Int iTheta, Int iPhi) const;
  template <typename S>
  void EvaluateInto(std::span<S> target) const;
};

struct StubTerminal : StubNode<1, ComplexValued> {};
struct StubExpression : StubNode<1, ComplexValued> {};

namespace GSHTrans {
template <>
struct IsTerminalTrait<StubTerminal> : std::true_type {};
}  // namespace GSHTrans

// Nodes that fail the concept, each in exactly one way.
struct MissingUpperIndex {
  using Value = ComplexValued;
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Complex;
  using GridType = GaussLegendreGrid<Real, All, All>;
  const GridType& Grid() const;
  Scalar operator[](Int, Int) const;
  template <typename S>
  void EvaluateInto(std::span<S>) const;
};

struct ReturnsReference : StubNode<1, ComplexValued> {
  const Scalar& operator[](Int, Int) const;
};

struct MissingEvaluateInto {
  static constexpr Int UpperIndex = 0;
  using Value = ComplexValued;
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Complex;
  using GridType = GaussLegendreGrid<Real, All, All>;
  const GridType& Grid() const;
  Scalar operator[](Int, Int) const;
};

struct MismatchedScalar {
  static constexpr Int UpperIndex = 0;
  using Value = ComplexValued;
  using Real = double;
  using Complex = std::complex<Real>;
  using Scalar = Real;  // must be Complex when Value is ComplexValued
  using GridType = GaussLegendreGrid<Real, All, All>;
  const GridType& Grid() const;
  Scalar operator[](Int, Int) const;
  template <typename S>
  void EvaluateInto(std::span<S>) const;
};

namespace {

//--------------------------------------------------------------------------//
//                              The node concept                             //
//--------------------------------------------------------------------------//

static_assert(SpinWeighted<StubNode<0, RealValued>>);
static_assert(SpinWeighted<StubNode<0, ComplexValued>>);
static_assert(SpinWeighted<StubNode<1, ComplexValued>>);
static_assert(SpinWeighted<StubNode<-2, ComplexValued>>);

// The reality constraint: no node may claim to be real-valued away from zero.
static_assert(!SpinWeighted<StubNode<1, RealValued>>);
static_assert(!SpinWeighted<StubNode<-1, RealValued>>);
static_assert(!SpinWeighted<StubNode<2, RealValued>>);

// Each of the interface requirements bites.
static_assert(!SpinWeighted<MissingUpperIndex>);
static_assert(!SpinWeighted<ReturnsReference>);
static_assert(!SpinWeighted<MissingEvaluateInto>);
static_assert(!SpinWeighted<MismatchedScalar>);
static_assert(!SpinWeighted<double>);
static_assert(!SpinWeighted<GaussLegendreGrid<double, All, All>>);

// Scalar follows from Value.
static_assert(std::same_as<ScalarFor<double, RealValued>, double>);
static_assert(std::same_as<ScalarFor<double, ComplexValued>,
                           std::complex<double>>);
static_assert(std::same_as<ScalarFor<long double, ComplexValued>,
                           std::complex<long double>>);

//--------------------------------------------------------------------------//
//                       Terminals, views and operands                       //
//--------------------------------------------------------------------------//

static_assert(IsTerminal<StubTerminal>);
static_assert(!IsTerminal<StubExpression>);

// The trait sees through references and cv-qualification, so a caller's value
// category cannot change what kind of thing an operand is.
static_assert(IsTerminal<StubTerminal&>);
static_assert(IsTerminal<const StubTerminal&>);
static_assert(IsTerminal<StubTerminal&&>);

// An lvalue terminal is referred to; everything else is owned.
static_assert(
    std::same_as<OperandStorage<StubTerminal&>, const StubTerminal&>);
static_assert(
    std::same_as<OperandStorage<const StubTerminal&>, const StubTerminal&>);

// The row a plain IsTerminal test gets wrong: an rvalue terminal must be moved
// in, not referred to, or `auto e = MakeField(...) + v;` dangles.
static_assert(std::same_as<OperandStorage<StubTerminal>, StubTerminal>);

// Expression nodes are held by value whichever way they arrive, so that a
// named expression may go out of scope before an expression built from it.
static_assert(std::same_as<OperandStorage<StubExpression&>, StubExpression>);
static_assert(std::same_as<OperandStorage<StubExpression>, StubExpression>);
static_assert(
    std::same_as<OperandStorage<const StubExpression&>, StubExpression>);

//--------------------------------------------------------------------------//
//                             Value propagation                             //
//--------------------------------------------------------------------------//

static_assert(std::same_as<CombinedValue<RealValued, RealValued>, RealValued>);
static_assert(
    std::same_as<CombinedValue<RealValued, ComplexValued>, ComplexValued>);
static_assert(
    std::same_as<CombinedValue<ComplexValued, RealValued>, ComplexValued>);
static_assert(
    std::same_as<CombinedValue<ComplexValued, ComplexValued>, ComplexValued>);

// A real scalar preserves; a complex one promotes.
static_assert(std::same_as<ValueAfterScalar<double, RealValued>, RealValued>);
static_assert(
    std::same_as<ValueAfterScalar<std::complex<double>, RealValued>,
                 ComplexValued>);
static_assert(
    std::same_as<ValueAfterScalar<double, ComplexValued>, ComplexValued>);

//--------------------------------------------------------------------------//
//                             The index algebra                             //
//--------------------------------------------------------------------------//

namespace R = IndexRules;

// Equal: addition and subtraction.
static_assert(R::Equal::Admissible<2, 2>);
static_assert(R::Equal::Admissible<-1, -1>);
static_assert(!R::Equal::Admissible<1, -1>);
static_assert(!R::Equal::Admissible<0, 2>);
static_assert(R::Equal::Apply<2, 2> == 2);
static_assert(R::Equal::Apply<-2, -2> == -2);

// Sum: multiplication. Always admissible, and the indices add.
static_assert(R::Sum::Admissible<2, -2>);
static_assert(R::Sum::Apply<1, 1> == 2);
static_assert(R::Sum::Apply<2, -2> == 0);
static_assert(R::Sum::Apply<-1, 0> == -1);

// FirstOnly: division, by a scalar field only.
static_assert(R::FirstOnly::Admissible<2, 0>);
static_assert(!R::FirstOnly::Admissible<2, 1>);
static_assert(!R::FirstOnly::Admissible<0, -1>);
static_assert(R::FirstOnly::Apply<2, 0> == 2);

// Same, Negate, Zero: the unary rules.
static_assert(R::Same::Apply<2> == 2);
static_assert(R::Same::Apply<-2> == -2);
static_assert(R::Negate::Apply<2> == -2);
static_assert(R::Negate::Apply<-2> == 2);
static_assert(R::Negate::Apply<0> == 0);
static_assert(R::Zero::Apply<2> == 0);
static_assert(R::Zero::Apply<-2> == 0);

// Conjugating twice returns the original index.
static_assert(R::Negate::Apply<R::Negate::Apply<3>> == 3);

//--------------------------------------------------------------------------//
//        The closure lemma: the reality constraint costs nothing            //
//--------------------------------------------------------------------------//

// Every rule preserves "real-valued implies upper index zero", so the concept
// checks it once and no node has to re-derive it. Each line below is one row
// of the closure table in plan section 3.4: given operands that satisfy the
// constraint, the result does too.
template <typename Value, Int N>
inline constexpr bool Lawful =
    std::same_as<Value, ComplexValued> or N == 0;

// Equal, on two real operands: both must be at zero, so the result is.
static_assert(Lawful<CombinedValue<RealValued, RealValued>,
                     R::Equal::Apply<0, 0>>);
// Sum, on two real operands: 0 + 0 = 0.
static_assert(
    Lawful<CombinedValue<RealValued, RealValued>, R::Sum::Apply<0, 0>>);
// Sum, where one operand is complex: the result is complex, so unconstrained.
static_assert(
    Lawful<CombinedValue<RealValued, ComplexValued>, R::Sum::Apply<0, 2>>);
// FirstOnly, on two real operands.
static_assert(
    Lawful<CombinedValue<RealValued, RealValued>, R::FirstOnly::Apply<0, 0>>);
// Same, with a real scalar: value and index both unchanged.
static_assert(Lawful<ValueAfterScalar<double, RealValued>, R::Same::Apply<0>>);
// Same, with a complex scalar: promoted, so the premise is discharged.
static_assert(Lawful<ValueAfterScalar<std::complex<double>, ComplexValued>,
                     R::Same::Apply<2>>);
// Negate: real implies zero implies -0 = 0.
static_assert(Lawful<RealValued, R::Negate::Apply<0>>);
// Zero: lands at zero unconditionally, so a real result is always lawful.
static_assert(Lawful<RealValued, R::Zero::Apply<2>>);
static_assert(Lawful<RealValued, R::Zero::Apply<-2>>);

//--------------------------------------------------------------------------//
//                    The terminal satisfies the concept                     //
//--------------------------------------------------------------------------//

using Grid = GaussLegendreGrid<double, All, All>;
using ScalarGrid = GaussLegendreGrid<double, NonNegative, All>;
using Complex = std::complex<double>;

static_assert(SpinWeighted<SpinField<0, Grid, RealValued>>);
static_assert(SpinWeighted<SpinField<0, Grid, ComplexValued>>);
static_assert(SpinWeighted<SpinField<2, Grid>>);
static_assert(SpinWeighted<SpinField<-2, Grid>>);

// ComplexValued is the default, since it is the only lawful choice away from
// zero.
static_assert(std::same_as<SpinField<2, Grid>::Value, ComplexValued>);
static_assert(std::same_as<SpinField<2, Grid>::Scalar, std::complex<double>>);
static_assert(std::same_as<SpinField<0, Grid, RealValued>::Scalar, double>);

static_assert(IsTerminal<SpinField<2, Grid>>);
static_assert(!std::is_default_constructible_v<SpinField<2, Grid>>);
static_assert(std::copy_constructible<SpinField<2, Grid>>);

// A real-valued field at nonzero upper index, and a negative upper index on a
// grid that stores only non-negative ones, are rejected by static_assert
// inside the class. Those cannot be written as static_assert(!...) here,
// because naming the type is not enough to fire them and instantiating it is a
// hard error by design. The equivalent facts are tested through the concept
// above, on stub nodes.
static_assert(SpinWeighted<SpinField<0, ScalarGrid, RealValued>>);

// EvaluateInto widens a real field into a complex destination and refuses the
// reverse, which would be a truncation. Written as a concept because a bare
// requires-expression naming a member template of a non-dependent type is a
// hard error on GCC rather than a false constraint.
template <typename Field, typename S>
concept EvaluatesInto = requires(const Field& field, std::span<S> target) {
  field.EvaluateInto(target);
};

static_assert(EvaluatesInto<SpinField<0, Grid, RealValued>, double>);
static_assert(EvaluatesInto<SpinField<0, Grid, RealValued>, Complex>);
static_assert(EvaluatesInto<SpinField<2, Grid>, Complex>);
static_assert(!EvaluatesInto<SpinField<2, Grid>, double>);
static_assert(!EvaluatesInto<SpinField<0, Grid, ComplexValued>, double>);

}  // namespace

TEST(SpinField, CompileTimeAlgebraIsPinned) { SUCCEED(); }

namespace {

constexpr Int lMax = 6;
constexpr Int nMax = 2;

auto TestGrid() { return Grid(lMax, nMax, FFTWpp::Estimate); }

}  // namespace

TEST(SpinField, StoresSamplesInTheCanonicalOrder) {
  auto grid = TestGrid();
  auto u = SpinField<2, Grid>(grid);

  ASSERT_EQ(u.Size(), grid.FieldSize());
  const auto nPhi = static_cast<Int>(grid.NumberOfLongitudes());

  // A freshly built field is zero.
  for (auto value : u) EXPECT_EQ(value, Complex{});

  // Mutable access writes where the flat layout says it should: phi fastest.
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      u[iTheta, iPhi] = Complex{static_cast<double>(iTheta),
                                static_cast<double>(iPhi)};
    }
  }
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      const auto flat = iTheta * nPhi + iPhi;
      EXPECT_EQ(u.Data()[flat], (Complex{static_cast<double>(iTheta),
                                         static_cast<double>(iPhi)}));
      EXPECT_EQ((u[iTheta, iPhi]), u.Data()[flat]);
    }
  }
}

TEST(SpinField, SamplesAFunctionOfPosition) {
  auto grid = TestGrid();
  auto f = [](auto theta, auto phi) {
    return Complex{std::cos(theta), std::sin(phi)};
  };
  auto u = SpinField<1, Grid>(grid, f);

  auto i = Int{0};
  for (auto [theta, phi] : grid.Points()) {
    EXPECT_EQ(u.Data()[i], f(theta, phi));
    ++i;
  }
  EXPECT_EQ(i, grid.FieldSize());
}

TEST(SpinField, EvaluateIntoAgreesWithTheElementLoop) {
  auto grid = TestGrid();
  auto u = SpinField<2, Grid>(grid, [](auto theta, auto phi) {
    return Complex{theta * phi, theta - phi};
  });

  auto target = std::vector<Complex>(u.Size());
  u.EvaluateInto(std::span(target));

  auto i = Int{0};
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      EXPECT_EQ(target[i], (u[iTheta, iPhi]));
      ++i;
    }
  }
}

TEST(SpinField, EvaluateIntoWidensRealToComplexButNotBack) {
  auto grid = TestGrid();
  auto u = SpinField<0, Grid, RealValued>(
      grid, [](auto theta, auto phi) { return theta + phi; });

  auto target = std::vector<Complex>(u.Size());
  u.EvaluateInto(std::span(target));
  for (auto i = Int{0}; i < u.Size(); ++i) {
    EXPECT_EQ(target[i].real(), u.Data()[i]);
    EXPECT_EQ(target[i].imag(), 0.0);
  }
}

TEST(SpinField, RejectsAMismatchedEvaluationTarget) {
  auto grid = TestGrid();
  auto u = SpinField<2, Grid>(grid);
  auto shortTarget = std::vector<Complex>(u.Size() - 1);
  auto longTarget = std::vector<Complex>(u.Size() + 1);

  EXPECT_THROW(u.EvaluateInto(std::span(shortTarget)), std::invalid_argument);
  EXPECT_THROW(u.EvaluateInto(std::span(longTarget)), std::invalid_argument);
}

TEST(SpinField, RejectsAnUpperIndexTheGridDoesNotCarry) {
  auto narrow = Grid(lMax, 1, FFTWpp::Estimate);
  EXPECT_THROW((SpinField<2, Grid>(narrow)), std::invalid_argument);
  EXPECT_THROW((SpinField<-2, Grid>(narrow)), std::invalid_argument);
  EXPECT_NO_THROW((SpinField<1, Grid>(narrow)));
  EXPECT_NO_THROW((SpinField<0, Grid>(narrow)));

  // The message says what the grid does carry.
  try {
    auto rejected = SpinField<2, Grid>(narrow);
    FAIL() << "expected a throw";
  } catch (const std::invalid_argument& error) {
    EXPECT_NE(std::string(error.what()).find("upper index 2"),
              std::string::npos)
        << error.what();
  }
}

//--------------------------------------------------------------------------//
//                  Family 1: the operator matrix over N                     //
//--------------------------------------------------------------------------//

namespace {

// Whether an operator is available at all. Written as concepts so that a
// negative case is a false constraint rather than a hard error.
template <typename L, typename R>
concept Addable = requires(L l, R r) { l + r; };
template <typename L, typename R>
concept Subtractable = requires(L l, R r) { l - r; };
template <typename L, typename R>
concept Multipliable = requires(L l, R r) { l * r; };
template <typename L, typename R>
concept Divisible = requires(L l, R r) { l / r; };
template <typename A>
concept HasRealPart = requires(A a) { real(a); };
template <typename A>
concept HasImagPart = requires(A a) { imag(a); };
template <typename A>
concept HasAbs = requires(A a) { abs(a); };
template <typename A>
concept Negatable = requires(A a) { -a; };
template <typename A>
concept HasConj = requires(A a) { conj(a); };

template <Int N>
using F = SpinField<N, Grid>;
using F0R = SpinField<0, Grid, RealValued>;
using OtherGrid = GaussLegendreGrid<long double, All, All>;

template <typename L, typename R>
using Add = decltype(std::declval<L>() + std::declval<R>());
template <typename L, typename R>
using Mul = decltype(std::declval<L>() * std::declval<R>());
template <typename L, typename R>
using Div = decltype(std::declval<L>() / std::declval<R>());

// Addition and subtraction need equal upper indices, over the whole range the
// library exposes.
static_assert(Addable<F<-2>&, F<-2>&>);
static_assert(Addable<F<-1>&, F<-1>&>);
static_assert(Addable<F<0>&, F<0>&>);
static_assert(Addable<F<1>&, F<1>&>);
static_assert(Addable<F<2>&, F<2>&>);
static_assert(!Addable<F<1>&, F<-1>&>);
static_assert(!Addable<F<0>&, F<2>&>);
static_assert(!Addable<F<2>&, F<1>&>);
static_assert(!Subtractable<F<1>&, F<-1>&>);
static_assert(Add<F<2>&, F<2>&>::UpperIndex == 2);
static_assert(Add<F<-2>&, F<-2>&>::UpperIndex == -2);

// Multiplication is always admissible and the indices add.
static_assert(Multipliable<F<2>&, F<-2>&>);
static_assert(Mul<F<1>&, F<1>&>::UpperIndex == 2);
static_assert(Mul<F<2>&, F<-2>&>::UpperIndex == 0);
static_assert(Mul<F<-1>&, F<0>&>::UpperIndex == -1);
static_assert(Mul<F<2>&, F<2>&>::UpperIndex == 4);

// Division needs the divisor at zero, and keeps the dividend's index.
static_assert(Divisible<F<2>&, F<0>&>);
static_assert(!Divisible<F<2>&, F<1>&>);
static_assert(!Divisible<F<0>&, F<-1>&>);
static_assert(Div<F<2>&, F<0>&>::UpperIndex == 2);

// Unary minus and conjugation exist everywhere; conj reverses the index and
// doing it twice returns the original.
static_assert(Negatable<F<2>&> && HasConj<F<2>&>);
static_assert(decltype(-std::declval<F<2>&>())::UpperIndex == 2);
static_assert(decltype(conj(std::declval<F<2>&>()))::UpperIndex == -2);
static_assert(decltype(conj(std::declval<F<-2>&>()))::UpperIndex == 2);
static_assert(
    decltype(conj(conj(std::declval<F<2>&>())))::UpperIndex == 2);

// abs and abs2 exist at every index and land at zero, real-valued.
static_assert(HasAbs<F<2>&>);
static_assert(decltype(abs(std::declval<F<2>&>()))::UpperIndex == 0);
static_assert(std::same_as<decltype(abs(std::declval<F<2>&>()))::Value,
                           RealValued>);
static_assert(std::same_as<decltype(abs2(std::declval<F<2>&>()))::Value,
                           RealValued>);
static_assert(std::same_as<decltype(abs2(std::declval<F<2>&>()))::Scalar,
                           double>);

// real and imag exist only at zero.
static_assert(HasRealPart<F<0>&> && HasImagPart<F<0>&>);
static_assert(!HasRealPart<F<1>&>);
static_assert(!HasImagPart<F<1>&>);
static_assert(!HasRealPart<F<-2>&>);
static_assert(std::same_as<decltype(real(std::declval<F<0>&>()))::Value,
                           RealValued>);

// Value propagation through the algebra.
static_assert(std::same_as<Add<F0R&, F0R&>::Value, RealValued>);
static_assert(std::same_as<Add<F0R&, F<0>&>::Value, ComplexValued>);
static_assert(std::same_as<Mul<F0R&, F0R&>::Value, RealValued>);
static_assert(std::same_as<Mul<F0R&, F<2>&>::Value, ComplexValued>);
static_assert(std::same_as<Div<F0R&, F0R&>::Value, RealValued>);

// conj on a real-valued field at zero is the identity, and stays real.
static_assert(std::same_as<decltype(conj(std::declval<F0R&>()))::Value,
                           RealValued>);
static_assert(decltype(conj(std::declval<F0R&>()))::UpperIndex == 0);

// Scalar multiplication: real preserves, complex promotes, index untouched.
static_assert(
    std::same_as<decltype(std::declval<F0R&>() * 2.0)::Value, RealValued>);
static_assert(std::same_as<
              decltype(std::declval<F0R&>() * Complex{0, 1})::Value,
              ComplexValued>);
static_assert(std::same_as<decltype(2.0 * std::declval<F<2>&>())::Value,
                           ComplexValued>);
static_assert(decltype(2.0 * std::declval<F<2>&>())::UpperIndex == 2);
static_assert(decltype(std::declval<F<2>&>() / 2.0)::UpperIndex == 2);

// A scalar over a field needs the field at zero.
static_assert(Divisible<double, F<0>&>);
static_assert(!Divisible<double, F<2>&>);

// One precision per tree: no mixed-precision promotion, and no mixing grids of
// different precision.
static_assert(!Addable<F<0>&, SpinField<0, OtherGrid>&>);
static_assert(!Multipliable<F<0>&, SpinField<0, OtherGrid>&>);
static_assert(!Multipliable<F<2>&, float>);
static_assert(!Multipliable<F<2>&, std::complex<float>>);

// Every lawful expression is itself a node, so expressions compose.
static_assert(SpinWeighted<Add<F<2>&, F<2>&>>);
static_assert(SpinWeighted<Mul<F<2>&, F<-2>&>>);
static_assert(SpinWeighted<decltype(abs(std::declval<F<2>&>()))>);
static_assert(Addable<Add<F<2>&, F<2>&>, F<2>&>);
static_assert(Multipliable<decltype(conj(std::declval<F<2>&>())), F<2>&>);

// Expression nodes are copyable and auto never slices, because operands are
// stored as const T& or T and never as a base reference.
static_assert(std::copy_constructible<Add<F<2>&, F<2>&>>);
static_assert(std::copy_constructible<Mul<F<2>&, F<-2>&>>);
static_assert(std::copy_constructible<decltype(abs2(std::declval<F<2>&>()))>);

// The showcase for the index algebra: conj(f) carries -N, so the product with
// g at +N lands at zero and can be integrated. It type-checks at every N.
static_assert(Multipliable<decltype(conj(std::declval<F<2>&>())), F<2>&>);
static_assert(Mul<decltype(conj(std::declval<F<2>&>())), F<2>&>::UpperIndex ==
              0);
static_assert(Mul<decltype(conj(std::declval<F<-1>&>())), F<-1>&>::UpperIndex ==
              0);

}  // namespace

//--------------------------------------------------------------------------//
//                 Family 3 (part): the algebra computes                     //
//--------------------------------------------------------------------------//

namespace {

auto MakeField(const Grid& grid, double a) {
  return SpinField<2, Grid>(grid, [a](auto theta, auto phi) {
    return Complex{a * std::cos(theta) + 0.3, a * std::sin(phi) - 0.2};
  });
}

auto MakeScalarField(const Grid& grid, double a) {
  return SpinField<0, Grid>(grid, [a](auto theta, auto phi) {
    return Complex{a + std::cos(theta) * std::cos(phi), 0.5 * a + theta};
  });
}

// Evaluate a node into a plain vector. Named Evaluated rather than
// Materialise, which is now the library's own and returns a SpinField.
template <typename NodeType>
auto Evaluated(const NodeType& node) {
  auto out = std::vector<typename NodeType::Scalar>(node.Grid().FieldSize());
  node.EvaluateInto(std::span(out));
  return out;
}

constexpr double tolerance = 1.0e-14;

void ExpectClose(Complex a, Complex b) {
  EXPECT_NEAR(a.real(), b.real(), tolerance);
  EXPECT_NEAR(a.imag(), b.imag(), tolerance);
}

}  // namespace

TEST(SpinField, PointwiseOperationsAgreeWithScalarArithmetic) {
  auto grid = TestGrid();
  auto u = MakeField(grid, 1.5);
  auto v = MakeField(grid, -0.75);
  auto s = MakeScalarField(grid, 2.0);

  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      const auto a = u[iTheta, iPhi];
      const auto b = v[iTheta, iPhi];
      const auto c = s[iTheta, iPhi];

      ExpectClose((u + v)[iTheta, iPhi], a + b);
      ExpectClose((u - v)[iTheta, iPhi], a - b);
      ExpectClose((u * v)[iTheta, iPhi], a * b);
      ExpectClose((u / s)[iTheta, iPhi], a / c);
      ExpectClose((-u)[iTheta, iPhi], -a);
      ExpectClose(conj(u)[iTheta, iPhi], std::conj(a));
      EXPECT_NEAR((abs(u)[iTheta, iPhi]), std::abs(a), tolerance);
      EXPECT_NEAR((abs2(u)[iTheta, iPhi]), std::norm(a), tolerance);
      EXPECT_NEAR((real(s)[iTheta, iPhi]), c.real(), tolerance);
      EXPECT_NEAR((imag(s)[iTheta, iPhi]), c.imag(), tolerance);
      ExpectClose((u * 2.0)[iTheta, iPhi], a * 2.0);
      ExpectClose((2.0 * u)[iTheta, iPhi], 2.0 * a);
      ExpectClose((u / 2.0)[iTheta, iPhi], a / 2.0);
      ExpectClose((Complex{0, 1} * u)[iTheta, iPhi], Complex{0, 1} * a);
      ExpectClose((2.0 / s)[iTheta, iPhi], 2.0 / c);
    }
  }
}

TEST(SpinField, EvaluateIntoAgreesWithTheElementLoopOnExpressions) {
  auto grid = TestGrid();
  auto u = MakeField(grid, 1.5);
  auto v = MakeField(grid, -0.75);

  const auto expression = conj(u) * v + u * conj(v);
  const auto evaluated = Evaluated(expression);

  auto i = Int{0};
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      ExpectClose(evaluated[i], expression[iTheta, iPhi]);
      ++i;
    }
  }
  EXPECT_EQ(i, grid.FieldSize());
}

TEST(SpinField, BinaryNodesRejectOperandsOnDifferentGrids) {
  auto grid = TestGrid();
  auto other = TestGrid();
  ASSERT_NE(grid.Identity(), other.Identity());

  auto u = MakeField(grid, 1.0);
  auto v = MakeField(other, 1.0);

  EXPECT_THROW(auto node = u + v, std::invalid_argument);
  EXPECT_THROW(auto node = u * v, std::invalid_argument);
  EXPECT_THROW(auto node = u - v, std::invalid_argument);

  // Equal parameters are not enough; it is handle identity that decides.
  EXPECT_EQ(u.Grid().MaxDegree(), v.Grid().MaxDegree());
  EXPECT_EQ(u.Size(), v.Size());

  try {
    auto node = u + v;
    FAIL() << "expected a throw";
  } catch (const std::invalid_argument& error) {
    EXPECT_NE(std::string(error.what()).find("different"), std::string::npos)
        << error.what();
  }
}

//--------------------------------------------------------------------------//
//                          Family 2: lifetime                               //
//--------------------------------------------------------------------------//
//
// These are the cases the layer being replaced got wrong: it stored operands
// as references to a CRTP base, so a nested expression bound to auto dangled.
// They are worth little in a plain build and everything under the sanitisers,
// which is where the suite also runs.

TEST(SpinField, NestedExpressionsBoundToAutoStayValid) {
  auto grid = TestGrid();
  auto u = MakeField(grid, 1.5);
  auto v = MakeField(grid, -0.75);
  auto w = MakeScalarField(grid, 2.0);

  auto inner = u + v;
  auto middle = inner * w;
  auto outer = conj(middle) + conj(inner) * w;

  const auto evaluated = Evaluated(outer);
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      const auto a = u[iTheta, iPhi];
      const auto b = v[iTheta, iPhi];
      const auto c = w[iTheta, iPhi];
      const auto expected = std::conj((a + b) * c) + std::conj(a + b) * c;
      ExpectClose(evaluated[iTheta * grid.NumberOfLongitudes() + iPhi],
                  expected);
    }
  }
}

TEST(SpinField, ExpressionsOwnRvalueTerminals) {
  auto grid = TestGrid();
  auto v = MakeField(grid, -0.75);

  // The terminal on the left is a temporary. A plain IsTerminal test would
  // store a reference to it and this would dangle.
  auto expression = MakeField(grid, 1.5) + v;

  const auto evaluated = Evaluated(expression);
  const auto reference = MakeField(grid, 1.5);
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(evaluated[i], reference.Data()[i] + v.Data()[i]);
  }
}

TEST(SpinField, NamedExpressionsMayDieBeforeWhatIsBuiltFromThem) {
  auto grid = TestGrid();
  auto u = MakeField(grid, 1.5);
  auto w = MakeScalarField(grid, 2.0);

  // `inner` is destroyed at the end of the lambda; `outer` holds a copy of it,
  // not a reference.
  auto outer = [&] {
    auto inner = u + u;
    return inner * w;
  }();

  const auto evaluated = Evaluated(outer);
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(evaluated[i],
                (u.Data()[i] + u.Data()[i]) * w.Data()[i]);
  }
}

TEST(SpinField, ExpressionsSurviveBeingReturnedAndStored) {
  auto grid = TestGrid();
  auto u = MakeField(grid, 1.5);
  auto v = MakeField(grid, -0.75);

  // Returned from a function: the operands are lvalue terminals held by
  // reference, and they outlive the expression, which is the documented
  // contract.
  auto make = [](const auto& a, const auto& b) { return a * conj(b); };
  auto returned = make(u, v);

  // Stored in a container, which needs copy construction.
  auto nodes = std::vector<decltype(returned)>{};
  nodes.push_back(returned);
  nodes.push_back(make(u, v));
  nodes.push_back(nodes.front());

  for (const auto& node : nodes) {
    const auto evaluated = Evaluated(node);
    for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
      ExpectClose(evaluated[i], u.Data()[i] * std::conj(v.Data()[i]));
    }
  }
}

//--------------------------------------------------------------------------//
//              Family 3: materialisation, assignment, aliasing              //
//--------------------------------------------------------------------------//

namespace {

template <typename F, typename E>
concept PlusAssignable = requires(F f, E e) { f += e; };
template <typename F, typename E>
concept TimesAssignable = requires(F f, E e) { f *= e; };
template <typename F, typename E>
concept DivideAssignable = requires(F f, E e) { f /= e; };
template <typename F, typename E>
concept AssignableFrom = requires(F f, E e) { f = e; };

using F0 = SpinField<0, Grid>;
using F2 = SpinField<2, Grid>;

// Assignment needs the same upper index; nothing is coerced.
static_assert(AssignableFrom<F2&, Add<F2&, F2&>>);
static_assert(!AssignableFrom<F2&, Add<F0&, F0&>>);
static_assert(!AssignableFrom<F0&, Mul<F2&, F2&>>);

// A real-valued expression may be assigned into a complex field; the reverse
// would be a truncation and does not compile.
static_assert(AssignableFrom<F0&, decltype(abs(std::declval<F2&>()))>);
static_assert(!AssignableFrom<F0R&, Add<F0&, F0&>>);
static_assert(AssignableFrom<F0R&, decltype(abs(std::declval<F2&>()))>);

// Compound assignment inherits the index rule of its binary form.
static_assert(PlusAssignable<F2&, F2&>);
static_assert(!PlusAssignable<F2&, F0&>);
static_assert(TimesAssignable<F2&, F0&>);
static_assert(!TimesAssignable<F2&, F2&>);   // would land at N = 4
static_assert(DivideAssignable<F2&, F0&>);
static_assert(!DivideAssignable<F2&, F2&>);
static_assert(TimesAssignable<F2&, double>);
static_assert(TimesAssignable<F2&, Complex>);
static_assert(DivideAssignable<F2&, double>);

// A complex scalar cannot multiply a real field in place: the field's value
// kind cannot change under assignment.
static_assert(TimesAssignable<F0R&, double>);
static_assert(!TimesAssignable<F0R&, Complex>);

// Materialise returns a field of the expression's own index and value kind.
static_assert(std::same_as<decltype(Materialise(std::declval<Add<F2&, F2&>>())),
                           SpinField<2, Grid, ComplexValued>>);
static_assert(
    std::same_as<decltype(Materialise(std::declval<
                          decltype(abs(std::declval<F2&>()))>())),
                 SpinField<0, Grid, RealValued>>);

}  // namespace

TEST(SpinField, ConstructionFromAnExpressionTakesTheExpressionsGrid) {
  auto grid = TestGrid();
  auto u = MakeField(grid, 1.5);
  auto v = MakeField(grid, -0.75);

  SpinField<2, Grid> w = u + v;
  EXPECT_EQ(w.Grid().Identity(), grid.Identity());
  EXPECT_EQ(w.Size(), grid.FieldSize());
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(w.Data()[i], u.Data()[i] + v.Data()[i]);
  }

  // And through Materialise, which is the same thing with the type deduced.
  auto m = Materialise(u * conj(v));
  static_assert(std::same_as<decltype(m), SpinField<0, Grid, ComplexValued>>);
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(m.Data()[i], u.Data()[i] * std::conj(v.Data()[i]));
  }
}

TEST(SpinField, RealExpressionsWidenIntoComplexFields) {
  auto grid = TestGrid();
  auto u = MakeField(grid, 1.5);

  SpinField<0, Grid> wide = abs2(u);
  auto narrow = Materialise(abs2(u));
  static_assert(std::same_as<decltype(narrow),
                             SpinField<0, Grid, RealValued>>);

  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    EXPECT_NEAR(wide.Data()[i].real(), std::norm(u.Data()[i]), tolerance);
    EXPECT_NEAR(wide.Data()[i].imag(), 0.0, tolerance);
    EXPECT_NEAR(narrow.Data()[i], std::norm(u.Data()[i]), tolerance);
  }
}

TEST(SpinField, AssignmentDoesNotRebindTheGrid) {
  auto grid = TestGrid();
  auto other = TestGrid();
  auto u = MakeField(grid, 1.5);
  auto v = MakeField(grid, -0.75);
  auto elsewhere = MakeField(other, 1.0);

  u = v + v;
  EXPECT_EQ(u.Grid().Identity(), grid.Identity());
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(u.Data()[i], v.Data()[i] + v.Data()[i]);
  }

  EXPECT_THROW(u = elsewhere + elsewhere, std::invalid_argument);
  // and the destination is untouched by the attempt
  EXPECT_EQ(u.Grid().Identity(), grid.Identity());
}

TEST(SpinField, CompoundAssignmentMatchesItsBinaryForm) {
  auto grid = TestGrid();
  auto u = MakeField(grid, 1.5);
  auto v = MakeField(grid, -0.75);
  auto s = MakeScalarField(grid, 2.0);
  const auto u0 = u;

  u += v;
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(u.Data()[i], u0.Data()[i] + v.Data()[i]);
  }

  u = u0;
  u -= v;
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(u.Data()[i], u0.Data()[i] - v.Data()[i]);
  }

  u = u0;
  u *= s;
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(u.Data()[i], u0.Data()[i] * s.Data()[i]);
  }

  u = u0;
  u /= s;
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(u.Data()[i], u0.Data()[i] / s.Data()[i]);
  }

  u = u0;
  u *= Complex{0.0, 2.0};
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(u.Data()[i], u0.Data()[i] * Complex{0.0, 2.0});
  }

  u = u0;
  u /= 4.0;
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(u.Data()[i], u0.Data()[i] / 4.0);
  }
}

// The aliasing theorem: every node is pointwise and index-preserving, so an
// assignment whose right-hand side mentions the destination reads element
// (iTheta, iPhi) only when writing that same element. No temporary is needed.
// This is the regression that fails if a re-indexing node ever enters the
// layer.
TEST(SpinField, InPlaceAssignmentIsSafeWhenTheDestinationAppears) {
  auto grid = TestGrid();
  auto u = SpinField<0, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::cos(theta) + 1.5, std::sin(phi) - 0.25};
  });
  auto v = MakeScalarField(grid, 0.5);
  const auto u0 = u;

  u = conj(u) * v + u;

  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    const auto expected =
        std::conj(u0.Data()[i]) * v.Data()[i] + u0.Data()[i];
    ExpectClose(u.Data()[i], expected);
  }

  // The same for compound assignment, and for a destination appearing twice.
  auto w = u0;
  w += w * v;
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    ExpectClose(w.Data()[i], u0.Data()[i] + u0.Data()[i] * v.Data()[i]);
  }

  auto x = u0;
  x = x * x - x;
  for (auto i = Int{0}; i < grid.FieldSize(); ++i) {
    const auto a = u0.Data()[i];
    ExpectClose(x.Data()[i], a * a - a);
  }
}

// A deep mixed tree, lazily evaluated and materialised, must agree exactly at
// every upper index the library exposes.
template <Int N>
void CheckLazyAgreesWithMaterialised() {
  auto grid = Grid(lMax, nMax, FFTWpp::Estimate);
  auto a = SpinField<N, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::cos(theta) + 1.25, std::sin(phi) - 0.3};
  });
  auto b = SpinField<N, Grid>(grid, [](auto theta, auto phi) {
    return Complex{std::sin(theta) * 0.5 + 1.0, std::cos(phi) + 0.4};
  });
  auto c = SpinField<0, Grid>(grid, [](auto theta, auto phi) {
    return Complex{2.0 + std::cos(theta * phi), 0.75};
  });

  // Upper index N throughout: conj(a) carries -N, so conj(a) * b lands at
  // zero, and multiplying that scalar back into a returns to N.
  const auto lazy = (a + b) * c - a * (conj(a) * b) / c + b * 2.0;

  auto stepOne = Materialise(a + b);
  auto stepTwo = Materialise(stepOne * c);
  auto stepThree = Materialise(conj(a) * b);
  auto stepFour = Materialise(a * stepThree);
  auto stepFive = Materialise(stepFour / c);
  auto stepSix = Materialise(b * 2.0);

  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      const auto flat = iTheta * grid.NumberOfLongitudes() + iPhi;
      const auto expected =
          stepTwo.Data()[flat] - stepFive.Data()[flat] + stepSix.Data()[flat];
      ExpectClose(lazy[iTheta, iPhi], expected);
    }
  }
}

TEST(SpinField, LazyAndMaterialisedAgreeAtEveryUpperIndex) {
  CheckLazyAgreesWithMaterialised<-2>();
  CheckLazyAgreesWithMaterialised<-1>();
  CheckLazyAgreesWithMaterialised<0>();
  CheckLazyAgreesWithMaterialised<1>();
  CheckLazyAgreesWithMaterialised<2>();
}

TEST(SpinField, CopiesDataButSharesTheGrid) {
  auto grid = TestGrid();
  auto u = SpinField<2, Grid>(grid);
  u[0, 0] = Complex{1.0, 2.0};

  auto v = u;
  EXPECT_EQ(u.Grid().Identity(), v.Grid().Identity());
  EXPECT_EQ((v[0, 0]), (Complex{1.0, 2.0}));

  // The data is a copy, not a share.
  v[0, 0] = Complex{3.0, 4.0};
  EXPECT_EQ((u[0, 0]), (Complex{1.0, 2.0}));

  // And the field keeps the grid alive on its own.
  auto escaped = [&grid] {
    auto local = Grid(lMax, nMax, FFTWpp::Estimate);
    return SpinField<2, Grid>(local);
  }();
  EXPECT_EQ(escaped.Size(), grid.FieldSize());
  EXPECT_NE(escaped.Grid().Identity(), grid.Identity());
}
