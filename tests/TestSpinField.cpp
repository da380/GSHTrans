#include <gtest/gtest.h>

#include <GSHTrans/All>
#include <complex>
#include <cstddef>
#include <span>
#include <type_traits>

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

}  // namespace

TEST(SpinField, CompileTimeAlgebraIsPinned) { SUCCEED(); }
