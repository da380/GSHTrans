#ifndef GSH_TRANS_SPIN_WEIGHTED_GUARD_H
#define GSH_TRANS_SPIN_WEIGHTED_GUARD_H

#include <complex>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "../Concepts.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                                 The node                                  //
//--------------------------------------------------------------------------//

// A spin-weighted field of definite upper index N on the sphere. Everything in
// this layer -- the owning terminal, non-owning views, and every expression
// node -- models SpinWeighted, and nothing dispatches through anything else:
// there are no virtual functions and no required base class.
//
// The unit is deliberately a single component, not a tensor. For rank >= 2 a
// collection labelled only by an upper index does not determine a tensor (see
// section 2 of the theory note, docs/canonical-components.tex), so tensors
// are built from these at the
// tensor layer rather than being what this layer is about.

// The scalar a node evaluates to: Real when the node is real-valued, Complex
// otherwise.
template <RealFloatingPoint Real, RealOrComplexValued Value>
using ScalarFor =
    std::conditional_t<std::same_as<Value, RealValued>, Real,
                       std::complex<Real>>;

template <typename T>
concept SpinWeighted = requires {
  // The spin weight, known at compile time.
  requires std::same_as<std::remove_cv_t<decltype(T::UpperIndex)>,
                        std::ptrdiff_t>;

  // The scalar kind, and the types that follow from it.
  typename T::Value;
  requires RealOrComplexValued<typename T::Value>;
  typename T::Real;
  requires RealFloatingPoint<typename T::Real>;
  typename T::Complex;
  requires std::same_as<typename T::Complex, std::complex<typename T::Real>>;
  typename T::Scalar;
  requires std::same_as<typename T::Scalar,
                        ScalarFor<typename T::Real, typename T::Value>>;

  typename T::GridType;

  // Theory note section 7 item 5, promoted into the concept: real-valuedness
  // is not preserved by the frame rotation e_{+-} -> e^{-+ i psi} e_{+-}, so
  // it is not a property any component of any tensor can have at N != 0. A
  // library able to represent one can represent something that does not
  // exist. The constraint is closed under every node in this layer, so it is
  // checked here once and never re-derived; what it forbids is a real-valued
  // *terminal* or *view* at nonzero upper index.
  requires std::same_as<typename T::Value, ComplexValued> or
               T::UpperIndex == 0;
} and requires(const T& node, std::ptrdiff_t iTheta, std::ptrdiff_t iPhi,
               std::span<typename T::Scalar> target) {
  // The grid is a value-semantic handle, so this is cheap to copy and two
  // nodes are on the same grid when Grid().Identity() agrees.
  { node.Grid() } -> std::convertible_to<const typename T::GridType&>;

  // By value, always, including on terminals: uniform value return is what
  // makes terminals, views and expressions interchangeable. A mutable
  // reference accessor exists on terminals and mutable views, outside the
  // concept.
  { node[iTheta, iPhi] } -> std::same_as<typename T::Scalar>;

  // Writes the field in the canonical layout: (iTheta, iPhi) with phi
  // fastest, flat index iTheta * nPhi + iPhi. This declaration is the
  // definition of record for that order.
  { node.EvaluateInto(target) } -> std::same_as<void>;
};

// Convenience for the common case of writing "the node type behind this
// possibly-reference operand type".
template <typename T>
using Node = std::remove_cvref_t<T>;

//--------------------------------------------------------------------------//
//                    What may be evaluated into a field                     //
//--------------------------------------------------------------------------//

// Whether an expression can be written into a field of upper index N,
// precision Real and scalar Scalar: the same upper index, the same precision,
// and a scalar that converts. The last condition is what permits a
// real-valued expression into a complex field and makes the reverse a compile
// error rather than a truncation.
//
// **A concept, and it has to be one.** This was a `static constexpr bool`
// member of SpinField, written as a chain of `and`s beginning with
// SpinWeighted. That is ill-formed for any operand that is not spin-weighted,
// because `and` short-circuits *evaluation* and not *well-formedness*: the
// initialiser still names Node<Expr>::UpperIndex, which does not exist. GCC
// accepted it and clang did not, so a grid passed where a field's constructor
// was being considered failed to compile there and nowhere else.
//
// Concept conjunction is the construct that actually short-circuits: an
// atomic constraint is only checked once the ones before it are satisfied, so
// the members below are never named for a type that is not spin-weighted.
template <typename Expr, std::ptrdiff_t N, typename Real, typename Scalar>
concept EvaluatesInto =
    SpinWeighted<Node<Expr>> && (Node<Expr>::UpperIndex == N) &&
    std::same_as<typename Node<Expr>::Real, Real> &&
    std::convertible_to<typename Node<Expr>::Scalar, Scalar>;

//--------------------------------------------------------------------------//
//                          What this layer needs of a grid                  //
//--------------------------------------------------------------------------//

// Deliberately much less than a grid offers. Stating it as a concept
// documents the coupling, gives a readable error when a grid is missing
// something, and is what the layered and tensor layers check against too.
template <typename G>
concept AngularGrid = requires(const G& grid) {
  typename G::Real;
  requires RealFloatingPoint<typename G::Real>;
  typename G::NRange;

  // Value-semantic: a terminal holds one by value and expressions copy it.
  requires std::copy_constructible<G>;

  // Identity, not structure: two grids are the same grid when they share an
  // implementation.
  { grid.Identity() } -> std::equality_comparable;

  { grid.FieldSize() } -> std::integral;
  { grid.MaxUpperIndex() } -> std::integral;
  grid.UpperIndices();

  // The point set, and the index ranges an evaluation loop runs over.
  grid.Points();
  grid.CoLatitudeIndices();
  grid.LongitudeIndices();

  // The two axes separately, which is more than an evaluation loop needs and
  // is required anyway.
  //
  // The expression layer only ever walks Points(), so this used to be left
  // out and Interpolate carried its own refinement asking for it -- a
  // rectilinear interpolant takes one abscissa range per axis, and the polar
  // padding has to build each of them. Grids on this library's terms are
  // separable, since SphericalGrid holds the two axes and forms Points() from
  // them, so the refinement was describing every grid there is. Asking here
  // instead means one concept rather than two, and Points() becomes a
  // convenience the grid already provides rather than a separate demand.
  //
  // Written against range_value_t rather than *begin(...): an axis accessor
  // returns a view by value, and views::all of a prvalue container is an
  // owning_view, which is not a borrowed range -- so ranges::begin on the
  // returned prvalue is ill-formed for such a grid even though it plainly has
  // the axis. That cost a debugging session when Interpolate's local concept
  // was first written.
  requires std::convertible_to<
      std::ranges::range_value_t<decltype(grid.CoLatitudes())>,
      typename G::Real>;
  requires std::convertible_to<
      std::ranges::range_value_t<decltype(grid.Longitudes())>,
      typename G::Real>;

  // Quadrature weights, kept as two factors rather than one product: the
  // longitude weights are uniform, so the sphere integral factorises.
  grid.CoLatitudeWeights();
  grid.LongitudeWeights();
};

//--------------------------------------------------------------------------//
//                         The generic evaluation loop                       //
//--------------------------------------------------------------------------//

// The default EvaluateInto: read every point through operator[] and write it
// in the canonical order, phi fastest. Terminals override this with a
// contiguous copy; expressions and views use it as it stands.
//
// Evaluation is const and stateless, so several threads may evaluate the same
// node concurrently into disjoint targets. That is a documented guarantee, not
// an accident: the layered layer parallelises over slices on it.
template <typename NodeType, typename S>
void EvaluateNodeInto(const NodeType& node, std::span<S> target) {
  const auto& grid = node.Grid();
  const auto size = static_cast<std::size_t>(grid.FieldSize());
  if (target.size() != size) {
    throw std::invalid_argument(
        "Evaluation target has size " + std::to_string(target.size()) +
        ", but this field has " + std::to_string(size) + " points");
  }
  auto iter = target.begin();
  for (auto iTheta : grid.CoLatitudeIndices()) {
    for (auto iPhi : grid.LongitudeIndices()) {
      *iter++ = static_cast<S>(node[iTheta, iPhi]);
    }
  }
}

//--------------------------------------------------------------------------//
//                          Terminals versus expressions                     //
//--------------------------------------------------------------------------//

// Specialised by the owning field and the view types. Deliberately a trait
// rather than a property inferred from the interface: whether a node owns
// storage is not something its accessors reveal, and guessing it wrong decides
// whether an expression stores a reference or a copy.
template <typename T>
struct IsTerminalTrait : std::false_type {};

template <typename T>
inline constexpr bool IsTerminal = IsTerminalTrait<std::remove_cvref_t<T>>::value;

// How an operand is stored inside a node, given the value category it was
// passed with. T is the deduced type of a forwarding reference, so it carries
// the category: an lvalue arrives as U&, an rvalue as U.
//
//   terminal or view, lvalue    ->  const U&    (cheap; caller owns it)
//   terminal or view, rvalue    ->  U           (moved in; nothing else owns it)
//   expression node,  either    ->  U           (small; copying is cheap)
//
// The rvalue-terminal row is the one a plain IsTerminal test gets wrong:
// `auto e = MakeField(...) + v;` would bind a reference to a temporary that
// dies at the end of the full expression. Expression nodes are always held by
// value so that `auto f = e * w;` stays valid when the named expression `e`
// goes out of scope first.
//
// Residual hazard, and it is not removable in C++: an lvalue terminal
// destroyed while an expression referring to it is still alive. Eigen has the
// same one. Owning terminals through a shared_ptr would fix it and would
// change the cost model of every field, so it is documented instead.
template <typename T>
using OperandStorage =
    std::conditional_t<IsTerminal<T> and std::is_lvalue_reference_v<T>,
                       const std::remove_cvref_t<T>&, std::remove_cvref_t<T>>;

//--------------------------------------------------------------------------//
//                              Value propagation                            //
//--------------------------------------------------------------------------//

// A binary node is real-valued exactly when both its operands are.
template <RealOrComplexValued LValue, RealOrComplexValued RValue>
using CombinedValue =
    std::conditional_t<std::same_as<LValue, RealValued> and
                           std::same_as<RValue, RealValued>,
                       RealValued, ComplexValued>;

// The value kind a scalar multiplier imposes: a complex scalar promotes a real
// field, a real one leaves it alone.
template <RealOrComplexFloatingPoint S, RealOrComplexValued Value>
using ValueAfterScalar =
    std::conditional_t<ComplexFloatingPoint<S>, ComplexValued, Value>;

//--------------------------------------------------------------------------//
//                              The index algebra                            //
//--------------------------------------------------------------------------//

// Each rule carries the compile-time computation of the result's upper index
// and the admissibility condition on its operands. The conditions are applied
// as requires-clauses on the free operators, so an unlawful combination is an
// overload-resolution failure at the call site rather than an error inside a
// node -- which also lets a negative test be written as
// static_assert(!requires { u + v; }).
namespace IndexRules {

using Int = std::ptrdiff_t;

// Binary rules.

// Addition and subtraction: the upper index must agree and is carried through.
struct Equal {
  template <Int NL, Int NR>
  static constexpr bool Admissible = NL == NR;
  template <Int NL, Int NR>
  static constexpr Int Apply = NL;
};

// Multiplication: upper indices add. Theory note eq:N.
struct Sum {
  template <Int NL, Int NR>
  static constexpr bool Admissible = true;
  template <Int NL, Int NR>
  static constexpr Int Apply = NL + NR;
};

// Division: only by a scalar field, so the divisor must be at upper index
// zero and the dividend's index survives.
struct FirstOnly {
  template <Int NL, Int NR>
  static constexpr bool Admissible = NR == 0;
  template <Int NL, Int NR>
  static constexpr Int Apply = NL;
};

// Unary rules.

// Negation, and multiplication or division by a scalar.
struct Same {
  template <Int N>
  static constexpr bool Admissible = true;
  template <Int N>
  static constexpr Int Apply = N;
};

// Conjugation reverses the upper index (theory note section 5). Getting this
// wrong was one of the three structural defects in the layer this replaces.
struct Negate {
  template <Int N>
  static constexpr bool Admissible = true;
  template <Int N>
  static constexpr Int Apply = -N;
};

// abs, abs2, real, imag and Map all land at zero -- which is precisely why
// real and imag are admissible only there, since that is the only place their
// results could be real.
struct Zero {
  template <Int N>
  static constexpr bool Admissible = true;
  template <Int N>
  static constexpr Int Apply = 0;
};

}  // namespace IndexRules

}  // namespace GSHTrans

#endif  // GSH_TRANS_SPIN_WEIGHTED_GUARD_H
