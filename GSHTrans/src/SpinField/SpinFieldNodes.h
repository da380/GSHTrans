#ifndef GSH_TRANS_SPIN_FIELD_NODES_GUARD_H
#define GSH_TRANS_SPIN_FIELD_NODES_GUARD_H

#include <complex>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "../Concepts.h"
#include "SpinWeighted.h"

namespace GSHTrans {

// The two expression node templates. Every operator in this layer builds one
// of these; there is nothing else.
//
// Both are pointwise *and index-preserving*: node (iTheta, iPhi) reads its
// operands only at (iTheta, iPhi). That is an invariant, not an accident. It
// is what makes in-place assignment safe without a temporary, and it is the
// property gradients, raising and lowering do not have -- which is why those
// live in the spectral layer instead. A re-indexing view, a phi shift or a
// transpose would be pointwise in the loose sense and would break the
// argument, so the stronger form is the one stated.

//--------------------------------------------------------------------------//
//                                  Binary                                   //
//--------------------------------------------------------------------------//

/// L and R are the *deduced* operand types and so carry the caller's value
/// category: an lvalue arrives as T&, an rvalue as T. OperandStorage turns that
/// into what is actually held.
template <typename Op, typename Rule, typename L, typename R>
class Binary {
  using LNode = Node<L>;
  using RNode = Node<R>;

 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.
  using Real = typename LNode::Real;  ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  using GridType = typename LNode::GridType;  ///< The angular grid this is defined on.

  /// The scalar is whatever the operation returns, and the value kind follows
  /// from it. Deriving it rather than declaring it is what makes the reality
  /// constraint hold by construction instead of by argument.
  using Scalar =
      std::invoke_result_t<Op, typename LNode::Scalar, typename RNode::Scalar>;
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;

  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex =
      Rule::template Apply<LNode::UpperIndex, RNode::UpperIndex>;

  static_assert(SpinWeighted<LNode> and SpinWeighted<RNode>);
  static_assert(std::same_as<typename LNode::GridType,
                             typename RNode::GridType>,
                "both operands must be on the same kind of grid");
  static_assert(std::same_as<typename LNode::Real, typename RNode::Real>,
                "one precision per expression tree: no mixed-precision "
                "promotion");
  static_assert(Rule::template Admissible<LNode::UpperIndex,
                                          RNode::UpperIndex>,
                "the upper indices of these operands do not satisfy this "
                "operation's index rule");
  static_assert(RealOrComplexFloatingPoint<Scalar>);
  static_assert(std::same_as<RemoveComplex<Scalar>, Real>,
                "the operation changed the precision of its operands");
  static_assert(std::same_as<Value, ComplexValued> or UpperIndex == 0,
                "closure lemma violated: a real-valued result away from upper "
                "index zero");

  Binary(L&& l, R&& r, Op op = Op{})
      : _l{std::forward<L>(l)}, _r{std::forward<R>(r)}, _op{std::move(op)} {
    // Handle identity, in all build modes. Two separately built grids with
    // equal parameters are not the same grid: their samples are different
    // arrays and pairing them elementwise is meaningless. Before this rewrite
    // the grid was taken from the left operand alone and a mismatch was
    // silent.
    if (_l.Grid().Identity() != _r.Grid().Identity()) {
      throw std::invalid_argument(
          "The operands of a binary spin-field expression are on different "
          "grids");
    }
  }

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _l.Grid(); }

  Scalar operator[](Int iTheta, Int iPhi) const {
    return _op(_l[iTheta, iPhi], _r[iTheta, iPhi]);
  }

  template <typename S>
  requires std::convertible_to<Scalar, S>
  void EvaluateInto(std::span<S> target) const {
    EvaluateNodeInto(*this, target);
  }

 private:
  OperandStorage<L> _l;
  OperandStorage<R> _r;
  [[no_unique_address]] Op _op;
};

//--------------------------------------------------------------------------//
//                                   Unary                                   //
//--------------------------------------------------------------------------//

template <typename Op, typename Rule, typename A>
class Unary {
  using ANode = Node<A>;

 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.
  using Real = typename ANode::Real;  ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  using GridType = typename ANode::GridType;  ///< The angular grid this is defined on.

  using Scalar = std::invoke_result_t<Op, typename ANode::Scalar>;  ///< The value type: Real when real-valued, Complex otherwise.
  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;

  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = Rule::template Apply<ANode::UpperIndex>;

  static_assert(SpinWeighted<ANode>);
  static_assert(Rule::template Admissible<ANode::UpperIndex>,
                "the upper index of this operand does not satisfy this "
                "operation's index rule");
  static_assert(RealOrComplexFloatingPoint<Scalar>);
  static_assert(std::same_as<RemoveComplex<Scalar>, Real>,
                "the operation changed the precision of its operand");
  static_assert(std::same_as<Value, ComplexValued> or UpperIndex == 0,
                "closure lemma violated: a real-valued result away from upper "
                "index zero");

  Unary(A&& a, Op op = Op{})
      : _a{std::forward<A>(a)}, _op{std::move(op)} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _a.Grid(); }

  Scalar operator[](Int iTheta, Int iPhi) const {
    return _op(_a[iTheta, iPhi]);
  }

  template <typename S>
  requires std::convertible_to<Scalar, S>
  void EvaluateInto(std::span<S> target) const {
    EvaluateNodeInto(*this, target);
  }

 private:
  OperandStorage<A> _a;
  [[no_unique_address]] Op _op;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SPIN_FIELD_NODES_GUARD_H
