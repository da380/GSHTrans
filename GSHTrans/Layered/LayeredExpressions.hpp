#pragma once

/**
 * @file LayeredExpressions.hpp
 * @brief Lazy expressions over layered fields: a layered node is a map from a
 * radius to a two-dimensional node.
 *
 * @details There is no second expression system. A layered expression has
 * one job, `Slice(i)`: the two-dimensional expression at radius `i`, formed
 * from its operands' slices with the two-dimensional operators. So the
 * upper-index rules, reality, precision and the angular-grid check are the
 * two-dimensional ones, applied as each slice is formed; aliasing is safe for
 * the two-dimensional reason (evaluation is pointwise and index-preserving),
 * slice by slice; and evaluation is the two-dimensional `EvaluateInto`, once
 * per slice, threaded over slices. What this layer adds is the radial axis --
 * operands must share their radial grid as well as their angular one -- and
 * the loop over radii.
 *
 * The operands of a layered expression are layered: owning stacks, views over
 * a caller's storage, other layered expressions, and scalars. A
 * two-dimensional field enters only explicitly, through AtEveryRadius, which
 * makes a layered node that is the same at every radius and stores nothing.
 */

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "../Parallel.hpp"
#include "../Policies.hpp"
#include "../SpinField/SpinFieldOverloads.hpp"
#include "../SpinField/SpinWeighted.hpp"
#include "LayeredSpinField.hpp"
#include "RadialGrid.hpp"

namespace GSHTrans {

// Every type that can stand in a layered expression: the owning stack, its
// views, and the nodes below. (The primary template is in
// LayeredSpinField.hpp, where the owning stack's assignment needs it.)
template <std::ptrdiff_t N, AngularGrid G, RealOrComplexValued V>
struct IsLayeredNodeTrait<LayeredSpinField<N, G, V>> : std::true_type {};
template <std::ptrdiff_t N, AngularGrid G, RealOrComplexValued V, typename E>
struct IsLayeredNodeTrait<LayeredSpinFieldView<N, G, V, E>> : std::true_type {};

/** @brief Whatever can stand in a layered expression. */
template <typename T>
concept LayeredSpinExpr = IsLayeredNodeTrait<std::remove_cvref_t<T>>::value;

namespace LayeredExprDetails {

template <typename T>
struct IsOwningLayered : std::false_type {};
template <std::ptrdiff_t N, AngularGrid G, RealOrComplexValued V>
struct IsOwningLayered<LayeredSpinField<N, G, V>> : std::true_type {};

// How an operand is held, as in the two-dimensional layer: an owning stack
// passed as an lvalue by reference (copying it is the whole stack), and
// anything else -- a view, a node, an rvalue stack -- by value.
template <typename T>
using Stored =
    std::conditional_t<std::is_lvalue_reference_v<T> &&
                           IsOwningLayered<std::remove_cvref_t<T>>::value,
                       const std::remove_cvref_t<T>&, std::remove_cvref_t<T>>;

// The same for the two-dimensional operand of AtEveryRadius.
template <typename T>
using Stored2D =
    std::conditional_t<std::is_lvalue_reference_v<T> &&
                           IsTerminal<std::remove_cvref_t<T>>,
                       const std::remove_cvref_t<T>&, std::remove_cvref_t<T>>;

/// What `Slice(i)` of a layered operand gives.
template <typename T>
using SliceOf = decltype(std::declval<const std::remove_cvref_t<T>&>().Slice(
    std::ptrdiff_t{}));

template <typename A, typename B>
void CheckSameGrids(const A& a, const B& b) {
  if (a.Radial().Identity() != b.Radial().Identity()) {
    throw std::invalid_argument(
        "The operands of a layered expression are on different radial grids");
  }
  if (a.Grid().Identity() != b.Grid().Identity()) {
    throw std::invalid_argument(
        "The operands of a layered expression are on different angular "
        "grids");
  }
}

}  // namespace LayeredExprDetails

/**
 * @brief A layered expression: `Slice(i)` is `make(operand.Slice(i)...)`.
 * @details `Make` builds the two-dimensional node at one radius from the
 * operands' slices; it is one of the two-dimensional operators, or one with a
 * scalar captured. The first operand supplies the grids, and with two
 * operands both grids are checked when the node is made.
 */
template <typename Make, typename... Operands>
class LayeredNode {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  /// The two-dimensional node at one radius.
  using SliceNode = std::remove_cvref_t<std::invoke_result_t<
      const Make&, LayeredExprDetails::SliceOf<Operands>...>>;

  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = SliceNode::UpperIndex;
  using Value = typename SliceNode::Value;        ///< Real- or complex-valued.
  using GridType = typename SliceNode::GridType;  ///< The angular grid type.
  using Real = typename SliceNode::Real;          ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  using Scalar = typename SliceNode::Scalar;  ///< The value type.
  using RadialGridType = RadialGrid<Real>;    ///< The radial grid type.

  /** @brief Combines @p operands, checking that they share their grids. */
  template <typename... Args>
  explicit LayeredNode(Make make, Args&&... operands)
      : make_{std::move(make)}, operands_{std::forward<Args>(operands)...} {
    if constexpr (sizeof...(Operands) == 2) {
      LayeredExprDetails::CheckSameGrids(std::get<0>(operands_),
                                         std::get<1>(operands_));
    }
  }

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return std::get<0>(operands_).Grid(); }
  /** @brief The radial grid this is defined on. */
  const RadialGridType& Radial() const {
    return std::get<0>(operands_).Radial();
  }
  /** @brief How many radii. */
  auto NumberOfRadii() const { return Radial().NumberOfRadii(); }
  /** @brief Indices of the radii. */
  auto RadiusIndices() const { return Radial().RadiusIndices(); }

  /// The two-dimensional expression at radius @p i.
  auto Slice(Int i) const {
    return std::apply(
        [&](const auto&... operand) { return make_(operand.Slice(i)...); },
        operands_);
  }

 private:
  Make make_;
  std::tuple<LayeredExprDetails::Stored<Operands>...> operands_;
};

template <typename Make, typename... Operands>
struct IsLayeredNodeTrait<LayeredNode<Make, Operands...>> : std::true_type {};

/**
 * @brief A two-dimensional expression as a layered one, the same at every
 * radius.
 * @details What a depth-independent weight is in a layered product. Nothing
 * is stored per radius: `Slice(i)` is the two-dimensional expression itself.
 * An owning field passed as an lvalue is held by reference and must outlive
 * the node.
 */
template <typename Expr>
class LayeredBroadcast {
  using E = std::remove_cvref_t<Expr>;

 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = E::UpperIndex;
  using Value = typename E::Value;        ///< Real- or complex-valued.
  using GridType = typename E::GridType;  ///< The angular grid type.
  using Real = typename E::Real;          ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  using Scalar = typename E::Scalar;   ///< The value type.
  using RadialGridType = RadialGrid<Real>;  ///< The radial grid type.

  /** @brief @p expr at every radius of @p radialGrid. */
  template <typename Arg>
  LayeredBroadcast(RadialGridType radialGrid, Arg&& expr)
      : radialGrid_{std::move(radialGrid)}, expr_{std::forward<Arg>(expr)} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return expr_.Grid(); }
  /** @brief The radial grid this is defined on. */
  const RadialGridType& Radial() const { return radialGrid_; }
  /** @brief How many radii. */
  auto NumberOfRadii() const { return radialGrid_.NumberOfRadii(); }
  /** @brief Indices of the radii. */
  auto RadiusIndices() const { return radialGrid_.RadiusIndices(); }

  /// The two-dimensional expression, which is the same at every radius.
  const E& Slice(Int) const { return expr_; }

 private:
  RadialGridType radialGrid_;
  LayeredExprDetails::Stored2D<Expr> expr_;
};

template <typename Expr>
struct IsLayeredNodeTrait<LayeredBroadcast<Expr>> : std::true_type {};

/**
 * @brief A two-dimensional expression, taken as the same at every radius of
 * @p radialGrid, for use in a layered expression.
 * @details The explicit way a two-dimensional operand enters a layered
 * expression; nothing is stored per radius. To make a stack from a
 * two-dimensional field, use Broadcast.
 */
template <SpinFieldExpr Expr>
auto AtEveryRadius(RadialGrid<typename Node<Expr>::Real> radialGrid,
                   Expr&& expr) {
  return LayeredBroadcast<Expr>(std::move(radialGrid),
                                std::forward<Expr>(expr));
}

namespace LayeredExprDetails {

template <typename Make, typename... Args>
auto MakeNode(Make make, Args&&... args) {
  return LayeredNode<Make, Args...>(std::move(make),
                                    std::forward<Args>(args)...);
}

}  // namespace LayeredExprDetails

//--------------------------------------------------------------------------//
//                                 Operators                                 //
//--------------------------------------------------------------------------//

// Each is offered exactly when the two-dimensional operator is defined on the
// slices, so a combination the two-dimensional layer refuses -- unequal upper
// indices in a sum, a real field at N != 0 -- is refused here at overload
// resolution too.

template <LayeredSpinExpr L, LayeredSpinExpr R>
requires requires(LayeredExprDetails::SliceOf<L> a,
                  LayeredExprDetails::SliceOf<R> b) { a + b; }
auto operator+(L&& l, R&& r) {
  return LayeredExprDetails::MakeNode(
      [](auto&& a, auto&& b) {
        return std::forward<decltype(a)>(a) + std::forward<decltype(b)>(b);
      },
      std::forward<L>(l), std::forward<R>(r));
}

template <LayeredSpinExpr L, LayeredSpinExpr R>
requires requires(LayeredExprDetails::SliceOf<L> a,
                  LayeredExprDetails::SliceOf<R> b) { a - b; }
auto operator-(L&& l, R&& r) {
  return LayeredExprDetails::MakeNode(
      [](auto&& a, auto&& b) {
        return std::forward<decltype(a)>(a) - std::forward<decltype(b)>(b);
      },
      std::forward<L>(l), std::forward<R>(r));
}

template <LayeredSpinExpr L, LayeredSpinExpr R>
requires requires(LayeredExprDetails::SliceOf<L> a,
                  LayeredExprDetails::SliceOf<R> b) { a* b; }
auto operator*(L&& l, R&& r) {
  return LayeredExprDetails::MakeNode(
      [](auto&& a, auto&& b) {
        return std::forward<decltype(a)>(a) * std::forward<decltype(b)>(b);
      },
      std::forward<L>(l), std::forward<R>(r));
}

template <LayeredSpinExpr L, LayeredSpinExpr R>
requires requires(LayeredExprDetails::SliceOf<L> a,
                  LayeredExprDetails::SliceOf<R> b) { a / b; }
auto operator/(L&& l, R&& r) {
  return LayeredExprDetails::MakeNode(
      [](auto&& a, auto&& b) {
        return std::forward<decltype(a)>(a) / std::forward<decltype(b)>(b);
      },
      std::forward<L>(l), std::forward<R>(r));
}

// With a scalar: the two-dimensional rule for what may scale a field.
template <LayeredSpinExpr A, typename S>
requires ScalarFor_<S, LayeredExprDetails::SliceOf<A>>
auto operator*(A&& a, S s) {
  return LayeredExprDetails::MakeNode(
      [s](auto&& x) { return std::forward<decltype(x)>(x) * s; },
      std::forward<A>(a));
}

template <typename S, LayeredSpinExpr A>
requires ScalarFor_<S, LayeredExprDetails::SliceOf<A>>
auto operator*(S s, A&& a) {
  return LayeredExprDetails::MakeNode(
      [s](auto&& x) { return s * std::forward<decltype(x)>(x); },
      std::forward<A>(a));
}

template <LayeredSpinExpr A, typename S>
requires ScalarFor_<S, LayeredExprDetails::SliceOf<A>>
auto operator/(A&& a, S s) {
  return LayeredExprDetails::MakeNode(
      [s](auto&& x) { return std::forward<decltype(x)>(x) / s; },
      std::forward<A>(a));
}

template <typename S, LayeredSpinExpr A>
requires ScalarFor_<S, LayeredExprDetails::SliceOf<A>> &&
         requires(S s, LayeredExprDetails::SliceOf<A> x) { s / x; }
auto operator/(S s, A&& a) {
  return LayeredExprDetails::MakeNode(
      [s](auto&& x) { return s / std::forward<decltype(x)>(x); },
      std::forward<A>(a));
}

// Unary: the two-dimensional ones, slice by slice.
template <LayeredSpinExpr A>
auto operator-(A&& a) {
  return LayeredExprDetails::MakeNode(
      [](auto&& x) { return -std::forward<decltype(x)>(x); },
      std::forward<A>(a));
}

/** @brief The complex conjugate, which reverses the upper index. */
template <LayeredSpinExpr A>
auto conj(A&& a) {
  return LayeredExprDetails::MakeNode(
      [](auto&& x) { return conj(std::forward<decltype(x)>(x)); },
      std::forward<A>(a));
}

/** @brief The modulus, at upper index zero. */
template <LayeredSpinExpr A>
requires requires(LayeredExprDetails::SliceOf<A> x) { abs(x); }
auto abs(A&& a) {
  return LayeredExprDetails::MakeNode(
      [](auto&& x) { return abs(std::forward<decltype(x)>(x)); },
      std::forward<A>(a));
}

/** @brief The squared modulus, at upper index zero. */
template <LayeredSpinExpr A>
requires requires(LayeredExprDetails::SliceOf<A> x) { abs2(x); }
auto abs2(A&& a) {
  return LayeredExprDetails::MakeNode(
      [](auto&& x) { return abs2(std::forward<decltype(x)>(x)); },
      std::forward<A>(a));
}

/** @brief The real part, at upper index zero. */
template <LayeredSpinExpr A>
requires requires(LayeredExprDetails::SliceOf<A> x) { real(x); }
auto real(A&& a) {
  return LayeredExprDetails::MakeNode(
      [](auto&& x) { return real(std::forward<decltype(x)>(x)); },
      std::forward<A>(a));
}

/** @brief The imaginary part, at upper index zero. */
template <LayeredSpinExpr A>
requires requires(LayeredExprDetails::SliceOf<A> x) { imag(x); }
auto imag(A&& a) {
  return LayeredExprDetails::MakeNode(
      [](auto&& x) { return imag(std::forward<decltype(x)>(x)); },
      std::forward<A>(a));
}

//--------------------------------------------------------------------------//
//                                Evaluation                                 //
//--------------------------------------------------------------------------//

/**
 * @brief Evaluates @p expr into @p target, every radius.
 * @details Threaded over radii under a parallel policy; each radius is the
 * two-dimensional evaluation. @p target may appear in @p expr: evaluation is
 * pointwise and index-preserving.
 * @param target An owning stack or a writable view, on @p expr's grids.
 * @param expr The layered expression.
 * @param policy Whether to thread over radii.
 * @throws std::invalid_argument if the two are not on the same grids.
 */
template <LayeredFieldStorage Target, LayeredSpinExpr Expr>
requires EvaluatesInto<LayeredExprDetails::SliceOf<Expr>,
                       std::remove_cvref_t<Target>::UpperIndex,
                       typename std::remove_cvref_t<Target>::Real,
                       typename std::remove_cvref_t<Target>::Scalar> &&
         requires(std::remove_cvref_t<Target>& t) {
           requires not std::is_const_v<
               typename decltype(t.Data())::element_type>;
         }
void Assign(Target&& target, const Expr& expr,
            Execution policy = Execution::Sequential()) {
  LayeredExprDetails::CheckSameGrids(target, expr);
  ForEach(policy, target.NumberOfRadii(), [&](std::ptrdiff_t i) {
    auto slice = target.Slice(i);
    expr.Slice(i).EvaluateInto(slice.Data());
  });
}

/**
 * @brief A new stack holding the value of @p expr.
 * @param expr The layered expression.
 * @param policy Whether to thread over radii.
 */
template <LayeredSpinExpr Expr>
auto Materialise(const Expr& expr, Execution policy = Execution::Sequential()) {
  using E = std::remove_cvref_t<Expr>;
  auto stack =
      LayeredSpinField<E::UpperIndex, typename E::GridType, typename E::Value>(
          expr.Radial(), expr.Grid());
  Assign(stack, expr, policy);
  return stack;
}

//--------------------------------------------------------------------------//
//                        Expanding a layered expression                     //
//--------------------------------------------------------------------------//

namespace LayeredExprDetails {

// How many radii an expression is evaluated for at a time before they are
// transformed: about this many, rounded to a multiple of the transform's own
// chunk.
inline constexpr std::ptrdiff_t ExpandBlockTarget = 64;

// The scratch a block is evaluated into: per calling thread, kept between
// calls and grown only, so a repeated Expand allocates nothing after the
// first.
template <typename Scalar>
std::vector<Scalar>& ExpandScratch(std::size_t size) {
  thread_local auto buffer = std::vector<Scalar>{};
  if (buffer.size() < size) buffer.resize(size);
  return buffer;
}

}  // namespace LayeredExprDetails

/**
 * @brief Expands a layered expression into @p out, to @p out's degree,
 * without its value ever being stored whole.
 * @details The expression is evaluated a block of radii at a time into
 * per-thread scratch and each block is transformed with one batched call, so
 * the scratch is a block of slices and not a stack. The block is a multiple
 * of the transform's own chunk (SphericalGrid::ForwardChunk), so every chunk
 * falls where one call over the stored value would put it, and the result is
 * that call's bit for bit. Evaluation is threaded over the block's radii and
 * the transform threads within itself, one after the other.
 * @param expr A layered expression.
 * @param out A layered expansion or writable view, on @p expr's grids.
 * @param policy Whether to thread.
 * @throws std::invalid_argument if the two are not on the same grids.
 */
template <LayeredSpinExpr Expr, LayeredExpansionStorage Out>
requires(!LayeredFieldStorage<Expr>) &&
        LayeredDetails::WritableLike<Out, std::remove_cvref_t<Expr>>
void Expand(const Expr& expr, Out&& out,
            Execution policy = Execution::Sequential()) {
  using E = std::remove_cvref_t<Expr>;
  using Scalar = typename E::Scalar;
  constexpr auto N = E::UpperIndex;
  LayeredDetails::CheckSameGrids(expr, out);

  const auto& grid = expr.Grid();
  const auto nR = static_cast<std::ptrdiff_t>(expr.NumberOfRadii());
  const auto fieldSize = static_cast<std::ptrdiff_t>(grid.FieldSize());
  const auto lMax = static_cast<std::ptrdiff_t>(out.MaxDegree());
  const auto coefficientSize =
      static_cast<std::ptrdiff_t>(out.CoefficientSize());
  const auto chunk = grid.template ForwardChunk<Scalar>(lMax, N, policy);
  const auto block = std::min(
      nR, chunk * std::max(std::ptrdiff_t{1},
                           LayeredExprDetails::ExpandBlockTarget / chunk));

  auto& scratch = LayeredExprDetails::ExpandScratch<Scalar>(
      static_cast<std::size_t>(block * fieldSize));
  auto target = out.Data();
  for (auto first = std::ptrdiff_t{0}; first < nR; first += block) {
    const auto count = std::min(block, nR - first);
    auto samples = std::span<Scalar>(
        scratch.data(), static_cast<std::size_t>(count * fieldSize));
    ForEach(policy, count, [&](std::ptrdiff_t i) {
      expr.Slice(first + i).EvaluateInto(
          samples.subspan(static_cast<std::size_t>(i * fieldSize),
                          static_cast<std::size_t>(fieldSize)));
    });
    auto coefficients =
        target.subspan(static_cast<std::size_t>(first * coefficientSize),
                       static_cast<std::size_t>(count * coefficientSize));
    grid.ForwardTransformation(
        lMax, N, std::span<const Scalar>(samples),
        Batch::Contiguous(count, fieldSize), coefficients,
        Batch::Contiguous(count, coefficientSize), policy);
  }
}

/** @brief Expands a layered expression to degree @p lMax. */
template <LayeredSpinExpr Expr>
requires(!LayeredFieldStorage<Expr>)
auto Expand(const Expr& expr, std::ptrdiff_t lMax,
            Execution policy = Execution::Sequential()) {
  using E = std::remove_cvref_t<Expr>;
  auto expansion =
      LayeredSpinExpansion<E::UpperIndex, typename E::GridType,
                           typename E::Value>(expr.Radial(), expr.Grid(), lMax);
  Expand(expr, expansion, policy);
  return expansion;
}

}  // namespace GSHTrans
