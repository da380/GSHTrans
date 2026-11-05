#pragma once

/**
 * @file LayeredTensorExpressions.hpp
 * @brief Tensor expressions over layered tensors: evaluation, materialisation
 * and expansion, and the scale by a layered scalar field.
 *
 * @details A tensor expression is a rank, a grid, an alphabet and a component
 * accessor (TensorExpr.hpp). A layered tensor's components are layered
 * expressions (LayeredTensorField::Component), and the tensor operations --
 * Permute, TensorProduct, Contract, Symmetrise -- form their components with
 * operators that layered expressions have, so applied to layered tensors they
 * give *layered tensor expressions* with no node of their own. What this file
 * adds is what is specific to the layered case: a layered scalar weight,
 * evaluating into a layered tensor, and expanding without storing the value.
 * Each works one stored component at a time through the layered spin layer,
 * so it threads over radii and expands in blocks as that layer does.
 */

#include <algorithm>
#include <complex>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "../Policies.hpp"
#include "../Tensor/MultiIndex.hpp"
#include "../Tensor/Orbits.hpp"
#include "../Tensor/TensorExpr.hpp"
#include "LayeredExpressions.hpp"
#include "LayeredSpinField.hpp"
#include "LayeredTensorField.hpp"

namespace GSHTrans {

/// A tensor expression whose components are layered expressions: a
/// LayeredTensorField, or a tensor operation on layered tensors.
template <typename T>
concept LayeredTensorExpr =
    TensorExpr<std::remove_cvref_t<T>> &&
    (TensorDetails::FirstRepresentedFlat<std::remove_cvref_t<T>>() >= 0) &&
    LayeredSpinExpr<TensorDetails::ComponentType<std::remove_cvref_t<T>>>;

namespace LayeredTensorDetails {

template <typename T>
struct IsLayeredTensorField : std::false_type {};
template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid, SlotAlphabet Slots,
          template <typename, LayeredDetails::Int> typename At>
struct IsLayeredTensorField<
    LayeredTensorField<Rank, Symmetry, Reality, Grid, Slots, At>>
    : std::true_type {};

// The radii of a layered tensor expression, from its first component.
template <typename Expr>
auto RadialOf(const Expr& expr) {
  using E = std::remove_cvref_t<Expr>;
  return TensorDetails::ComponentOf<TensorDetails::FirstRepresented<E>>(expr)
      .Radial();
}

// A component as the stack it is written into wants it. A component the
// reality condition pins to one real number, but which the expression gives
// as complex, contributes the part that survives -- the real part, or the
// coefficient of i for an imaginary orbit -- as flat Materialise does; that is
// the caller's assertion that the expression is a real tensor. A real
// component going into a complex stack is promoted only where a transform
// needs the value kinds to match.
template <ComponentConstraint constraint, typename TargetValue,
          bool PromoteReal, typename Source>
auto AsTargetValue(Source&& source) {
  using S = std::remove_cvref_t<Source>;
  using Real = typename S::Real;
  if constexpr (std::same_as<TargetValue, RealValued> &&
                std::same_as<typename S::Value, ComplexValued>) {
    if constexpr (constraint == ComponentConstraint::Imaginary) {
      return imag(std::forward<Source>(source));
    } else {
      return real(std::forward<Source>(source));
    }
  } else if constexpr (PromoteReal &&
                       std::same_as<TargetValue, ComplexValued> &&
                       std::same_as<typename S::Value, RealValued>) {
    return std::complex<Real>{1, 0} * std::forward<Source>(source);
  } else {
    return std::forward<Source>(source);
  }
}

// The stack of the component at the given slot letters.
template <auto Indices, typename Target, std::size_t... I>
auto& StackOf(Target& target, std::index_sequence<I...>) {
  return target.template ComponentStack<Indices[I]...>();
}

// Run `apply(indices, stack)` for every stored component of a layered tensor
// field or expansion, with the stack and the component's slot letters.
template <typename Target, typename Apply>
void ForEachStoredSlot(Target& target, Apply&& apply) {
  using T = std::remove_cvref_t<Target>;
  [&]<std::size_t... Slot>(std::index_sequence<Slot...>) {
    (
        [&] {
          constexpr auto flat = T::ComponentLayout.flatOfSlot[Slot];
          constexpr auto indices =
              MultiIndex<T::Rank, typename T::SlotSet>::FromFlat(flat).Slots();
          constexpr auto constraint = T::Orbits.constraint[flat];
          auto& stack = StackOf<indices>(
              target,
              std::make_index_sequence<static_cast<std::size_t>(T::Rank)>{});
          apply.template operator()<indices, constraint>(stack);
        }(),
        ...);
  }(std::make_index_sequence<static_cast<std::size_t>(T::StoredComponents)>{});
}

}  // namespace LayeredTensorDetails

/// A layered tensor times a layered scalar field (upper index zero).
template <typename W, typename T>
requires LayeredSpinExpr<W> &&
         (std::remove_cvref_t<W>::UpperIndex == 0) && LayeredTensorExpr<T>
auto operator*(W&& weight, T&& tensor) {
  return ScaleNode<W, T>(std::forward<W>(weight), std::forward<T>(tensor));
}

/// The same, written the other way round.
template <typename T, typename W>
requires LayeredSpinExpr<W> &&
         (std::remove_cvref_t<W>::UpperIndex == 0) && LayeredTensorExpr<T>
auto operator*(T&& tensor, W&& weight) {
  return ScaleNode<W, T>(std::forward<W>(weight), std::forward<T>(tensor));
}

/**
 * @brief Evaluates a layered tensor expression into @p target.
 * @details Every component @p target stores is evaluated from @p expr, one
 * stored component at a time and threaded over radii; a component @p expr
 * does not represent is zero. The symmetry and reality are @p target's: as in
 * flat Materialise, choosing them asserts them of the value.
 * @throws std::invalid_argument if the grids differ.
 */
template <LayeredTensorExpr Expr, std::ptrdiff_t Rank,
          TensorSymmetry<Rank> Symmetry, TensorReality Reality,
          AngularGrid Grid, SlotAlphabet Slots,
          template <typename, LayeredDetails::Int> typename At>
requires(std::remove_cvref_t<Expr>::Rank == Rank) &&
        std::same_as<typename std::remove_cvref_t<Expr>::SlotSet, Slots>
void Assign(
    LayeredTensorField<Rank, Symmetry, Reality, Grid, Slots, At>& target,
    const Expr& expr, Execution policy = Execution::Sequential()) {
  using E = std::remove_cvref_t<Expr>;
  LayeredTensorDetails::ForEachStoredSlot(
      target, [&]<auto indices, ComponentConstraint constraint>(auto& stack) {
        using Stack = std::remove_cvref_t<decltype(stack)>;
        if constexpr (TensorDetails::Represents<indices, E>()) {
          Assign(
              stack,
              LayeredTensorDetails::AsTargetValue<constraint,
                                                  typename Stack::Value, false>(
                  TensorDetails::ComponentOf<indices>(expr)),
              policy);
        } else {
          std::ranges::fill(stack.Data(), typename Stack::Scalar{});
        }
      });
}

/**
 * @brief A new layered tensor holding the value of @p expr, with the symmetry
 * and reality the caller states (none and complex by default).
 */
template <typename Symmetry = void, TensorReality Reality = ComplexTensor,
          LayeredTensorExpr Expr>
auto Materialise(const Expr& expr, Execution policy = Execution::Sequential()) {
  using E = std::remove_cvref_t<Expr>;
  using Chosen = std::conditional_t<std::same_as<Symmetry, void>,
                                    NoSymmetry<E::Rank>, Symmetry>;
  auto field = LayeredTensorField<E::Rank, Chosen, Reality,
                                  typename E::GridType, typename E::SlotSet>(
      LayeredTensorDetails::RadialOf(expr), expr.Grid());
  Assign(field, expr, policy);
  return field;
}

/**
 * @brief Expands a layered tensor expression into @p out, to @p out's degree,
 * without storing its value.
 * @details One expansion per component @p out stores, each a layered spin
 * Expand: a block of radii evaluated and transformed at a time, bitwise equal
 * to expanding the stored value. A component @p expr does not represent is
 * zero.
 * @throws std::invalid_argument if the grids differ.
 */
template <LayeredTensorExpr Expr, std::ptrdiff_t Rank,
          TensorSymmetry<Rank> Symmetry, TensorReality Reality,
          AngularGrid Grid, SlotAlphabet Slots,
          template <typename, LayeredDetails::Int> typename At>
requires(!LayeredTensorDetails::IsLayeredTensorField<
            std::remove_cvref_t<Expr>>::value) &&
        (std::remove_cvref_t<Expr>::Rank == Rank) &&
        std::same_as<typename std::remove_cvref_t<Expr>::SlotSet, Slots>
void Expand(
    const Expr& expr,
    LayeredTensorExpansion<Rank, Symmetry, Reality, Grid, Slots, At>& out,
    Execution policy = Execution::Sequential()) {
  using E = std::remove_cvref_t<Expr>;
  LayeredTensorDetails::ForEachStoredSlot(
      out, [&]<auto indices, ComponentConstraint constraint>(auto& stack) {
        using Stack = std::remove_cvref_t<decltype(stack)>;
        if constexpr (TensorDetails::Represents<indices, E>()) {
          Expand(
              LayeredTensorDetails::AsTargetValue<constraint,
                                                  typename Stack::Value, true>(
                  TensorDetails::ComponentOf<indices>(expr)),
              stack, policy);
        } else {
          std::ranges::fill(stack.Data(), std::complex<typename Grid::Real>{});
        }
      });
}

/**
 * @brief Expands a layered tensor expression to degree @p lMax, with the
 * symmetry and reality the caller states (none and complex by default).
 */
template <typename Symmetry = void, TensorReality Reality = ComplexTensor,
          LayeredTensorExpr Expr>
requires(!LayeredTensorDetails::IsLayeredTensorField<
         std::remove_cvref_t<Expr>>::value)
auto Expand(const Expr& expr, std::ptrdiff_t lMax,
            Execution policy = Execution::Sequential()) {
  using E = std::remove_cvref_t<Expr>;
  using Chosen = std::conditional_t<std::same_as<Symmetry, void>,
                                    NoSymmetry<E::Rank>, Symmetry>;
  auto expansion =
      LayeredTensorExpansion<E::Rank, Chosen, Reality, typename E::GridType,
                             typename E::SlotSet>(
          LayeredTensorDetails::RadialOf(expr), expr.Grid(), lMax);
  Expand(expr, expansion, policy);
  return expansion;
}

}  // namespace GSHTrans
