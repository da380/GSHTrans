#ifndef GSH_TRANS_LAYERED_TENSOR_FIELD_GUARD_H
#define GSH_TRANS_LAYERED_TENSOR_FIELD_GUARD_H

#include <complex>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

#include "../Concepts.h"
#include "../Expansion/SpinExpansion.h"
#include "../Expansion/TensorExpansion.h"
#include "../Policies.h"
#include "../Tensor/Orbits.h"
#include "../Tensor/TensorField.h"
#include "../Utility.h"
#include "LayeredSpinField.h"
#include "RadialGrid.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                         Layered tensors: the layout                       //
//--------------------------------------------------------------------------//

// A tensor field on the radial-times-angular product, held as one radius-major
// stack per stored component.
//
// The flat tensor packs every component into one buffer, ordered by upper index
// so that the components sharing one are a batch. A layered tensor does not
// need that: the radial axis already supplies `nR` fields per component, so
// each component is its own batch and there is nothing to gain by interleaving
// them. Keeping the components apart instead buys three things.
//
//   - Each component's radial stack is exactly a LayeredSpinField, contiguous
//     and directly transformable, with no offset arithmetic and no second copy
//     of the flat type's block layout to keep in step with it.
//   - It is the layout the applications this exists for already use, which
//     hold one radius-major array per component. Handing back the stack is
//     handing back their array.
//   - The real and complex components separate for free. A component the
//     reality condition pins to one real number per point lives in a
//     RealValued stack, which is the same split the flat tensor makes and for
//     the same reason.
//
// What it costs is that the tensor is not one buffer, so there is no single
// span over the whole thing. Nothing here wants one: every operation is per
// component, because a transform is.
//
// Components are reached by multi-index as on the flat type, with the same two
// accessors and the same conditions. Slice(r) returning a *tensor* view is
// deliberately not offered -- there is no tensor view type, and reaching the
// component directly is the smaller of the two ways to do without one.

namespace LayeredDetails {

using Int = std::ptrdiff_t;

template <typename Flat, Int Slot>
using ValueAt = std::conditional_t<Flat::ComponentLayout.realOfSlot[Slot],
                                   RealValued, ComplexValued>;

template <typename Flat, Int Slot>
using FieldStackAt =
    LayeredSpinField<Flat::ComponentLayout.upperIndexOfSlot[Slot],
                     typename Flat::GridType, ValueAt<Flat, Slot>>;

template <typename Flat, Int Slot>
using ExpansionStackAt =
    LayeredSpinExpansion<Flat::ComponentLayout.upperIndexOfSlot[Slot],
                         typename Flat::GridType, ValueAt<Flat, Slot>>;

template <template <typename, Int> typename At, typename Flat, typename Seq>
struct StacksImpl;

template <template <typename, Int> typename At, typename Flat,
          std::size_t... Slots>
struct StacksImpl<At, Flat, std::index_sequence<Slots...>> {
  using type = std::tuple<At<Flat, static_cast<Int>(Slots)>...>;
};

template <template <typename, Int> typename At, typename Flat>
using Stacks =
    typename StacksImpl<At, Flat,
                        std::make_index_sequence<static_cast<std::size_t>(
                            Flat::StoredComponents)>>::type;

}  // namespace LayeredDetails

//--------------------------------------------------------------------------//
//                          The layered tensor field                         //
//--------------------------------------------------------------------------//

/**
 * @brief A tensor field on a radial stack of spheres: one LayeredSpinField
 * per stored component.
 *
 * @details The combinatorics are the flat tensor's, which this type does not
 * repeat. It adds a radial axis and nothing else, so a component is reached
 * by multi-index exactly as on the flat type and under the same conditions.
 *
 * @tparam _Rank The tensor rank.
 * @tparam _Symmetry The permutation symmetry of the slots.
 * @tparam _Reality Whether the tensor is real or complex.
 * @tparam _Grid The angular grid.
 * @tparam _Slots The alphabet the slots are drawn from.
 */
template <std::ptrdiff_t _Rank, TensorSymmetry<_Rank> _Symmetry,
          TensorReality _Reality, AngularGrid _Grid,
          SlotAlphabet _Slots = AllSlots>
class LayeredTensorField {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.

  /** @brief The tensor rank. */
  static constexpr Int Rank = _Rank;
  using Symmetry = _Symmetry;  ///< The permutation symmetry of the slots.
  using Reality = _Reality;    ///< Whether the tensor is real or complex.
  using GridType = _Grid;      ///< The angular grid this is defined on.
  using Real = typename _Grid::Real;   ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  using RadialGridType = RadialGrid<Real>;  ///< The radial grid type.
  using SlotSet = _Slots;  ///< The alphabet the slots are drawn from.

  /// The flat tensor of the same shape, which owns the combinatorics: the
  /// orbits, the slot layout and the conditions on the accessors. This type
  /// adds a radial axis and nothing else, so it defines none of that itself --
  /// including which multi-indices exist, which is why the alphabet is passed
  /// through rather than consulted here. The accessors inherit the letter
  /// check with everything else, since they are conditioned on Flat's
  /// Represents and Writable.
  using Flat =
      TensorField<Rank, Symmetry, Reality, GridType, ComponentMajor, SlotSet>;

  /** @brief The orbits of the symmetry group, and what each pins. */
  static constexpr auto& Orbits = Flat::Orbits;
  /** @brief Where each stored component sits, taken from the flat tensor. */
  static constexpr auto& ComponentLayout = Flat::ComponentLayout;
  /** @brief How many components the tensor has, stored or not. */
  static constexpr Int Components = Flat::Components;
  /** @brief How many are actually stored, one per orbit. */
  static constexpr Int StoredComponents = Flat::StoredComponents;

  /** @brief Whether this holds a component at those slot letters. */
  template <Int... Alphas>
  static constexpr bool Represents = Flat::template Represents<Alphas...>;

  /** @brief Whether that component may be written through. */
  template <Int... Alphas>
  static constexpr bool Writable = Flat::template Writable<Alphas...>;

  /** @brief Whether that component is identically zero. */
  template <Int... Alphas>
  static constexpr bool Vanishes = Flat::template Vanishes<Alphas...>;

  LayeredTensorField() = delete;

  /**
   * @brief A zero tensor field on the product of the two grids.
   * @param radialGrid The radii.
   * @param grid The angular grid, which must carry every upper index from
   * -Rank to Rank.
   */
  LayeredTensorField(RadialGridType radialGrid, GridType grid)
      : _stacks{Build(radialGrid, grid,
                      std::make_index_sequence<static_cast<std::size_t>(
                          StoredComponents)>{})},
        _radialGrid{std::move(radialGrid)},
        _grid{std::move(grid)} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _grid; }
  /** @brief The radial grid this is defined on. */
  const RadialGridType& Radial() const { return _radialGrid; }

  /** @brief How many radii the stack holds. */
  auto NumberOfRadii() const { return _radialGrid.NumberOfRadii(); }
  /** @brief Indices of the stored radii. */
  auto RadiusIndices() const { return _radialGrid.RadiusIndices(); }
  /** @brief How many samples one angular field holds. */
  auto FieldSize() const { return static_cast<Int>(_grid.FieldSize()); }

  /// The whole radial stack of one stored component: a LayeredSpinField, and
  /// therefore something Expand, ApplyRadially and IntegrateRadially already
  /// take. This is the accessor most callers want.
  template <Int... Alphas>
  requires Writable<Alphas...>
  auto& ComponentStack() {
    return std::get<static_cast<std::size_t>(SlotOf<Alphas...>())>(_stacks);
  }

  /// The same, read-only.
  template <Int... Alphas>
  requires Writable<Alphas...>
  const auto& ComponentStack() const {
    return std::get<static_cast<std::size_t>(SlotOf<Alphas...>())>(_stacks);
  }

  /// One angular slice of a stored component, as a writable view.
  template <Int... Alphas>
  requires Writable<Alphas...>
  auto Component(Int i) {
    return ComponentStack<Alphas...>().Slice(i);
  }

  /// One angular slice of *any* representable component. A derived one is its
  /// representative's, with the sign and conjugation the orbit table records --
  /// the same relation the flat type applies, applied to a slice.
  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component(Int i) const {
    constexpr auto flat = Flat::template FlatOf<Alphas...>;
    constexpr auto sign = Orbits.sign[flat];
    constexpr auto conjugated = Orbits.conjugate[flat];
    constexpr auto constraint = Orbits.constraint[flat];
    constexpr auto scale = static_cast<Real>(sign);

    auto view = SliceOfSlot<SlotOfFlat(Orbits.representative[flat])>(i);

    if constexpr (constraint == ComponentConstraint::None) {
      if constexpr (conjugated) {
        return scale * conj(std::move(view));
      } else if constexpr (sign == 1) {
        return view;
      } else {
        return -std::move(view);
      }
    } else {
      // A pinned component is one real number per point: Real means the value
      // is that number, Imaginary that it is i times it. Conjugation turns the
      // imaginary case's sign, as on the flat type.
      constexpr auto turn = conjugated ? -scale : scale;
      if constexpr (constraint == ComponentConstraint::Real) {
        return scale * view;
      } else {
        return turn * view;
      }
    }
  }

 private:
  LayeredDetails::Stacks<LayeredDetails::FieldStackAt, Flat> _stacks;
  RadialGridType _radialGrid;
  GridType _grid;

  template <Int... Alphas>
  static constexpr Int SlotOf() {
    return Flat::SlotOfFlat(Flat::template FlatOf<Alphas...>);
  }

  static constexpr Int SlotOfFlat(Int flat) { return Flat::SlotOfFlat(flat); }

  template <Int Slot>
  auto SliceOfSlot(Int i) const {
    return std::get<static_cast<std::size_t>(Slot)>(_stacks).Slice(i);
  }

  template <std::size_t... Slots>
  static auto Build(const RadialGridType& radialGrid, const GridType& grid,
                    std::index_sequence<Slots...>) {
    return LayeredDetails::Stacks<LayeredDetails::FieldStackAt, Flat>{
        LayeredDetails::FieldStackAt<Flat, static_cast<Int>(Slots)>(radialGrid,
                                                                    grid)...};
  }
};

//--------------------------------------------------------------------------//
//                        The layered tensor expansion                       //
//--------------------------------------------------------------------------//

/**
 * @brief The spectral counterpart of LayeredTensorField: one
 * LayeredSpinExpansion per stored component.
 *
 * @copydetails LayeredTensorField
 */
template <std::ptrdiff_t _Rank, TensorSymmetry<_Rank> _Symmetry,
          TensorReality _Reality, AngularGrid _Grid,
          SlotAlphabet _Slots = AllSlots>
class LayeredTensorExpansion {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.

  /** @brief The tensor rank. */
  static constexpr Int Rank = _Rank;
  using Symmetry = _Symmetry;  ///< The permutation symmetry of the slots.
  using Reality = _Reality;    ///< Whether the tensor is real or complex.
  using GridType = _Grid;      ///< The angular grid this is defined on.
  using Real = typename _Grid::Real;   ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  using RadialGridType = RadialGrid<Real>;  ///< The radial grid type.
  using SlotSet = _Slots;  ///< The alphabet the slots are drawn from.

  /** @brief The flat tensor of the same shape, which owns the
   * combinatorics. */
  using Flat =
      TensorField<Rank, Symmetry, Reality, GridType, ComponentMajor, SlotSet>;
  /** @brief The spatial tensor this expands. */
  using FieldType =
      LayeredTensorField<Rank, Symmetry, Reality, GridType, SlotSet>;

  /** @brief The orbits of the symmetry group, and what each pins. */
  static constexpr auto& Orbits = Flat::Orbits;
  /** @brief Where each stored component sits, taken from the flat tensor. */
  static constexpr auto& ComponentLayout = Flat::ComponentLayout;
  /** @brief How many components are actually stored, one per orbit. */
  static constexpr Int StoredComponents = Flat::StoredComponents;

  /** @brief Whether this holds a component at those slot letters. */
  template <Int... Alphas>
  static constexpr bool Represents = Flat::template Represents<Alphas...>;

  /** @brief Whether that component may be written through. */
  template <Int... Alphas>
  static constexpr bool Writable = Flat::template Writable<Alphas...>;

  LayeredTensorExpansion() = delete;

  /**
   * @brief A zero expansion on the product of the two grids.
   * @param radialGrid The radii.
   * @param grid The angular grid the coefficients belong to.
   * @param lMax The largest degree stored.
   */
  LayeredTensorExpansion(RadialGridType radialGrid, GridType grid, Int lMax)
      : _stacks{Build(radialGrid, grid, lMax,
                      std::make_index_sequence<static_cast<std::size_t>(
                          StoredComponents)>{})},
        _radialGrid{std::move(radialGrid)},
        _grid{std::move(grid)},
        _lMax{lMax} {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return _grid; }
  /** @brief The radial grid this is defined on. */
  const RadialGridType& Radial() const { return _radialGrid; }

  /** @brief How many radii the stack holds. */
  auto NumberOfRadii() const { return _radialGrid.NumberOfRadii(); }
  /** @brief Indices of the stored radii. */
  auto RadiusIndices() const { return _radialGrid.RadiusIndices(); }
  /** @brief The largest degree stored. */
  auto MaxDegree() const { return _lMax; }

  /// The whole radial stack of one stored component's coefficients, as a
  /// LayeredSpinExpansion.
  template <Int... Alphas>
  requires Writable<Alphas...>
  auto& ComponentStack() {
    return std::get<static_cast<std::size_t>(SlotOf<Alphas...>())>(_stacks);
  }

  /// The same, read-only.
  template <Int... Alphas>
  requires Writable<Alphas...>
  const auto& ComponentStack() const {
    return std::get<static_cast<std::size_t>(SlotOf<Alphas...>())>(_stacks);
  }

  /// The coefficient of any representable component, at one radius.
  ///
  /// This is TensorExpansion::Coefficient with a radius index threaded
  /// through it, and the same four cases: stored, permutation-relative,
  /// reality-relative, or an orbit that vanishes. Degrees below the
  /// component's own |N| return zero, which is what lets the gradient read
  /// shifted components without special-casing the ones that do not exist.
  template <Int... Alphas>
  requires(sizeof...(Alphas) == Rank)
  Complex Coefficient(Int i, Int l, Int m) const {
    constexpr auto flat = Flat::template FlatOf<Alphas...>;
    constexpr auto n = Flat::template UpperIndexOf<Alphas...>;
    constexpr auto constraint = Orbits.constraint[flat];

    if constexpr (constraint == ComponentConstraint::Zero) {
      return Complex{};
    } else {
      if (l < (n < 0 ? -n : n) || l > _lMax || m < -l || m > l) {
        return Complex{};
      }

      constexpr auto rep = Orbits.representative[flat];
      constexpr auto slot = Flat::SlotOfFlat(rep);
      constexpr auto repN = ComponentLayout.upperIndexOfSlot[slot];
      constexpr auto sign = static_cast<Real>(Orbits.sign[flat]);
      constexpr auto conjugated = Orbits.conjugate[flat];
      constexpr auto real = ComponentLayout.realOfSlot[slot];

      const auto& block = std::get<static_cast<std::size_t>(slot)>(_stacks);

      // A real block holds only m >= 0; the rest is f_{l,-m} = (-1)^m
      // conj(f_{lm}).
      const auto stored = [&](Int order) {
        if constexpr (real) {
          if (order < 0) {
            return static_cast<Real>(MinusOneToPower(order)) *
                   std::conj(block[i, l, -order]);
          }
        }
        return Complex{block[i, l, order]};
      };

      constexpr auto turn = constraint == ComponentConstraint::Imaginary
                                ? Complex{0, 1}
                                : Complex{1, 0};

      if constexpr (conjugated) {
        return sign * turn * static_cast<Real>(MinusOneToPower(m + repN)) *
               std::conj(stored(-m));
      } else {
        return sign * turn * stored(m);
      }
    }
  }

 private:
  LayeredDetails::Stacks<LayeredDetails::ExpansionStackAt, Flat> _stacks;
  RadialGridType _radialGrid;
  GridType _grid;
  Int _lMax;

  template <Int... Alphas>
  static constexpr Int SlotOf() {
    return Flat::SlotOfFlat(Flat::template FlatOf<Alphas...>);
  }

  template <std::size_t... Slots>
  static auto Build(const RadialGridType& radialGrid, const GridType& grid,
                    Int lMax, std::index_sequence<Slots...>) {
    return LayeredDetails::Stacks<LayeredDetails::ExpansionStackAt, Flat>{
        LayeredDetails::ExpansionStackAt<Flat, static_cast<Int>(Slots)>(
            radialGrid, grid, lMax)...};
  }
};

//--------------------------------------------------------------------------//
//                          Between the two domains                          //
//--------------------------------------------------------------------------//

namespace LayeredDetails {

// Expand a compile-time multi-index into the pack the component accessors
// want. The same trick ContravariantDerivative.h uses, and for the same
// reason: the traversal over slots is over `constexpr` arrays, and the
// accessors take template parameter packs.
template <auto Indices, typename Result, typename Field, std::size_t... I>
void ExpandOne(Result& result, const Field& field, Int lMax, Execution policy,
               std::index_sequence<I...>) {
  const auto& in = field.template ComponentStack<Indices[I]...>();
  auto& out = result.template ComponentStack<Indices[I]...>();
  constexpr auto N = std::remove_cvref_t<decltype(in)>::UpperIndex;
  auto target = out.Data();
  in.Grid().ForwardTransformation(lMax, N, in.Data(), in.Batch(), target,
                                  out.Batch(), policy);
}

template <auto Indices, typename Result, typename Expansion, std::size_t... I>
void EvaluateOne(Result& result, const Expansion& expansion, Execution policy,
                 std::index_sequence<I...>) {
  const auto& in = expansion.template ComponentStack<Indices[I]...>();
  auto& out = result.template ComponentStack<Indices[I]...>();
  constexpr auto N = std::remove_cvref_t<decltype(in)>::UpperIndex;
  auto target = out.Data();
  in.Grid().InverseTransformation(in.MaxDegree(), N, in.Data(), in.Batch(),
                                  target, out.Batch(), policy);
}

}  // namespace LayeredDetails

// One batched transform per stored component: nR fields at a time, rather than
// one field at a time as a loop over radii and components would give.
template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid, SlotAlphabet SlotSet>
auto Expand(
    const LayeredTensorField<Rank, Symmetry, Reality, Grid, SlotSet>& field,
    std::ptrdiff_t lMax, Execution policy = Execution::Sequential()) {
  using Result = LayeredTensorExpansion<Rank, Symmetry, Reality, Grid, SlotSet>;
  using Flat = typename Result::Flat;

  auto result = Result(field.Radial(), field.Grid(), lMax);
  [&]<std::size_t... Slots>(std::index_sequence<Slots...>) {
    (
        [&] {
          constexpr auto flat = Flat::ComponentLayout.flatOfSlot[Slots];
          constexpr auto indices =
              MultiIndex<Rank, SlotSet>::FromFlat(flat).Slots();
          LayeredDetails::ExpandOne<indices>(
              result, field, lMax, policy,
              std::make_index_sequence<static_cast<std::size_t>(Rank)>{});
        }(),
        ...);
  }(std::make_index_sequence<static_cast<std::size_t>(
        Result::StoredComponents)>{});
  return result;
}

template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid, SlotAlphabet SlotSet>
auto Evaluate(const LayeredTensorExpansion<Rank, Symmetry, Reality, Grid,
                                           SlotSet>& expansion,
              Execution policy = Execution::Sequential()) {
  using Result = LayeredTensorField<Rank, Symmetry, Reality, Grid, SlotSet>;
  using Flat = typename Result::Flat;

  auto result = Result(expansion.Radial(), expansion.Grid());
  [&]<std::size_t... Slots>(std::index_sequence<Slots...>) {
    (
        [&] {
          constexpr auto flat = Flat::ComponentLayout.flatOfSlot[Slots];
          constexpr auto indices =
              MultiIndex<Rank, SlotSet>::FromFlat(flat).Slots();
          LayeredDetails::EvaluateOne<indices>(
              result, expansion, policy,
              std::make_index_sequence<static_cast<std::size_t>(Rank)>{});
        }(),
        ...);
  }(std::make_index_sequence<static_cast<std::size_t>(
        Result::StoredComponents)>{});
  return result;
}

//--------------------------------------------------------------------------//
//                             The usual ranks                               //
//--------------------------------------------------------------------------//

// Names for what the applications actually hold. A rank-0 tensor is a scalar
// and a rank-1 one is a vector; spelling them out at every use adds a symmetry
// argument that has only one possible value and a rank that the name already
// says.
template <AngularGrid Grid, TensorReality Reality = RealTensor>
using LayeredScalarField = LayeredTensorField<0, NoSymmetry<0>, Reality, Grid>;

template <AngularGrid Grid, TensorReality Reality = RealTensor>
using LayeredScalarExpansion =
    LayeredTensorExpansion<0, NoSymmetry<0>, Reality, Grid>;

template <AngularGrid Grid, TensorReality Reality = RealTensor>
using LayeredVectorField = LayeredTensorField<1, NoSymmetry<1>, Reality, Grid>;

template <AngularGrid Grid, TensorReality Reality = RealTensor>
using LayeredVectorExpansion =
    LayeredTensorExpansion<1, NoSymmetry<1>, Reality, Grid>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_LAYERED_TENSOR_FIELD_GUARD_H
