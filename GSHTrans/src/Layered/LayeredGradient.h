#ifndef GSH_TRANS_LAYERED_GRADIENT_GUARD_H
#define GSH_TRANS_LAYERED_GRADIENT_GUARD_H

#include <complex>
#include <cstddef>
#include <stdexcept>
#include <utility>

#include "../Concepts.h"
#include "../Policies.h"
#include "../Expansion/ContravariantDerivative.h"
#include "../Tensor/MultiIndex.h"
#include "LayeredTensorField.h"
#include "RadialOperator.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                    The gradient of a three-dimensional field              //
//--------------------------------------------------------------------------//

// D&T's gradient in the canonical basis (C.151)-(C.153) is
//
//     grad = e_0 d_r + r^{-1} grad_1,
//
// and the two halves live on opposite sides of the seam. The angular half,
// grad_1, is `SurfaceGradient` -- already written, already tested, and the same
// formula whether or not there is a radial axis. The radial half is d_r, which
// this library does not and will not own: it belongs to the application's
// discretisation and arrives as a RadialOperator.
//
// So the full gradient of a rank-q field is a rank-(q+1) field whose blocks
// split by the *new* leading index sigma:
//
//     sigma = +-1 :  r^{-1} (grad_1 T)^{sigma a_1...a_q}
//     sigma =   0 :  d_r T^{a_1...a_q}
//
// The e_0 basis vector is constant along r and the e_{+-} depend only on
// (theta, phi), so the radial derivative carries no connection terms and the
// sigma = 0 block really is the componentwise radial derivative. All the
// connection terms are in grad_1, including the ones that move a slot between
// e_0 and e_{+-}, which is why grad_1 has to be applied to a tensor whose slots
// already run over all three canonical directions and not only the tangential
// two.
//
// Everything happens in the spectral domain, which is where the model
// application does it and where it is cheapest: grad_1 is a multiplication by
// Omega and a subtraction at fixed (l, m), the radial operator commutes with
// the angular transform, and a coefficient buffer is complex whatever the
// field's reality, so a radial operator here only ever sees complex data.

namespace LayeredGradientDetails {

using Int = std::ptrdiff_t;

// The operand seen at one radius: exactly the accessor the gradient formula
// wants, with the radius bound. Nothing is copied.
template <typename Layered>
class OperandAtRadius {
 public:
  OperandAtRadius(const Layered& layered, Int i) : _layered{layered}, _i{i} {}

  template <Int... Alphas>
  auto Coefficient(Int l, Int m) const {
    return _layered.template Coefficient<Alphas...>(_i, l, m);
  }

 private:
  const Layered& _layered;
  Int _i;
};

// One result block at one radius, presented as the gradient formula expects a
// block to look: degrees, orders, and a writable coefficient.
template <typename Stack>
class BlockAtRadius {
 public:
  BlockAtRadius(Stack& stack, Int i) : _stack{stack}, _i{i} {}

  auto Degrees() const { return _stack.Degrees(); }
  auto Orders(Int l) const { return _stack.Orders(l); }
  auto& operator[](Int l, Int m) { return _stack[_i, l, m]; }

 private:
  Stack& _stack;
  Int _i;
};

template <typename Layered>
class ResultAtRadius {
 public:
  using Real = typename Layered::Real;
  using Complex = typename Layered::Complex;
  static constexpr auto& Orbits = Layered::Orbits;

  ResultAtRadius(Layered& layered, Int i) : _layered{layered}, _i{i} {}

  template <Int... Alphas>
  auto Component() {
    return BlockAtRadius(_layered.template ComponentStack<Alphas...>(), _i);
  }

 private:
  Layered& _layered;
  Int _i;
};

// Apply a radial operator to every stored component of a layered tensor
// expansion. The derived components follow from the stored ones by relations
// with constant coefficients, so differentiating the representatives and
// deriving is the same as deriving and differentiating.
template <auto Indices, typename Stacks, typename Op, std::size_t... I>
void ApplyOne(const Stacks& in, Stacks& out, const Op& op, Execution policy,
              std::index_sequence<I...>) {
  ApplyRadially(in.template ComponentStack<Indices[I]...>(),
                out.template ComponentStack<Indices[I]...>(), op, policy);
}

// The sigma = 0 block of the result: the radial derivative, read through the
// general accessor so that an operand component which is derived rather than
// stored costs nothing extra here.
template <auto Indices, typename Result, typename Derivative>
void FillRadialComponent(Result& result, const Derivative& derivative, Int i) {
  using Complex = typename Result::Complex;

  constexpr auto Rank = static_cast<Int>(Indices.size()) - 1;
  constexpr auto source = ContravariantDetails::DropFirst<Rank>(Indices);
  constexpr auto flat = MultiIndex<Rank + 1>(Indices).Flat();
  constexpr auto constraint = Result::Orbits.constraint[flat];

  auto block = [&]<std::size_t... J>(std::index_sequence<J...>) {
    return BlockAtRadius(result.template ComponentStack<Indices[J]...>(), i);
  }(std::make_index_sequence<static_cast<std::size_t>(Rank + 1)>{});

  for (auto l : block.Degrees()) {
    for (auto m : block.Orders(l)) {
      const auto value = [&]<std::size_t... J>(std::index_sequence<J...>) {
        return derivative.template Coefficient<source[J]...>(i, l, m);
      }(std::make_index_sequence<static_cast<std::size_t>(Rank)>{});

      // As in the surface gradient, a pinned-imaginary component would store
      // the real field whose i-multiple it is. The branch is unreachable while
      // the result has NoSymmetry -- under negation alone the only self-paired
      // multi-index is the all-zero one, which is pinned real -- and is here so
      // that it stays correct if that changes.
      if constexpr (constraint == ComponentConstraint::Imaginary) {
        block[l, m] = Complex{0, -1} * value;
      } else {
        block[l, m] = value;
      }
    }
  }
}

}  // namespace LayeredGradientDetails

//--------------------------------------------------------------------------//

// The tangential part of the gradient: r^{-1} grad_1, with the sigma = 0 block
// present and zero.
//
// This is the surface gradient of the flat algebra with the r^{-1} reinstated,
// which is what makes it a gradient on the ball rather than on the unit
// sphere. It is complete on its own for a caller who wants only the tangential
// part, and it is the angular half of `Gradient`.
//
// Like its flat counterpart it takes a general tensor and not a tangential
// one: grad_1 moves slots between e_0 and e_+-, so it does not close on the
// tangential bundle, and a tangential operand is embedded first
// (field-algebra-plan.md section 18.2 [D9]). Every MultiIndex<Rank + 1> below
// is therefore over AllSlots by right rather than by oversight -- the result
// lands in the general bundle whatever the operand was.
template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid>
auto SurfaceGradient(
    const LayeredTensorExpansion<Rank, Symmetry, Reality, Grid>& operand) {
  using Result =
      LayeredTensorExpansion<Rank + 1, NoSymmetry<Rank + 1>, Reality, Grid>;
  using Real = typename Result::Real;
  using Flat = typename Result::Flat;

  const auto& radial = operand.Radial();
  if (radial.Radius(0) <= 0) {
    throw std::invalid_argument(
        "The gradient in canonical components carries an explicit r^{-1}, so "
        "it is not defined at r = 0; the origin is a coordinate singularity "
        "of the basis and not of the field");
  }

  auto result = Result(radial, operand.Grid(), operand.MaxDegree());

  for (auto i : radial.RadiusIndices()) {
    const auto scale = static_cast<Real>(1) / radial.Radius(i);
    auto at = LayeredGradientDetails::ResultAtRadius(result, i);
    const auto from = LayeredGradientDetails::OperandAtRadius(operand, i);

    [&]<std::size_t... Slots>(std::index_sequence<Slots...>) {
      (
          [&] {
            constexpr auto flat = Flat::ComponentLayout.flatOfSlot[Slots];
            constexpr auto indices =
                MultiIndex<Rank + 1>::FromFlat(flat).Slots();
            ContravariantDetails::FillComponent<indices>(at, from, scale);
          }(),
          ...);
    }(std::make_index_sequence<
        static_cast<std::size_t>(Result::StoredComponents)>{});
  }
  return result;
}

// The full gradient: rank q in, rank q + 1 out, with the radial derivative
// supplied by the caller.
//
// The operator is applied to the coefficients, which is legitimate because it
// acts along r alone and so commutes with the angular transform -- the property
// TestLayered pins down rather than assumes. It sees complex spans of length
// nR, one per (component, degree, order).
template <std::ptrdiff_t Rank, TensorSymmetry<Rank> Symmetry,
          TensorReality Reality, AngularGrid Grid, typename Op>
auto Gradient(const LayeredTensorExpansion<Rank, Symmetry, Reality, Grid>&
                  operand,
              const Op& radialOperator,
              Execution policy = Execution::Sequential()) {
  using Operand = LayeredTensorExpansion<Rank, Symmetry, Reality, Grid>;
  using OperandFlat = typename Operand::Flat;

  auto result = SurfaceGradient(operand);
  using Result = decltype(result);
  using ResultFlat = typename Result::Flat;

  // d_r T, componentwise on the operand's own stored set.
  auto derivative = Operand(operand.Radial(), operand.Grid(),
                            operand.MaxDegree());
  [&]<std::size_t... Slots>(std::index_sequence<Slots...>) {
    (
        [&] {
          constexpr auto flat = OperandFlat::ComponentLayout.flatOfSlot[Slots];
          constexpr auto indices = MultiIndex<Rank>::FromFlat(flat).Slots();
          LayeredGradientDetails::ApplyOne<indices>(
              operand, derivative, radialOperator, policy,
              std::make_index_sequence<static_cast<std::size_t>(Rank)>{});
        }(),
        ...);
  }(std::make_index_sequence<
      static_cast<std::size_t>(Operand::StoredComponents)>{});

  // and into the sigma = 0 block of the result.
  for (auto i : operand.Radial().RadiusIndices()) {
    [&]<std::size_t... Slots>(std::index_sequence<Slots...>) {
      (
          [&] {
            constexpr auto flat = ResultFlat::ComponentLayout.flatOfSlot[Slots];
            constexpr auto indices =
                MultiIndex<Rank + 1>::FromFlat(flat).Slots();
            if constexpr (indices[0] == 0) {
              LayeredGradientDetails::FillRadialComponent<indices>(
                  result, derivative, i);
            }
          }(),
          ...);
    }(std::make_index_sequence<
        static_cast<std::size_t>(Result::StoredComponents)>{});
  }
  return result;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_LAYERED_GRADIENT_GUARD_H
