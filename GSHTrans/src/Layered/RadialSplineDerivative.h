#ifndef GSH_TRANS_RADIAL_SPLINE_DERIVATIVE_GUARD_H
#define GSH_TRANS_RADIAL_SPLINE_DERIVATIVE_GUARD_H

// d/dr of the cubic spline through a radial line.
//
// The whole of this header is conditional on GSHTRANS_WITH_INTERPOLATION,
// which is on by default. Without it a caller has FiniteDifferenceDerivative
// and LagrangeDerivative, which are the two that need nothing outside the
// standard library.
//
// **This used to be self-contained and is not any more, and that is a
// deliberate reversal** (field-algebra-plan.md section 21.2 [R5]). It held its
// own natural-spline system and its own Thomas sweep, for one reason: the
// factorisation inside Interpolation was not reachable, so building on that
// library would have meant a spline constructed per radial line. It is
// reachable now -- `CubicSplineSystem` is the matrix, factorised once, with a
// `Solve` that allocates nothing and is `const`, which is exactly the contract
// RadialOperator.h asks an operator to meet. So sixty lines of duplicated
// spline go, and what is left is assembly.
//
// The exchange is worth more than the lines. This gains the endpoint
// conditions the hand-written version never had -- and NotAKnot in particular
// stays fourth order right up to the ends, where Natural costs an order, which
// is precisely where a spline derivative was least trustworthy and precisely
// where a boundary condition gets applied.

#ifdef GSHTRANS_HAVE_INTERPOLATION

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Interpolation/CubicSplineSystem.hpp>

#include "../Concepts.h"
#include "RadialGrid.h"

namespace GSHTrans {

// The endpoint conditions, taken from Interpolation rather than mirrored, so
// that there is one vocabulary for them and not two.
using Interpolation::BoundaryCondition;

// The middle ground between the two operators of RadialDerivatives.h: global
// like the differentiation matrix, but linear in the number of radii rather
// than quadratic, and unlike a global polynomial it does not fall apart as the
// nodes multiply.
//
// A convenience rather than a workhorse. The discretisations these codes
// actually run on are finite differences, a finite-element basis, or a radial
// spectral basis; a spline is what a caller reaches for to process a field
// rather than to solve on one. Which is why this is written for clarity and
// correctness and not tuned.
template <RealFloatingPoint _Real>
class SplineDerivative {
 public:
  using Int = std::ptrdiff_t;
  using Real = _Real;
  using System = Interpolation::CubicSplineSystem<std::span<const Real>>;

  SplineDerivative() = delete;

  explicit SplineDerivative(
      RadialGrid<Real> radial,
      BoundaryCondition left = BoundaryCondition::Natural,
      BoundaryCondition right = BoundaryCondition::Natural)
      : _radial{std::move(radial)} {
    if (left == BoundaryCondition::Clamped ||
        right == BoundaryCondition::Clamped) {
      throw std::invalid_argument(
          "A clamped spline needs the end slopes of the data it is fitted to, "
          "and a radial operator is handed one line at a time with no way to "
          "say what they are -- so this operator offers Natural and NotAKnot "
          "and a caller wanting clamped ends supplies their own operator");
    }

    // One system per element where the grid knows them, and one over the whole
    // grid where it does not. That is the whole of what the partition changes
    // here: a spline never spans an interface, so a layered model is an
    // ordinary case rather than a refusal.
    const auto radii = _radial.Radii();
    if (_radial.HasElements()) {
      for (auto k : _radial.ElementIndices()) {
        _block.emplace_back(_radial.ElementStart(k), _radial.ElementSize(k));
      }
    } else {
      for (std::size_t i = 0; i + 1 < radii.size(); i++) {
        if (!(radii[i] < radii[i + 1])) {
          throw std::invalid_argument(
              "A spline derivative needs strictly increasing radii, and this "
              "grid repeats one without saying what it means -- which is what "
              "RadialGrid::WithElements is for, and with a partition this "
              "operator fits each element separately instead of refusing");
        }
      }
      _block.emplace_back(Int{0}, static_cast<Int>(radii.size()));
    }

    for (const auto& [first, count] : _block) {
      _systems.emplace_back(radii.subspan(static_cast<std::size_t>(first),
                                          static_cast<std::size_t>(count)),
                            left, right);
    }
  }

  const RadialGrid<Real>& Radial() const { return _radial; }

  // How many splines a line is carried by: one per element, or one.
  Int PieceCount() const { return static_cast<Int>(_systems.size()); }

  // The factorised system of one piece, for a caller who wants the curvatures
  // themselves rather than the slopes -- which is the point of the upstream
  // type being public, so it would be odd to hide it again here.
  const System& SplineSystem(Int piece = 0) const {
    return _systems[static_cast<std::size_t>(piece)];
  }

  template <typename Scalar>
  void operator()(std::span<const Scalar> in, std::span<Scalar> out) const {
    const auto n = static_cast<std::size_t>(_radial.NumberOfRadii());
    if (in.size() != n || out.size() != n) {
      throw std::invalid_argument(
          "A radial operator acts on a line of one value per radius");
    }

    // Scratch, not state: Solve is const and allocates nothing, so one
    // operator serves every thread and only the curvatures are per-call.
    thread_local auto curvature = std::vector<Scalar>{};
    if (curvature.size() < n) curvature.resize(n);

    for (std::size_t p = 0; p < _systems.size(); p++) {
      const auto [first, count] = _block[p];
      const auto lo = static_cast<std::size_t>(first);
      const auto size = static_cast<std::size_t>(count);
      const auto work = std::span<Scalar>(curvature.data(), size);

      _systems[p].Solve(in.subspan(lo, size), work);
      _systems[p].template EvaluateAtNodes<1>(
          in.subspan(lo, size), std::span<const Scalar>(curvature.data(), size),
          out.subspan(lo, size));
    }
  }

 private:
  RadialGrid<Real> _radial;
  std::vector<std::pair<Int, Int>> _block;
  std::vector<System> _systems;
};

}  // namespace GSHTrans

#endif  // GSHTRANS_HAVE_INTERPOLATION

#endif  // GSH_TRANS_RADIAL_SPLINE_DERIVATIVE_GUARD_H
