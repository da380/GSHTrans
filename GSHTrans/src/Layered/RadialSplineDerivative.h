#ifndef GSH_TRANS_RADIAL_SPLINE_DERIVATIVE_GUARD_H
#define GSH_TRANS_RADIAL_SPLINE_DERIVATIVE_GUARD_H

/**
 * @file RadialSplineDerivative.h
 * @brief @f$d/dr@f$ of the cubic spline through a radial line.
 *
 * @details The whole of this header is conditional on the interpolation
 * dependency, which is on by default. Without it a caller has
 * FiniteDifferenceDerivative and LagrangeDerivative, which are the two that
 * need nothing outside the standard library.
 *
 * **This is deliberately not self-contained.** The spline system comes from
 * that dependency rather than being written here: `CubicSplineSystem` is the
 * matrix, factorised once, with a `Solve` that allocates nothing and is
 * `const`, which is exactly the contract RadialOperator.h asks an operator to
 * meet. It also carries the endpoint conditions — and NotAKnot in particular
 * stays fourth order right up to the ends, where Natural costs an order,
 * which is precisely where a spline derivative is least trustworthy and
 * precisely where a boundary condition gets applied.
 */

#ifdef GSHTRANS_HAVE_INTERPOLATION

#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Interpolation/CubicSplineSystem.hpp>

#include "../Concepts.h"
#include "RadialGrid.h"

namespace GSHTrans {

// The endpoint conditions, taken from Interpolation rather than mirrored, so
// that there is one vocabulary for them and not two.
using Interpolation::BoundaryCondition;

/// The middle ground between the two operators of RadialDerivatives.h: global
/// like the differentiation matrix, but linear in the number of radii rather
/// than quadratic, and unlike a global polynomial it does not fall apart as the
/// nodes multiply.
///
/// A convenience rather than a workhorse. The discretisations these codes
/// actually run on are finite differences, a finite-element basis, or a radial
/// spectral basis; a spline is what a caller reaches for to process a field
/// rather than to solve on one. Which is why this is written for clarity and
/// correctness and not tuned.
template <RealFloatingPoint _Real>
class SplineDerivative {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  using Real = _Real;  ///< The precision.
  using System = Interpolation::CubicSplineSystem<
      std::span<const Real>>;  ///< The factorised spline system.

  SplineDerivative() = delete;

  /**
   * @brief Factorises one spline system per element.
   * @param radial The radii; if they know their elements there is one system
   * per element, otherwise one for the whole grid.
   * @param left The condition at the first node.
   * @param right The condition at the last.
   */
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

  /** @brief The radial grid this is defined on. */
  const RadialGrid<Real>& Radial() const { return _radial; }

  /// How many splines a line is carried by: one per element, or one.
  Int PieceCount() const { return static_cast<Int>(_systems.size()); }

  /// The factorised system of one piece, for a caller who wants the curvatures
  /// themselves rather than the slopes -- which is the point of the upstream
  /// type being public, so it would be odd to hide it again here.
  const System& SplineSystem(Int piece = 0) const {
    return _systems[static_cast<std::size_t>(piece)];
  }

  /**
   * @brief Differentiates one radial line.
   * @param in One value per radius.
   * @param out Where the derivative goes, of the same length.
   * @throws std::invalid_argument if either is not one value per radius.
   */
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
