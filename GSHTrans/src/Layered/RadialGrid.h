#ifndef GSH_TRANS_RADIAL_GRID_GUARD_H
#define GSH_TRANS_RADIAL_GRID_GUARD_H

#include <algorithm>
#include <cstddef>
#include <memory>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../Concepts.h"

namespace GSHTrans {

// The radial half of a three-dimensional field: a set of radii, and weights
// for integrating over them.
//
// Deliberately thin. It carries nodes, weights and identity and nothing else:
// element connectivity, the spectral-element basis, differentiation matrices
// and any factorisation belong to the application that built them, and the
// library only ever needs to know how many radii there are, what they are, and
// how to integrate over them.
//
// What the radii are actually used for is the quadrature over the ball and the
// r^{-1} in the gradient. What the *count* is used for is the stack, and that
// is most of it.
//
// A value-semantic handle over shared immutable state, for the same reason
// GaussLegendreGrid is one: two layered fields are on the same radial grid when
// their handles agree, and copying one is a pointer copy.
template <RealFloatingPoint _Real>
class RadialGrid {
 public:
  using Int = std::ptrdiff_t;
  using Real = _Real;

  RadialGrid() = delete;

  // Nodes in increasing order, with the weights of whatever rule the caller
  // means to integrate with. The weights may be empty, in which case the grid
  // carries points alone and Integrate is unavailable.
  RadialGrid(std::vector<Real> radii, std::vector<Real> weights = {})
      : _impl{std::make_shared<const Impl>(std::move(radii),
                                           std::move(weights))} {}

  auto NumberOfRadii() const { return static_cast<Int>(_impl->radii.size()); }
  auto RadiusIndices() const {
    return std::ranges::views::iota(Int{0}, NumberOfRadii());
  }

  auto Radii() const { return std::span<const Real>(_impl->radii); }
  auto Weights() const { return std::span<const Real>(_impl->weights); }
  auto HasWeights() const { return !_impl->weights.empty(); }

  Real Radius(Int i) const { return _impl->radii[i]; }
  Real Weight(Int i) const { return _impl->weights[i]; }

  // Identity, not structure: two grids with equal radii built separately are
  // different grids, exactly as for the angular grid.
  auto Identity() const { return _impl.get(); }

 private:
  struct Impl {
    Impl(std::vector<Real> r, std::vector<Real> w)
        : radii{std::move(r)}, weights{std::move(w)} {
      if (radii.empty()) {
        throw std::invalid_argument("A radial grid needs at least one radius");
      }
      if (!weights.empty() && weights.size() != radii.size()) {
        throw std::invalid_argument(
            "A radial grid has " + std::to_string(radii.size()) +
            " radii but " + std::to_string(weights.size()) + " weights");
      }
      if (!std::ranges::is_sorted(radii)) {
        throw std::invalid_argument("Radii must be in increasing order");
      }
      if (radii.front() < 0) {
        throw std::invalid_argument("Radii must be non-negative");
      }
    }

    std::vector<Real> radii;
    std::vector<Real> weights;
  };

  std::shared_ptr<const Impl> _impl;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_RADIAL_GRID_GUARD_H
