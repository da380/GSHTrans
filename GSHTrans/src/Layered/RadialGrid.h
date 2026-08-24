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
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.
  using Real = _Real;  ///< The precision.

  RadialGrid() = delete;

  // Nodes in increasing order, with the weights of whatever rule the caller
  // means to integrate with. The weights may be empty, in which case the grid
  // carries points alone and Integrate is unavailable.
  RadialGrid(std::vector<Real> radii, std::vector<Real> weights = {})
      : _impl{std::make_shared<const Impl>(std::move(radii),
                                           std::move(weights),
                                           std::vector<Int>{})} {}

  // The same, knowing which radii belong to which element.
  //
  // This is the one piece of *structure* the grid carries, and it is here
  // because it is the fact more than one thing needs and nothing can infer:
  // a block-diagonal derivative needs the blocks, and an interpolant must not
  // cross a material interface. Without it a repeated radius is
  // indistinguishable from a caller's mistake -- the grid has always accepted
  // one, and has never been able to say what it meant.
  //
  // The blocks are disjoint and given by their `nElements + 1` starts, so
  // every node belongs to exactly one element and two elements meet at a
  // *repeated radius* (field-algebra-plan.md section 20.5 [E1]). That is what
  // makes a derivative well defined at an interface without anyone having to
  // choose between averaging the two sides and picking one: both exist, and
  // they are the pair `Interpolation::Piecewise::Limits` returns on the other
  // side of the join.
  //
  // The price, and it is the caller's: a continuous mesh duplicates its
  // interior element boundaries and keeps the duplicates equal. That is the
  // same contract `Piecewise` sets when it declines to check continuity, and
  // for the same reason -- enforcing it would only start an argument about
  // tolerance.
  //
  // A named constructor rather than a third defaulted parameter, so that a
  // caller with elements and no weights does not have to pass an empty vector
  // to reach it.
  static RadialGrid WithElements(std::vector<Real> radii,
                                 std::vector<Int> elementStarts,
                                 std::vector<Real> weights = {}) {
    return RadialGrid(std::make_shared<const Impl>(std::move(radii),
                                                   std::move(weights),
                                                   std::move(elementStarts)));
  }

  /** @brief How many radii the stack holds. */
  auto NumberOfRadii() const { return static_cast<Int>(_impl->radii.size()); }
  /** @brief Indices of the stored radii. */
  auto RadiusIndices() const {
    return std::ranges::views::iota(Int{0}, NumberOfRadii());
  }

  auto Radii() const { return std::span<const Real>(_impl->radii); }
  /** @brief The quadrature weight at every point. */
  auto Weights() const { return std::span<const Real>(_impl->weights); }
  auto HasWeights() const { return !_impl->weights.empty(); }

  Real Radius(Int i) const { return _impl->radii[i]; }
  Real Weight(Int i) const { return _impl->weights[i]; }

  //------------------------------------------------------------------------//
  //                             The elements                                //
  //------------------------------------------------------------------------//

  // Whether this grid knows its elements at all. A grid built without them is
  // exactly what it was before, and everything that reads them says what it
  // does when there are none.
  auto HasElements() const { return !_impl->starts.empty(); }

  auto ElementCount() const {
    return HasElements() ? static_cast<Int>(_impl->starts.size()) - 1 : Int{0};
  }

  auto ElementIndices() const {
    return std::ranges::views::iota(Int{0}, ElementCount());
  }

  // Element k holds the radii [ElementStart(k), ElementEnd(k)), which is the
  // half-open form everything else here uses.
  Int ElementStart(Int k) const { return _impl->starts[k]; }
  Int ElementEnd(Int k) const { return _impl->starts[k + 1]; }
  Int ElementSize(Int k) const { return ElementEnd(k) - ElementStart(k); }

  auto ElementRadiusIndices(Int k) const {
    return std::ranges::views::iota(ElementStart(k), ElementEnd(k));
  }

  // Which element a radius index belongs to. Exactly one does, the blocks
  // being disjoint, which is the whole point of [E1].
  Int ElementOf(Int i) const {
    for (auto k = Int{0}; k < ElementCount(); k++) {
      if (i < ElementEnd(k)) return k;
    }
    return ElementCount() - 1;
  }

  // The breakpoints, of which there are ElementCount() + 1. Read off the
  // radii rather than stored beside them, so that there is one source of
  // truth and no way for the two to disagree ([E3]).
  Real Breakpoint(Int k) const {
    return k == 0 ? _impl->radii.front()
                  : _impl->radii[static_cast<std::size_t>(ElementEnd(k - 1) - 1)];
  }

  // Identity, not structure: two grids with equal radii built separately are
  // different grids, exactly as for the angular grid.
  /** @brief Identity, for deciding whether two share an implementation. */
  auto Identity() const { return _impl.get(); }

 private:
  struct Impl {
    Impl(std::vector<Real> r, std::vector<Real> w, std::vector<Int> s)
        : radii{std::move(r)}, weights{std::move(w)}, starts{std::move(s)} {
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
      if (!starts.empty()) ValidateElements();
    }

    // Everything the partition has to satisfy for the blocks to be elements
    // rather than an arbitrary grouping. Each check is a different way of
    // being wrong, so each has its own message.
    void ValidateElements() const {
      const auto nR = static_cast<Int>(radii.size());
      if (starts.size() < 2) {
        throw std::invalid_argument(
            "An element partition needs at least one element, so at least two "
            "starts");
      }
      if (starts.front() != 0 || starts.back() != nR) {
        throw std::invalid_argument(
            "An element partition must cover every radius, so it starts at 0 "
            "and ends at " + std::to_string(nR));
      }
      for (std::size_t k = 0; k + 1 < starts.size(); k++) {
        // At least two nodes: a block spanning no interval is not an element.
        if (starts[k + 1] - starts[k] < 2) {
          throw std::invalid_argument(
              "Element " + std::to_string(k) +
              " holds fewer than two radii, so it spans no interval");
        }
        // Strictly increasing *within* a block, which is the statement that a
        // repeated radius means an interface and never anything else.
        for (auto i = starts[k]; i + 1 < starts[k + 1]; i++) {
          if (!(radii[static_cast<std::size_t>(i)] <
                radii[static_cast<std::size_t>(i + 1)])) {
            throw std::invalid_argument(
                "Radii repeat inside element " + std::to_string(k) +
                ", and a repeated radius is how two elements meet rather than "
                "something that happens within one");
          }
        }
      }
      // And the blocks tile: where one ends the next begins, at the same
      // radius. A gap would leave the field undefined between them.
      for (std::size_t k = 1; k + 1 < starts.size(); k++) {
        const auto last = static_cast<std::size_t>(starts[k] - 1);
        const auto first = static_cast<std::size_t>(starts[k]);
        if (radii[last] != radii[first]) {
          throw std::invalid_argument(
              "Elements " + std::to_string(k - 1) + " and " +
              std::to_string(k) +
              " do not meet: one ends at a different radius from the one the "
              "next begins at, which leaves a gap the field is not defined "
              "on");
        }
      }
    }

    std::vector<Real> radii;
    std::vector<Real> weights;
    std::vector<Int> starts;
  };

  explicit RadialGrid(std::shared_ptr<const Impl> impl)
      : _impl{std::move(impl)} {}

  std::shared_ptr<const Impl> _impl;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_RADIAL_GRID_GUARD_H
