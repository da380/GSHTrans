#ifndef GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
#define GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H

#include <FFTWpp/Core>
#include <GaussQuad/All>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Concepts.h"
#include "Policies.h"
#include "SphericalGrid.h"

namespace GSHTrans {

// A Gauss-Legendre grid: the quadrature, and nothing else.
//
// Everything a transform needs lives in SphericalGrid, which holds the Wigner
// values, both Legendre kernels, the Fourier stages, batching, chunking,
// threading and the plan cache. This class supplies the nodes and the weights
// and adds no data of its own -- so a copy of one into its base loses ForBand
// and nothing more ([C26]).
//
// The nodes are the Gauss-Legendre points of degree lMax + 1 mapped by
// theta = acos(-x), which puts nTheta = lMax + 1 colatitudes strictly inside
// (0, pi), symmetric about the equator. Both facts are load-bearing: the
// interior condition is what the transform's contract asks for and what
// Interpolate's polar padding rests on, and the symmetry is what lets the
// matrix kernel store half the Wigner values.
//
// Exactness is the reason to prefer it. Gauss-Legendre integrates a
// degree-(2n-1) polynomial exactly with n nodes, so lMax + 1 colatitudes
// resolve a band-lMax field with no quadrature error at all, on half the
// latitudes an equiangular rule would need for the same band.
template <RealFloatingPoint _Real, OrderIndexRange _MRange, IndexRange _NRange>
class GaussLegendreGrid : public SphericalGrid<_Real, _MRange, _NRange> {
  using Base = SphericalGrid<_Real, _MRange, _NRange>;

 public:
  using Int = std::ptrdiff_t;
  using Real = _Real;
  using Complex = std::complex<Real>;
  using MRange = _MRange;
  using NRange = _NRange;

  GaussLegendreGrid() = delete;

  // The chunking policy is a property of the machine rather than of the call,
  // which is why it is set here alongside the planner flag and not on every
  // transform. `Automatic` assumes a modest cache; a caller who knows their
  // machine passes `Chunking::ForCache(bytes)`, and one who has measured
  // their own optimum passes `Chunking::Fixed(count)`.
  GaussLegendreGrid(Int lMax, Int nMax, FFTWpp::Flag flag = FFTWpp::Measure,
                    Chunking chunking = Chunking::Automatic(),
                    WignerValues values = WignerValues::Stored(),
                    TransformKernel kernel = TransformKernel::Loop())
      : GaussLegendreGrid{Nodes(lMax), lMax,   nMax,  flag,
                          chunking,    values, kernel} {}

  // A grid for working with fields of maximum degree lBand, with quadrature
  // headroom for degree oversampling * lBand.
  //
  // The distinction this expresses is between the *band* of the fields being
  // worked with and the *resolution* of the grid, which the library could not
  // previously say: the transforms already take lMax per call, so an
  // oversampled grid with truncated transforms works, but there was no way to
  // ask for one. Headroom is wanted whenever a quantity of higher degree than
  // the fields themselves is formed on the grid -- a product of two band-L
  // fields has band 2L, and integrating |f|^2 for band-L f integrates a
  // degree-2L quantity. The 3/2 rule is oversampling = 1.5; oversampling = 2
  // is exact for a single product. A local interpolation scheme wants rather
  // more, four to eight (field-algebra-plan.md section 22).
  static auto ForBand(Int lBand, Int nMax, Real oversampling = 1,
                      FFTWpp::Flag flag = FFTWpp::Measure,
                      Chunking chunking = Chunking::Automatic(),
                      WignerValues values = WignerValues::Stored(),
                      TransformKernel kernel = TransformKernel::Loop()) {
    if (lBand < 0) {
      throw std::invalid_argument("Band must be non-negative");
    }
    if (!(oversampling >= 1)) {
      throw std::invalid_argument("Oversampling factor must be at least one");
    }
    const auto lGrid = static_cast<Int>(
        std::ceil(oversampling * static_cast<Real>(lBand)));
    return GaussLegendreGrid(lGrid, nMax, flag, chunking, values, kernel);
  }

  GaussLegendreGrid(const GaussLegendreGrid&) = default;
  GaussLegendreGrid(GaussLegendreGrid&&) = default;
  GaussLegendreGrid& operator=(const GaussLegendreGrid&) = default;
  GaussLegendreGrid& operator=(GaussLegendreGrid&&) = default;

  // The base's With() returns a base, which would lose ForBand. These keep the
  // derived type, and cost nothing since the derived class adds no data.
  auto With(Chunking chunking) const {
    auto grid = *this;
    static_cast<Base&>(grid) = Base::With(chunking);
    return grid;
  }

  auto With(FFTWpp::Flag flag) const {
    auto grid = *this;
    static_cast<Base&>(grid) = Base::With(flag);
    return grid;
  }

 private:
  // Delegated to, so that Nodes() is evaluated once rather than once per
  // argument. The nodes come first because a member initialiser cannot name
  // a later parameter.
  GaussLegendreGrid(std::pair<std::vector<Real>, std::vector<Real>> nodes,
                    Int lMax, Int nMax, FFTWpp::Flag flag, Chunking chunking,
                    WignerValues values, TransformKernel kernel)
      : Base{lMax,
             nMax,
             std::move(nodes.first),
             std::move(nodes.second),
             flag,
             chunking,
             values,
             kernel} {}

  // The Gauss-Legendre nodes and weights, as colatitudes.
  //
  // GaussQuad gives the points on [-1, 1]; theta = acos(-x) maps them to
  // (0, pi) in increasing order, and the weights are unchanged because the
  // transform's Legendre stage sums against dtheta-free weights -- the
  // Jacobian is already in the definition of the quadrature on the cosine.
  static std::pair<std::vector<Real>, std::vector<Real>> Nodes(Int lMax) {
    // GaussQuadrature returns the nodes and weights as a pair; Quadrature1D
    // is what carries the Transform below, so it is built explicitly rather
    // than relied on to convert.
    auto quadrature = GaussQuad::Quadrature1D<Real>(
        GaussQuad::LegendrePolynomial<Real>{}.GaussQuadrature(lMax + 1));
    quadrature.Transform([](auto x) { return std::acos(-x); },
                         [](auto) -> Real { return 1; });
    auto coLatitudes = std::vector<Real>();
    auto weights = std::vector<Real>();
    coLatitudes.reserve(static_cast<std::size_t>(lMax + 1));
    weights.reserve(static_cast<std::size_t>(lMax + 1));
    for (auto theta : quadrature.Points()) coLatitudes.push_back(theta);
    for (auto weight : quadrature.Weights()) weights.push_back(weight);
    return {std::move(coLatitudes), std::move(weights)};
  }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
