#ifndef GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
#define GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H

/**
 * @file GaussLegendreGrid.h
 * @brief A Gauss-Legendre grid: the quadrature, and nothing else.
 */

#include <FFTWpp/Core>
#include <GaussQuad/All>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Concepts.h"
#include "Policies.h"
#include "SphericalGrid.h"

namespace GSHTrans {

/**
 * @brief A Gauss-Legendre grid: the quadrature, and nothing else.
 *
 * @details Everything a transform needs lives in SphericalGrid, which holds
 * the Wigner values, both Legendre kernels, the Fourier stages, batching,
 * chunking, threading and the plan cache. This class supplies the nodes and
 * the weights and adds no data of its own, so slicing a grid to its base
 * loses ForBand() and nothing more.
 *
 * The nodes are the Gauss-Legendre points of degree @f$l_{\max} + 1@f$ mapped
 * by @f$\theta = \arccos(-x)@f$, which puts @f$l_{\max} + 1@f$ colatitudes
 * strictly inside @f$(0, \pi)@f$, symmetric about the equator. Both facts are
 * load-bearing: the interior condition is what the transform's contract asks
 * for and what the interpolant's polar padding rests on, and the symmetry is
 * what lets the matrix kernel store half the Wigner values.
 *
 * Exactness is the reason to prefer it. Gauss-Legendre integrates a degree
 * @f$2n-1@f$ polynomial exactly with @f$n@f$ nodes, so @f$l_{\max} + 1@f$
 * colatitudes resolve a band-@f$l_{\max}@f$ field with no quadrature error at
 * all, on half the latitudes an equiangular rule would need for the same
 * band.
 *
 * @tparam _Real The precision.
 * @tparam _MRange Whether all orders are stored, or only the non-negative
 * ones.
 * @tparam _NRange Which upper indices the grid covers.
 */
template <RealFloatingPoint _Real, OrderIndexRange _MRange, IndexRange _NRange>
class GaussLegendreGrid : public SphericalGrid<_Real, _MRange, _NRange> {
  /// The base this derives from.
  using Base = SphericalGrid<_Real, _MRange, _NRange>;

 public:
  using Int = std::ptrdiff_t;          ///< Signed index type used throughout.
  using Real = _Real;                  ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  /// Whether all orders are stored, or only the non-negative ones.
  using MRange = _MRange;
  using NRange = _NRange;  ///< Which upper indices are covered.

  GaussLegendreGrid() = delete;

  /**
   * @brief A grid resolving degrees up to @p lMax.
   *
   * @details The chunking policy is a property of the machine rather than of
   * the call, which is why it is set here alongside the planner flag and not
   * on every transform.
   *
   * @param lMax The largest degree the grid resolves, giving
   * @f$l_{\max} + 1@f$ colatitudes.
   * @param nMax The largest upper index the grid covers.
   * @param flag How hard FFTW should work at planning.
   * @param chunking How many fields of a batch the inner loop takes at once.
   * @param values Whether the Wigner table is stored or generated.
   * @param kernel Which Legendre kernel the transforms use.
   */
  GaussLegendreGrid(Int lMax, Int nMax, FFTWpp::Flag flag = FFTWpp::Measure,
                    Chunking chunking = Chunking::Automatic(),
                    WignerValues values = WignerValues::Stored(),
                    TransformKernel kernel = TransformKernel::Loop())
      : GaussLegendreGrid{Nodes(lMax), lMax,   nMax,  flag,
                          chunking,    values, kernel} {}

  /**
   * @brief A grid for fields of band @p lBand, with quadrature headroom.
   *
   * @details The distinction this expresses is between the *band* of the
   * fields being worked with and the *resolution* of the grid. Headroom is
   * wanted whenever a quantity of higher degree than the fields themselves is
   * formed on the grid: a product of two band-@f$L@f$ fields has band
   * @f$2L@f$, and integrating @f$|f|^2@f$ for band-@f$L@f$ @f$f@f$ integrates
   * a degree-@f$2L@f$ quantity. The 3/2 rule is `oversampling = 1.5`;
   * `oversampling = 2` is exact for a single product; a local interpolation
   * scheme wants rather more, four to eight.
   *
   * @param lBand The band of the fields to be worked with.
   * @param nMax The largest upper index the grid covers.
   * @param oversampling The headroom factor, at least one; the grid resolves
   * degrees up to @f$\lceil \mathrm{oversampling} \times l_{band} \rceil@f$.
   * @param flag How hard FFTW should work at planning.
   * @param chunking How many fields of a batch the inner loop takes at once.
   * @param values Whether the Wigner table is stored or generated.
   * @param kernel Which Legendre kernel the transforms use.
   */
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
    const auto lGrid =
        static_cast<Int>(std::ceil(oversampling * static_cast<Real>(lBand)));
    return GaussLegendreGrid(lGrid, nMax, flag, chunking, values, kernel);
  }

  /** @brief Copy constructor. */
  GaussLegendreGrid(const GaussLegendreGrid&) = default;
  /** @brief Move constructor. */
  GaussLegendreGrid(GaussLegendreGrid&&) = default;
  /** @brief Copy assignment. */
  GaussLegendreGrid& operator=(const GaussLegendreGrid&) = default;
  /** @brief Move assignment. */
  GaussLegendreGrid& operator=(GaussLegendreGrid&&) = default;

  /**
   * @brief The same grid, sharing one table, under a different chunking
   * policy.
   * @details The base's `With()` returns a base, which would lose ForBand().
   * This keeps the derived type, and costs nothing since the derived class
   * adds no data.
   * @param chunking The policy to use instead.
   */
  auto With(Chunking chunking) const {
    auto grid = *this;
    static_cast<Base&>(grid) = Base::With(chunking);
    return grid;
  }

  /**
   * @brief The same grid, sharing one table, under a different planner flag.
   * @param flag The flag to plan with instead.
   */
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
