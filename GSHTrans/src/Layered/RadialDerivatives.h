#ifndef GSH_TRANS_RADIAL_DERIVATIVES_GUARD_H
#define GSH_TRANS_RADIAL_DERIVATIVES_GUARD_H

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../Concepts.h"
#include "RadialGrid.h"

namespace GSHTrans {

// Ready-made radial derivatives, offered and not imposed.
//
// The seam of RadialOperator.h is where the library stops: it applies whatever
// callable it is given along the radial axis and owns no discretisation. That
// is still true. What these add is that `Gradient` can be used without the
// caller writing a differentiation matrix first, which until now they had to.
//
// They are also meant to be *read*. Each is a type satisfying RadialOperator
// and nothing more, so a caller whose discretisation is none of these has a
// worked example rather than a framework -- and in particular has an example
// of the two obligations RadialOperator.h states. Every one of these does all
// of its node-dependent work in the constructor and allocates nothing in the
// call, and every one keeps its per-line scratch in thread_local storage so
// that one operator can serve every thread.
//
// All of them are linear maps that depend on the radii alone, which is what
// makes that possible: what varies from line to line is the data, and what
// costs anything to compute does not.
//
// The two here need nothing outside the standard library. SplineDerivative is
// the third and lives in RadialSplineDerivative.h, because it is built on
// Interpolation's factorised spline system and so exists only when that
// dependency does.

namespace RadialDetails {

using Int = std::ptrdiff_t;

// Fornberg's algorithm (1988) for finite-difference weights on arbitrary
// spacing: the weights of the first derivative at `z` of a function sampled at
// `nodes`. Exact for polynomials of degree less than the stencil size.
//
// The recurrence builds the weights for every derivative order up to the one
// wanted while it walks the nodes, which is why the table is two columns wide
// and only the second is returned. Written out rather than taken from a
// library because it is fifteen lines and the alternative is a dependency for
// fifteen lines.
template <RealFloatingPoint Real>
auto FirstDerivativeWeights(Real z, std::span<const Real> nodes) {
  const auto n = static_cast<Int>(nodes.size());
  auto c = std::vector<Real>(static_cast<std::size_t>(2 * n), Real{0});
  const auto at = [n](Int i, Int k) {
    return static_cast<std::size_t>(i * 2 + k);
  };

  auto c1 = Real{1};
  auto c4 = nodes[0] - z;
  c[at(0, 0)] = 1;

  for (auto i = Int{1}; i < n; i++) {
    const auto mn = i < 1 ? i : Int{1};
    auto c2 = Real{1};
    const auto c5 = c4;
    c4 = nodes[static_cast<std::size_t>(i)] - z;

    for (auto j = Int{0}; j < i; j++) {
      const auto c3 = nodes[static_cast<std::size_t>(i)] -
                      nodes[static_cast<std::size_t>(j)];
      c2 *= c3;
      if (j == i - 1) {
        for (auto k = mn; k >= 1; k--) {
          c[at(i, k)] =
              c1 * (static_cast<Real>(k) * c[at(i - 1, k - 1)] -
                    c5 * c[at(i - 1, k)]) / c2;
        }
        c[at(i, 0)] = -c1 * c5 * c[at(i - 1, 0)] / c2;
      }
      for (auto k = mn; k >= 1; k--) {
        c[at(j, k)] =
            (c4 * c[at(j, k)] - static_cast<Real>(k) * c[at(j, k - 1)]) / c3;
      }
      c[at(j, 0)] = c4 * c[at(j, 0)] / c3;
    }
    c1 = c2;
  }

  auto weights = std::vector<Real>(static_cast<std::size_t>(n));
  for (auto i = Int{0}; i < n; i++) {
    weights[static_cast<std::size_t>(i)] = c[at(i, 1)];
  }
  return weights;
}

// The differentiation matrix of a set of nodes, row-major: D(j, i) is the
// derivative at node j of the cardinal polynomial through node i.
//
// Built barycentrically, which is the numerically sound way to form it, and
// shared between LagrangeDerivative -- which uses it over the whole grid --
// and ElementDerivative, which uses one per element. A spectral element *is*
// this matrix over its own nodes, so writing it twice would have been writing
// the same thing twice.
//
// The diagonal is minus the row sum: the negative-sum trick. It is the
// statement that a constant differentiates to zero, imposed rather than
// evaluated from a closed form -- so the error in differentiating a constant
// is the rounding of one sum rather than the accuracy of that form, which is
// what would otherwise grow with the number of nodes. Applying the matrix sums
// in a different order, so the answer is a few epsilon and not identically
// zero.
template <RealFloatingPoint Real>
auto DifferentiationMatrix(std::span<const Real> nodes) {
  const auto n = static_cast<Int>(nodes.size());

  // Barycentric weights, w_i = 1 / prod_{k != i} (r_i - r_k).
  auto w = std::vector<Real>(static_cast<std::size_t>(n), Real{1});
  for (auto i = Int{0}; i < n; i++) {
    for (auto k = Int{0}; k < n; k++) {
      if (k == i) continue;
      w[static_cast<std::size_t>(i)] *= nodes[static_cast<std::size_t>(i)] -
                                        nodes[static_cast<std::size_t>(k)];
    }
    w[static_cast<std::size_t>(i)] = Real{1} / w[static_cast<std::size_t>(i)];
  }

  auto d = std::vector<Real>(static_cast<std::size_t>(n * n), Real{0});
  for (auto j = Int{0}; j < n; j++) {
    auto diagonal = Real{0};
    for (auto i = Int{0}; i < n; i++) {
      if (i == j) continue;
      const auto value =
          (w[static_cast<std::size_t>(i)] / w[static_cast<std::size_t>(j)]) /
          (nodes[static_cast<std::size_t>(j)] -
           nodes[static_cast<std::size_t>(i)]);
      d[static_cast<std::size_t>(j * n + i)] = value;
      diagonal -= value;
    }
    d[static_cast<std::size_t>(j * n + j)] = diagonal;
  }
  return d;
}

}  // namespace RadialDetails

//--------------------------------------------------------------------------//
//                             Finite differences                            //
//--------------------------------------------------------------------------//

/// d/dr by finite differences of a chosen order of accuracy, on whatever
/// spacing the radial grid has.
///
/// The stencil is centred where there is room and one-sided at the ends, which
/// is what makes it usable on a grid that has ends -- a centred rule alone
/// would leave the first and last radii undefined, and those are exactly the
/// radii a boundary condition is applied at.
///
/// The weights are Fornberg's, so unequal spacing costs nothing extra and the
/// rule is exact for polynomials of degree at most `order`. This is the one to
/// reach for by default: it is banded, so a line costs `nR * (order + 1)`
/// multiplications, and it is the only one here whose cost does not grow with
/// the number of radii.
template <RealFloatingPoint _Real>
class FiniteDifferenceDerivative {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.
  using Real = _Real;  ///< The precision.

  FiniteDifferenceDerivative() = delete;

  /**
   * @brief Builds the stencils for a grid.
   * @param radial The radii, which need not be uniformly spaced.
   * @param order The order of accuracy, giving a stencil of `order + 1`
   * points.
   */
  explicit FiniteDifferenceDerivative(RadialGrid<Real> radial, Int order = 2)
      : _radial{std::move(radial)}, _width{order + 1} {
    const auto nR = _radial.NumberOfRadii();
    if (order < 1) {
      throw std::invalid_argument(
          "A finite-difference derivative needs an order of at least one");
    }
    if (_width > nR) {
      throw std::invalid_argument(
          "A finite-difference derivative of order " + std::to_string(order) +
          " needs a stencil of " + std::to_string(_width) +
          " radii, but this grid has " + std::to_string(nR));
    }

    // Everything below depends on the radii alone, so it is done once.
    const auto radii = _radial.Radii();
    _first.resize(static_cast<std::size_t>(nR));
    _weights.resize(static_cast<std::size_t>(nR * _width));
    for (auto j = Int{0}; j < nR; j++) {
      // Centred where there is room, shifted just enough at the ends.
      auto first = j - _width / 2;
      if (first < 0) first = 0;
      if (first > nR - _width) first = nR - _width;
      _first[static_cast<std::size_t>(j)] = first;

      const auto stencil = radii.subspan(static_cast<std::size_t>(first),
                                         static_cast<std::size_t>(_width));
      const auto w = RadialDetails::FirstDerivativeWeights<Real>(
          radii[static_cast<std::size_t>(j)], stencil);
      for (auto k = Int{0}; k < _width; k++) {
        _weights[static_cast<std::size_t>(j * _width + k)] =
            w[static_cast<std::size_t>(k)];
      }
    }
  }

  /** @brief The radial grid this is defined on. */
  const RadialGrid<Real>& Radial() const { return _radial; }
  /** @brief The order of accuracy. */
  Int Order() const { return _width - 1; }

  /**
   * @brief Differentiates one radial line.
   * @param in One value per radius.
   * @param out Where the derivative goes, of the same length.
   * @throws std::invalid_argument if either is not one value per radius.
   */
  template <typename Scalar>
  void operator()(std::span<const Scalar> in, std::span<Scalar> out) const {
    const auto nR = _radial.NumberOfRadii();
    CheckLine(in.size(), out.size(), nR);
    for (auto j = Int{0}; j < nR; j++) {
      const auto first = _first[static_cast<std::size_t>(j)];
      auto sum = Scalar{};
      for (auto k = Int{0}; k < _width; k++) {
        sum += _weights[static_cast<std::size_t>(j * _width + k)] *
               in[static_cast<std::size_t>(first + k)];
      }
      out[static_cast<std::size_t>(j)] = sum;
    }
  }

 private:
  RadialGrid<Real> _radial;
  Int _width;
  std::vector<Int> _first;
  std::vector<Real> _weights;

  static void CheckLine(std::size_t in, std::size_t out, Int nR) {
    if (in != static_cast<std::size_t>(nR) ||
        out != static_cast<std::size_t>(nR)) {
      throw std::invalid_argument(
          "A radial operator acts on a line of one value per radius");
    }
  }
};

//--------------------------------------------------------------------------//
//                       The global polynomial derivative                    //
//--------------------------------------------------------------------------//

/// d/dr of the polynomial of degree nR - 1 through all the samples: the
/// differentiation matrix of the given nodes.
///
/// Exact for polynomials of degree at most nR - 1, which is the highest any
/// operator on these nodes can be, and the right thing on a *few* nodes --
/// the Gauss-Lobatto points of one spectral element, say. It is the wrong
/// thing on many: the cost is nR^2 per line where the finite-difference rule
/// is nR * (order + 1), and global polynomial interpolation through equally
/// spaced points diverges as their number grows.
///
/// Built barycentrically, which is the numerically sound way to form the
/// matrix, and once: the entries depend on the nodes alone.
template <RealFloatingPoint _Real>
class LagrangeDerivative {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.
  using Real = _Real;  ///< The precision.

  LagrangeDerivative() = delete;

  /**
   * @brief Builds the differentiation matrix for a grid.
   * @param radial The radii, of which there must be at least two.
   */
  explicit LagrangeDerivative(RadialGrid<Real> radial)
      : _radial{std::move(radial)} {
    const auto nR = _radial.NumberOfRadii();
    if (nR < 2) {
      throw std::invalid_argument(
          "A differentiation matrix needs at least two radii");
    }
    const auto radii = _radial.Radii();

    _d = RadialDetails::DifferentiationMatrix<Real>(radii);
  }

  /** @brief The radial grid this is defined on. */
  const RadialGrid<Real>& Radial() const { return _radial; }

  /// The matrix itself, row-major, for a caller who wants to compose it with
  /// something rather than apply it.
  std::span<const Real> Matrix() const { return std::span<const Real>(_d); }

  /**
   * @brief Differentiates one radial line.
   * @param in One value per radius.
   * @param out Where the derivative goes, of the same length.
   * @throws std::invalid_argument if either is not one value per radius.
   */
  template <typename Scalar>
  void operator()(std::span<const Scalar> in, std::span<Scalar> out) const {
    const auto nR = _radial.NumberOfRadii();
    if (in.size() != static_cast<std::size_t>(nR) ||
        out.size() != static_cast<std::size_t>(nR)) {
      throw std::invalid_argument(
          "A radial operator acts on a line of one value per radius");
    }
    for (auto j = Int{0}; j < nR; j++) {
      auto sum = Scalar{};
      for (auto i = Int{0}; i < nR; i++) {
        sum += _d[static_cast<std::size_t>(j * nR + i)] *
               in[static_cast<std::size_t>(i)];
      }
      out[static_cast<std::size_t>(j)] = sum;
    }
  }

 private:
  RadialGrid<Real> _radial;
  std::vector<Real> _d;
};

//--------------------------------------------------------------------------//
//                          The element derivative                           //
//--------------------------------------------------------------------------//

/// d/dr element by element: block-diagonal, each block the differentiation
/// matrix of its own nodes.
///
/// This is what a spectral element does. Within an element the field is the
/// polynomial through its nodes -- the Gauss-Lobatto-Legendre points, in
/// practice, though nothing here requires that -- so the derivative there is
/// that polynomial's, which is `LagrangeDerivative` restricted to the element.
/// The blocks do not talk to each other, and that is the point rather than an
/// approximation: the field is not assumed differentiable across an interface,
/// because at a material interface it is not.
///
/// **What makes this well defined is the element structure.** The blocks
/// are disjoint, so two
/// elements meet at a repeated radius and each owns one of the pair. The
/// derivative at an interface is therefore two numbers, one per side, each
/// stored at its own index -- which is what a discontinuity *is*, and is the
/// same pair `Interpolation::Piecewise::Limits` hands back. Had the elements
/// shared a node this operator would have had to choose between them, and that
/// choice belongs to the discretisation rather than here.
///
/// Needs a grid built by `RadialGrid::WithElements`; a grid that does not know
/// its elements cannot say what the blocks are, and guessing is precisely what
/// the partition exists to stop.
template <RealFloatingPoint _Real>
class ElementDerivative {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout the library.
  using Real = _Real;  ///< The precision.

  ElementDerivative() = delete;

  /**
   * @brief Builds one differentiation matrix per element.
   * @param radial The radii, which must know their elements.
   * @throws std::invalid_argument if the grid has no elements.
   */
  explicit ElementDerivative(RadialGrid<Real> radial)
      : _radial{std::move(radial)} {
    if (!_radial.HasElements()) {
      throw std::invalid_argument(
          "An element derivative needs a grid that knows its elements, which "
          "is what RadialGrid::WithElements builds; a plain grid is a list of "
          "radii and cannot say where one element ends and the next begins");
    }

    // One matrix per element, all node-dependent and so all built here.
    // Stored end to end, since the elements need not be the same size.
    const auto radii = _radial.Radii();
    _offset.resize(static_cast<std::size_t>(_radial.ElementCount() + 1));
    _offset[0] = 0;
    for (auto k : _radial.ElementIndices()) {
      const auto size = _radial.ElementSize(k);
      const auto nodes = radii.subspan(
          static_cast<std::size_t>(_radial.ElementStart(k)),
          static_cast<std::size_t>(size));
      const auto block = RadialDetails::DifferentiationMatrix<Real>(nodes);
      _d.insert(_d.end(), block.begin(), block.end());
      _offset[static_cast<std::size_t>(k + 1)] =
          static_cast<Int>(_d.size());
    }
  }

  /** @brief The radial grid this is defined on. */
  const RadialGrid<Real>& Radial() const { return _radial; }

  /// One element's matrix, row-major over its own nodes.
  std::span<const Real> Matrix(Int k) const {
    const auto first = static_cast<std::size_t>(_offset[static_cast<std::size_t>(k)]);
    const auto last = static_cast<std::size_t>(_offset[static_cast<std::size_t>(k + 1)]);
    return std::span<const Real>(_d).subspan(first, last - first);
  }

  /**
   * @brief Differentiates one radial line.
   * @param in One value per radius.
   * @param out Where the derivative goes, of the same length.
   * @throws std::invalid_argument if either is not one value per radius.
   */
  template <typename Scalar>
  void operator()(std::span<const Scalar> in, std::span<Scalar> out) const {
    const auto nR = _radial.NumberOfRadii();
    if (in.size() != static_cast<std::size_t>(nR) ||
        out.size() != static_cast<std::size_t>(nR)) {
      throw std::invalid_argument(
          "A radial operator acts on a line of one value per radius");
    }

    for (auto k : _radial.ElementIndices()) {
      const auto first = _radial.ElementStart(k);
      const auto size = _radial.ElementSize(k);
      const auto block = Matrix(k);
      for (auto j = Int{0}; j < size; j++) {
        auto sum = Scalar{};
        for (auto i = Int{0}; i < size; i++) {
          sum += block[static_cast<std::size_t>(j * size + i)] *
                 in[static_cast<std::size_t>(first + i)];
        }
        out[static_cast<std::size_t>(first + j)] = sum;
      }
    }
  }

 private:
  RadialGrid<Real> _radial;
  std::vector<Real> _d;
  std::vector<Int> _offset;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_RADIAL_DERIVATIVES_GUARD_H
