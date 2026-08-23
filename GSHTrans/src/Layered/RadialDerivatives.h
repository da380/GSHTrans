#ifndef GSH_TRANS_RADIAL_DERIVATIVES_GUARD_H
#define GSH_TRANS_RADIAL_DERIVATIVES_GUARD_H

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
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

}  // namespace RadialDetails

//--------------------------------------------------------------------------//
//                             Finite differences                            //
//--------------------------------------------------------------------------//

// d/dr by finite differences of a chosen order of accuracy, on whatever
// spacing the radial grid has.
//
// The stencil is centred where there is room and one-sided at the ends, which
// is what makes it usable on a grid that has ends -- a centred rule alone
// would leave the first and last radii undefined, and those are exactly the
// radii a boundary condition is applied at.
//
// The weights are Fornberg's, so unequal spacing costs nothing extra and the
// rule is exact for polynomials of degree at most `order`. This is the one to
// reach for by default: it is banded, so a line costs `nR * (order + 1)`
// multiplications, and it is the only one here whose cost does not grow with
// the number of radii.
template <RealFloatingPoint _Real>
class FiniteDifferenceDerivative {
 public:
  using Int = std::ptrdiff_t;
  using Real = _Real;

  FiniteDifferenceDerivative() = delete;

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

  const RadialGrid<Real>& Radial() const { return _radial; }
  Int Order() const { return _width - 1; }

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

// d/dr of the polynomial of degree nR - 1 through all the samples: the
// differentiation matrix of the given nodes.
//
// Exact for polynomials of degree at most nR - 1, which is the highest any
// operator on these nodes can be, and the right thing on a *few* nodes --
// the Gauss-Lobatto points of one spectral element, say. It is the wrong
// thing on many: the cost is nR^2 per line where the finite-difference rule
// is nR * (order + 1), and global polynomial interpolation through equally
// spaced points diverges as their number grows.
//
// Built barycentrically, which is the numerically sound way to form the
// matrix, and once: the entries depend on the nodes alone.
template <RealFloatingPoint _Real>
class LagrangeDerivative {
 public:
  using Int = std::ptrdiff_t;
  using Real = _Real;

  LagrangeDerivative() = delete;

  explicit LagrangeDerivative(RadialGrid<Real> radial)
      : _radial{std::move(radial)} {
    const auto nR = _radial.NumberOfRadii();
    if (nR < 2) {
      throw std::invalid_argument(
          "A differentiation matrix needs at least two radii");
    }
    const auto radii = _radial.Radii();

    // Barycentric weights, w_i = 1 / prod_{k != i} (r_i - r_k).
    auto w = std::vector<Real>(static_cast<std::size_t>(nR), Real{1});
    for (auto i = Int{0}; i < nR; i++) {
      for (auto k = Int{0}; k < nR; k++) {
        if (k == i) continue;
        w[static_cast<std::size_t>(i)] *=
            radii[static_cast<std::size_t>(i)] - radii[static_cast<std::size_t>(k)];
      }
      w[static_cast<std::size_t>(i)] = Real{1} / w[static_cast<std::size_t>(i)];
    }

    // D(j, i) = (w_i / w_j) / (r_j - r_i) off the diagonal, and the diagonal
    // is minus the row sum: the negative-sum trick. It is the statement that a
    // constant differentiates to zero, imposed rather than evaluated from a
    // closed form -- so the error in differentiating a constant is the
    // rounding of one sum rather than the accuracy of that form, which is
    // what would otherwise grow with the number of nodes. Applying the matrix
    // sums in a different order, so the answer is a few epsilon and not
    // identically zero.
    _d.assign(static_cast<std::size_t>(nR * nR), Real{0});
    for (auto j = Int{0}; j < nR; j++) {
      auto diagonal = Real{0};
      for (auto i = Int{0}; i < nR; i++) {
        if (i == j) continue;
        const auto value =
            (w[static_cast<std::size_t>(i)] / w[static_cast<std::size_t>(j)]) /
            (radii[static_cast<std::size_t>(j)] -
             radii[static_cast<std::size_t>(i)]);
        _d[static_cast<std::size_t>(j * nR + i)] = value;
        diagonal -= value;
      }
      _d[static_cast<std::size_t>(j * nR + j)] = diagonal;
    }
  }

  const RadialGrid<Real>& Radial() const { return _radial; }

  // The matrix itself, row-major, for a caller who wants to compose it with
  // something rather than apply it.
  std::span<const Real> Matrix() const { return std::span<const Real>(_d); }

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
//                            The spline derivative                          //
//--------------------------------------------------------------------------//

// d/dr of the natural cubic spline through the samples.
//
// The middle ground between the two above, and the one to reach for when the
// radii are many and unevenly spaced: it is global, like the differentiation
// matrix, but its cost is linear in nR rather than quadratic, and unlike a
// global polynomial it does not fall apart as the nodes multiply.
//
// **Why this is written out rather than built on Interpolation, which has a
// perfectly good CubicSpline.** That class computes its coefficients at
// construction and holds them in a vector, so one costs about four
// allocations -- and a radial operator is called once per line, of order
// 66,000 times per application at lMax = 256. Ten milliseconds of allocation
// against a two-millisecond transform is not a trade worth making in the
// inner loop of a matrix-free solve.
//
// The way out is what RadialOperator.h's contract asks for. The natural
// spline's second derivatives satisfy a tridiagonal system whose *matrix
// depends on the radii alone*; only its right-hand side carries the data. So
// the matrix is factorised once, here, and a line costs one Thomas sweep into
// thread_local scratch and no allocation at all. That is a thing upstream's
// interface cannot express, which is why this is not duplication -- and
// `Interpolation::CubicSpline` is the oracle the tests check it against,
// which is a stronger check than any property this file could asssert about
// itself (field-algebra-plan.md section 19.5 [R4]).
//
// Natural conditions, S'' = 0 at both ends. That is a choice and it shows: it
// is wrong for data that is genuinely curved at the boundary, and the first
// and last intervals are where a spline derivative is least trustworthy
// whatever conditions are imposed. A caller who knows the end slopes has a
// better operator available in three lines of their own, which is the point
// of the seam.
template <RealFloatingPoint _Real>
class SplineDerivative {
 public:
  using Int = std::ptrdiff_t;
  using Real = _Real;

  SplineDerivative() = delete;

  explicit SplineDerivative(RadialGrid<Real> radial)
      : _radial{std::move(radial)} {
    const auto nR = _radial.NumberOfRadii();
    if (nR < 2) {
      throw std::invalid_argument("A spline derivative needs at least two radii");
    }
    const auto radii = _radial.Radii();
    const auto n = static_cast<std::size_t>(nR);

    _h.resize(n - 1);
    for (std::size_t i = 0; i + 1 < n; i++) {
      _h[i] = radii[i + 1] - radii[i];
      if (_h[i] <= Real{0}) {
        throw std::invalid_argument(
            "A spline derivative needs strictly increasing radii, and this "
            "grid repeats one -- which is how a material interface is "
            "written, and is a case for a piecewise operator rather than "
            "this one");
      }
    }

    // The system for the second derivatives m, with m = 0 at both ends. Only
    // the three diagonals live here; the right-hand side is the data and is
    // built per line.
    _sub.assign(n, Real{0});
    _diag.assign(n, Real{1});
    _super.assign(n, Real{0});
    for (std::size_t i = 1; i + 1 < n; i++) {
      _sub[i] = _h[i - 1] / 6;
      _diag[i] = (_h[i - 1] + _h[i]) / 3;
      _super[i] = _h[i] / 6;
    }
  }

  const RadialGrid<Real>& Radial() const { return _radial; }

  template <typename Scalar>
  void operator()(std::span<const Scalar> in, std::span<Scalar> out) const {
    const auto n = static_cast<std::size_t>(_radial.NumberOfRadii());
    if (in.size() != n || out.size() != n) {
      throw std::invalid_argument(
          "A radial operator acts on a line of one value per radius");
    }

    // Scratch, not state: one operator serves every thread, so anything the
    // call writes to lives here (RadialOperator.h).
    thread_local auto diagonal = std::vector<Real>{};
    thread_local auto m = std::vector<Scalar>{};
    if (diagonal.size() < n) diagonal.resize(n);
    if (m.size() < n) m.resize(n);

    for (std::size_t i = 0; i < n; i++) diagonal[i] = _diag[i];
    m[0] = Scalar{};
    m[n - 1] = Scalar{};
    for (std::size_t i = 1; i + 1 < n; i++) {
      m[i] = (in[i + 1] - in[i]) / _h[i] - (in[i] - in[i - 1]) / _h[i - 1];
    }

    // Thomas, in place. No pivoting, and none needed: an interior row has
    // off-diagonal magnitude (h[i-1] + h[i]) / 6 against a diagonal of
    // (h[i-1] + h[i]) / 3, so the system is strictly diagonally dominant.
    for (std::size_t i = 1; i < n; i++) {
      const auto factor = _sub[i] / diagonal[i - 1];
      diagonal[i] -= factor * _super[i - 1];
      m[i] -= factor * m[i - 1];
    }
    m[n - 1] = m[n - 1] / diagonal[n - 1];
    for (auto i = static_cast<Int>(n) - 2; i >= 0; i--) {
      const auto k = static_cast<std::size_t>(i);
      m[k] = (m[k] - _super[k] * m[k + 1]) / diagonal[k];
    }

    // S'(r_j) from the left end of the interval that starts there, and from
    // the right end of the last interval for the final node.
    for (std::size_t j = 0; j + 1 < n; j++) {
      out[j] = (in[j + 1] - in[j]) / _h[j] -
               _h[j] * (Real{2} * m[j] + m[j + 1]) / Real{6};
    }
    const auto last = n - 2;
    out[n - 1] = (in[n - 1] - in[last]) / _h[last] +
                 _h[last] * (m[last] + Real{2} * m[n - 1]) / Real{6};
  }

 private:
  RadialGrid<Real> _radial;
  std::vector<Real> _h;
  std::vector<Real> _sub, _diag, _super;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_RADIAL_DERIVATIVES_GUARD_H
