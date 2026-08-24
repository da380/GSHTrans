#ifndef GSH_TRANS_INTERPOLATE_GUARD_H
#define GSH_TRANS_INTERPOLATE_GUARD_H

// Interpolation of a field: a callable of the two angles.
//
// The plan is field-algebra-plan.md section 22, and the decisions it records
// are referred to below as [I1] to [I9]. The short version of what makes this
// more than a forwarding call to Interpolation:
//
//   -- the longitudes run 0 ... 2pi - dphi, so the wrap is not represented and
//      a query in the last cell would interpolate against nothing;
//   -- neither pole is a grid point, so every local scheme extrapolates there,
//      and does so exactly where a spin-weighted field's sin^|m|(theta)
//      behaviour is most delicate.
//
// Both are fixed by handing the upstream scheme a *padded* grid rather than
// the field's own, and the padding is exact: the wrap column is column zero
// copied, and the two polar rows are computed from the expansion.

#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <memory>
#include <numbers>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "../Concepts.h"
#include "../Policies.h"
#include "../SpinField/SpinWeighted.h"

namespace GSHTrans {

namespace InterpolateDetails {

//--------------------------------------------------------------------------//
//                      What interpolation needs of a grid                   //
//--------------------------------------------------------------------------//

// AngularGrid asks only for Points(), because the expression layer evaluates
// point by point and never needs to know that the point set is a product. An
// interpolant does: a rectilinear scheme takes one abscissa range per axis,
// and the padding has to build each of them separately. Stating the extra
// requirement as a concept documents the coupling, in the same way
// SpinWeighted.h states the base one.
//
// Written against range_value_t rather than against *begin(...) deliberately:
// an axis accessor returns a view by value, and views::all of a prvalue
// container is an owning_view, which is not a borrowed range -- so
// ranges::begin on the returned prvalue is ill-formed and the concept would be
// unsatisfiable for a grid that in fact has both axes. Found by the concept
// rejecting GaussLegendreGrid.
template <typename G>
concept SeparableAngularGrid = AngularGrid<G> && requires(const G& grid) {
  { grid.CoLatitudes() } -> std::ranges::range;
  { grid.Longitudes() } -> std::ranges::range;
  requires std::convertible_to<
      std::ranges::range_value_t<decltype(grid.CoLatitudes())>,
      typename G::Real>;
  requires std::convertible_to<
      std::ranges::range_value_t<decltype(grid.Longitudes())>,
      typename G::Real>;
};

//--------------------------------------------------------------------------//
//                              The padded grid                              //
//--------------------------------------------------------------------------//

// The field's samples on a grid that covers the closed sphere: colatitudes
// running 0 ... pi and longitudes running 0 ... 2pi, both inclusive, with the
// values row-major and theta-major so that element (i, j) sits at
// i * Columns() + j -- which is both this library's layout (SpinField.h) and
// Interpolation's (Bilinear.hpp), and is why no repack is needed.
//
// It owns its arrays. [I1]: an interpolant is a snapshot, and here there is
// nothing to borrow anyway, since every array below is new storage whatever
// the caller passed.
template <RealFloatingPoint _Real, RealOrComplexFloatingPoint _Scalar>
struct Padded {
  using Real = _Real;
  using Scalar = _Scalar;

  std::vector<Real> theta;     // nTheta + 2, running 0 ... pi
  std::vector<Real> phi;       // nPhi + 1,   running 0 ... 2 pi
  std::vector<Scalar> values;  // Rows() * Columns(), row-major

  auto Rows() const { return theta.size(); }
  auto Columns() const { return phi.size(); }

  Scalar At(std::size_t i, std::size_t j) const {
    return values[i * Columns() + j];
  }
};

// Build one from a field's samples and the two polar rows.
//
// The polar rows are arguments rather than being computed here, so that the
// padding can be tested without an expansion -- which is what makes P1 of
// section 22.3 a step of its own. Each is nPhi values at the grid's own
// longitudes; the wrap column is added to them exactly as it is to every
// interior row, which is right because exp(i N 2pi) = exp(i N 0) for integer
// N and so the polar row closes on itself like any other.
template <SeparableAngularGrid GridType, RealOrComplexFloatingPoint Scalar>
auto Pad(const GridType& grid, std::span<const Scalar> samples,
         std::span<const Scalar> north, std::span<const Scalar> south) {
  using Real = typename GridType::Real;
  constexpr auto pi = std::numbers::pi_v<Real>;

  auto padded = Padded<Real, Scalar>{};

  for (auto t : grid.CoLatitudes()) padded.theta.push_back(t);
  for (auto p : grid.Longitudes()) padded.phi.push_back(p);

  const auto nTheta = padded.theta.size();
  const auto nPhi = padded.phi.size();

  if (nTheta == 0 || nPhi == 0) {
    throw std::invalid_argument("Interpolate: the grid has no points");
  }
  if (samples.size() != nTheta * nPhi) {
    throw std::invalid_argument(
        "Interpolate: expected " + std::to_string(nTheta * nPhi) +
        " samples, but was given " + std::to_string(samples.size()));
  }
  if (north.size() != nPhi || south.size() != nPhi) {
    throw std::invalid_argument(
        "Interpolate: a polar row must hold one value per longitude");
  }

  // Checked rather than assumed. Gauss-Legendre nodes are interior to
  // (0, pi), so this holds -- but if it ever did not, the padded axis would
  // stop being strictly increasing and upstream would refuse it with a
  // message about abscissae rather than about poles.
  if (!(padded.theta.front() > 0) || !(padded.theta.back() < pi)) {
    throw std::invalid_argument(
        "Interpolate: the grid's colatitudes must lie strictly inside "
        "(0, pi), so that the poles can be added to them");
  }
  if (!(padded.phi.front() == 0)) {
    throw std::invalid_argument(
        "Interpolate: the grid's longitudes must start at zero");
  }

  padded.theta.insert(padded.theta.begin(), Real{0});
  padded.theta.push_back(pi);
  padded.phi.push_back(2 * pi);

  const auto columns = padded.phi.size();
  padded.values.reserve(padded.theta.size() * columns);

  auto pushRow = [&padded](std::span<const Scalar> row) {
    padded.values.insert(padded.values.end(), row.begin(), row.end());
    padded.values.push_back(row[0]);  // the wrap: phi = 2 pi is phi = 0
  };

  pushRow(north);
  for (std::size_t i = 0; i < nTheta; i++) {
    pushRow(samples.subspan(i * nPhi, nPhi));
  }
  pushRow(south);

  return padded;
}

}  // namespace InterpolateDetails

}  // namespace GSHTrans

#endif  // GSH_TRANS_INTERPOLATE_GUARD_H
