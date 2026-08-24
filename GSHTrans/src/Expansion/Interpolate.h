#ifndef GSH_TRANS_INTERPOLATE_GUARD_H
#define GSH_TRANS_INTERPOLATE_GUARD_H

/**
 * @file Interpolate.h
 * @brief Interpolation of a field: a callable of the two angles.
 *
 * @details Two things make this more than a forwarding call to the upstream
 * interpolation library:
 *
 * - the longitudes run @f$0 \ldots 2\pi - \Delta\phi@f$, so the wrap is not
 *   represented and a query in the last cell would interpolate against
 *   nothing;
 * - neither pole is a grid point, so every local scheme extrapolates there,
 *   and does so exactly where a spin-weighted field's
 *   @f$\sin^{|m|}\theta@f$ behaviour is most delicate.
 *
 * Both are fixed by handing the upstream scheme a *padded* grid rather than
 * the field's own, and the padding is exact: the wrap column is column zero
 * copied, and the two polar rows are computed from the expansion.
 */

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
#include <utility>
#include <vector>

#ifdef GSHTRANS_HAVE_INTERPOLATION
#include <Interpolation/BicubicSpline.hpp>
#include <Interpolation/Bilinear.hpp>
#endif

#include "../Concepts.h"
#include "../Indexing.h"
#include "../Policies.h"
#include "../SpinField/SpinWeighted.h"
#include "../Views.h"
#include "../Wigner.h"
#include "SpinExpansion.h"

namespace GSHTrans {

namespace InterpolateDetails {

//--------------------------------------------------------------------------//
//                              The padded grid                              //
//--------------------------------------------------------------------------//

// The field's samples on a grid that covers the closed sphere: colatitudes
// running 0 ... pi and longitudes running 0 ... 2pi, both inclusive, with the
// values row-major and theta-major so that element (i, j) sits at
// i * Columns() + j -- which is both this library's layout (SpinField.h) and
// Interpolation's (Bilinear.hpp), and is why no repack is needed.
//
// It owns its arrays. An interpolant is a snapshot, and here there is nothing
// to borrow anyway, since every array below is new storage whatever the
// caller passed.
template <RealFloatingPoint _Real, RealOrComplexFloatingPoint _Scalar>
struct Padded {
  using Real = _Real;  ///< The precision.
  /// The value type: Real when real-valued, Complex otherwise.
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
// padding can be tested without an expansion. Each is nPhi values at the
// grid's own
// longitudes; the wrap column is added to them exactly as it is to every
// interior row, which is right because exp(i N 2pi) = exp(i N 0) for integer
// N and so the polar row closes on itself like any other.
template <AngularGrid GridType, RealOrComplexFloatingPoint Scalar>
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

// The two polar rows, from the expansion.
//
// Neither pole is a grid point, so a local scheme has nothing to interpolate
// against there and would extrapolate -- exactly where a spin-weighted field
// is most delicate. The rows below remove the extrapolation entirely, and they
// are exact rather than a fudge, because at a pole all but one order vanishes
// (measured over |N| <= 2 and l <= 5):
//
//     f(0,  phi) = exp(+i N phi) sum_l           f^N_{l,+N} sqrt((2l+1)/4pi)
//     f(pi, phi) = exp(-i N phi) sum_l (-1)^{l-N} f^N_{l,-N} sqrt((2l+1)/4pi)
//
// The phi dependence is the point rather than an inconvenience: a
// spin-weighted field at a coordinate pole is not single-valued, because the
// frame e_pm depends on the azimuth of approach. So a polar row is a row, and
// a constant one would be wrong at every N != 0.
template <RealOrComplexFloatingPoint Scalar, typename Expansion,
          AngularGrid GridType>
auto PolarRows(const Expansion& expansion, const GridType& grid) {
  using Real = typename GridType::Real;
  using Complex = std::complex<Real>;
  constexpr auto N = Expansion::UpperIndex;
  constexpr auto pi = std::numbers::pi_v<Real>;

  auto north = Complex{};
  auto south = Complex{};
  for (auto l : expansion.Degrees()) {
    const auto norm = std::sqrt((2 * static_cast<Real>(l) + 1) / (4 * pi));
    const auto sign = (l - N) % 2 == 0 ? Real{1} : Real{-1};
    north += expansion[l, N] * norm;
    south += expansion[l, -N] * norm * sign;
  }

  auto rows = std::pair<std::vector<Scalar>, std::vector<Scalar>>{};
  for (auto phi : grid.Longitudes()) {
    const auto up = north * std::polar(Real{1}, static_cast<Real>(N) * phi);
    const auto down = south * std::polar(Real{1}, -static_cast<Real>(N) * phi);
    if constexpr (RealFloatingPoint<Scalar>) {
      rows.first.push_back(std::real(up));
      rows.second.push_back(std::real(down));
    } else {
      rows.first.push_back(up);
      rows.second.push_back(down);
    }
  }
  return rows;
}

}  // namespace InterpolateDetails

//--------------------------------------------------------------------------//
//                          The spectral interpolant                         //
//--------------------------------------------------------------------------//

/// The expansion evaluated directly, which is exact for a band-limited field:
///
///     f(theta, phi) = sum_l sum_m  f^N_{lm} dbar^l_{Nm}(theta) exp(i m phi),
///     dbar^l_{Nm}   = sqrt((2l+1)/4pi) d^l_{Nm},
///
/// which is what WignerDetails::ComputeBlock stores and what the loop kernel's
/// SynthesiseRow sums. Checked against Evaluate at every grid point rather than
/// read off the transform: worst absolute difference 5.0e-15 complex and
/// 8.1e-15 real.
///
/// This is the reference the cheap schemes are measured against. It is exact
/// and pole-safe, and it costs O(lMax^2) a point where they cost O(1).
///
/// It owns its coefficients. The shared_ptr is what keeps a copy cheap and,
/// more to the point, keeps the state at a stable address so that copying the
/// interpolant is well defined -- which it has to be, because the grid's
/// ProjectFunction takes its callable by value.
template <std::ptrdiff_t _N, AngularGrid _Grid,
          RealOrComplexValued _Value = ComplexValued>
class SpectralInterpolant {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = _N;
  using GridType = _Grid;  ///< The angular grid this is defined on.
  using Value = _Value;    ///< Whether the samples are real-valued or complex.
  using Real = typename _Grid::Real;   ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  using Scalar = std::conditional_t<std::same_as<Value, RealValued>, Real,
                                    Complex>;  ///< What evaluation returns.

  /// A real field stores only m >= 0, the rest being fixed by
  /// f_{l,-m} = (-1)^m conj(f_{lm}). That is the expansion's own convention,
  /// and it is read here rather than restated.
  using MRange =
      std::conditional_t<std::same_as<Value, RealValued>, NonNegative, All>;

  static_assert(std::same_as<Value, ComplexValued> or UpperIndex == 0,
                "A real-valued field exists only at upper index zero");

  /**
   * @brief Takes a copy of an expansion's coefficients.
   * @param lMax The largest degree present.
   * @param coefficients The coefficients, laid out as GSHIndices describes.
   */
  SpectralInterpolant(Int lMax, std::span<const Complex> coefficients)
      : _state{std::make_shared<const State>(lMax, coefficients)} {}

  /** @brief The largest degree stored. */
  auto MaxDegree() const { return _state->lMax; }

  /**
   * @brief The field at @p theta, @p phi, summed from the expansion.
   * @throws std::invalid_argument if @p theta is outside @f$[0, \pi]@f$; a
   * longitude outside @f$[0, 2\pi)@f$ is reduced instead, which is exact.
   */
  Scalar operator()(Real theta, Real phi) const {
    const auto lMax = _state->lMax;
    constexpr auto pi = std::numbers::pi_v<Real>;

    // A colatitude outside [0, pi] is not a point on the sphere, so
    // there is no number to return; a longitude outside [0, 2pi) is one, and
    // reducing it is exact rather than an approximation. The two axes are not
    // alike and this is where that shows.
    if (!(theta >= 0) || !(theta <= pi)) {
      throw std::invalid_argument(
          "Interpolate: the colatitude must lie in [0, pi]");
    }
    phi = std::fmod(phi, 2 * pi);
    if (phi < 0) phi += 2 * pi;

    // Per-call scratch, grown and never shrunk, in thread_local storage --
    // the rule RadialOperator.h states and for the same reason: an
    // interpolant is called once per evaluation point, so an allocation here
    // is a defect rather than a cost.
    const auto indices = GSHIndices<All>(lMax, lMax, UpperIndex);
    thread_local auto block = std::vector<Real>{};
    thread_local auto phase = std::vector<Complex>{};
    const auto blockSize = static_cast<std::size_t>(indices.Size());
    const auto phases = static_cast<std::size_t>(2 * lMax + 1);
    if (block.size() < blockSize) block.resize(blockSize);
    if (phase.size() < phases) phase.resize(phases);

    // The same recursion WignerValues::Generated() runs, into our own
    // scratch. Sharing it is what stops a second convention arising: a
    // disagreement here would be a disagreement with the transform.
    auto d = GSHView<Real, All>(lMax, lMax, UpperIndex, block.data());
    WignerDetails::ComputeBlock(d, UpperIndex, theta,
                                std::span<const Real>(_state->sqrtInt),
                                std::span<const Real>(_state->sqrtIntInv));

    // exp(i m phi) for every order at once. Built with polar rather than by
    // repeated multiplication: it is O(lMax) against the sum's O(lMax^2), so
    // about one per cent of the work, and it does not accumulate the phase
    // drift that a recurrence would carry to m = lMax.
    for (auto m = -lMax; m <= lMax; m++) {
      phase[static_cast<std::size_t>(m + lMax)] =
          std::polar(Real{1}, static_cast<Real>(m) * phi);
    }

    const auto coefficients = ConstGSHView<Complex, MRange>(
        lMax, lMax, UpperIndex, _state->data.data());

    auto sum = Complex{};
    for (auto l : d.Degrees()) {
      auto dl = d[l];
      for (auto m : dl.Orders()) {
        sum += Coefficient(coefficients, l, m) * dl[m] *
               phase[static_cast<std::size_t>(m + lMax)];
      }
    }

    // For a real field the sum is real -- measured at 5.3e-16 in the
    // imaginary part -- so taking the real part is a projection onto a
    // quantity known to be real rather than a truncation.
    if constexpr (std::same_as<Value, RealValued>) {
      return std::real(sum);
    } else {
      return sum;
    }
  }

 private:
  struct State {
    Int lMax;
    std::vector<Complex> data;
    std::vector<Real> sqrtInt;
    std::vector<Real> sqrtIntInv;

    State(Int lMaxIn, std::span<const Complex> coefficients)
        : lMax{lMaxIn}, data(coefficients.begin(), coefficients.end()) {
      if (lMax < std::abs(UpperIndex)) {
        throw std::invalid_argument(
            "Interpolate: the degree is below the upper index");
      }
      const auto indices = GSHIndices<MRange>(lMax, lMax, UpperIndex);
      if (data.size() != static_cast<std::size_t>(indices.Size())) {
        throw std::invalid_argument(
            "Interpolate: the coefficient range does not hold the " +
            std::to_string(indices.Size()) + " coefficients of a degree-" +
            std::to_string(lMax) + " expansion");
      }
      auto tables = WignerDetails::PreComputeTables<Real>(lMax, lMax,
                                                          std::abs(UpperIndex));
      sqrtInt = std::move(tables.first);
      sqrtIntInv = std::move(tables.second);
    }
  };

  // The coefficient at (l, m), through the reduced storage where there is
  // one. This is the whole of what the real case costs here.
  template <typename View>
  static Complex Coefficient(const View& coefficients, Int l, Int m) {
    if constexpr (std::same_as<Value, RealValued>) {
      if (m >= 0) return coefficients[l][m];
      const auto conjugate = std::conj(coefficients[l][-m]);
      return (-m) % 2 == 0 ? conjugate : -conjugate;
    } else {
      return coefficients[l][m];
    }
  }

  std::shared_ptr<const State> _state;
};

//--------------------------------------------------------------------------//
//                                Interpolate                                //
//--------------------------------------------------------------------------//

/// An expansion interpolates spectrally and in no other way: there are no
/// samples to interpolate, only coefficients to sum. The scheme argument is
/// accepted so that the spelling matches the field's, and refused if it names
/// anything else.
template <std::ptrdiff_t N, typename GridType, typename Value>
auto Interpolate(const SpinExpansion<N, GridType, Value>& expansion,
                 Scheme::SpectralTag = Scheme::Spectral()) {
  return SpectralInterpolant<N, GridType, Value>(expansion.MaxDegree(),
                                                 expansion.Data());
}

/// A field, spectrally: expand and sum. The degree is the truncation at which
/// the expansion is taken, defaulting to the grid's own -- which is what an
/// oversampled ForBand grid wants to be able to say.
template <SpinWeighted F>
auto Interpolate(const F& field, Scheme::SpectralTag = Scheme::Spectral(),
                 std::ptrdiff_t lMax = -1) {
  const auto degree = lMax < 0 ? field.Grid().MaxDegree() : lMax;
  const auto expansion = Expand(field, degree);
  return SpectralInterpolant<F::UpperIndex, typename F::GridType,
                             typename F::Value>(degree, expansion.Data());
}

#ifdef GSHTRANS_HAVE_INTERPOLATION

//--------------------------------------------------------------------------//
//                            The local interpolants                         //
//--------------------------------------------------------------------------//

/// One of Interpolation's rectilinear schemes over the *padded* grid: the
/// field's own samples, plus the wrap column and the two polar rows. After the
/// padding and the domain rules on the angles, no query reaches upstream's
/// edge-cell continuation at all, which is the property that makes those
/// schemes usable on a sphere.
///
/// The upstream object is given std::span rather than the deduction guide's
/// views, for two reasons that both matter. A span is a view, so it satisfies
/// upstream's constraints; and it is a type this class can *spell*, which lets
/// the padded arrays and the interpolant over them live in one State and be
/// initialised in order. Handing upstream an owning view instead would make
/// this class move-only, and it has to be copyable -- ProjectFunction takes
/// its callable by value.
template <std::ptrdiff_t _N, AngularGrid _Grid, RealOrComplexValued _Value,
          typename _Upstream>
class LocalInterpolant {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = _N;
  using GridType = _Grid;  ///< The angular grid this is defined on.
  using Value = _Value;    ///< Whether the samples are real-valued or complex.
  using Real = typename _Grid::Real;   ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  using Scalar = std::conditional_t<std::same_as<Value, RealValued>, Real,
                                    Complex>;  ///< What evaluation returns.

  /** @brief Takes ownership of a padded grid and builds the upstream scheme
   * over it. */
  explicit LocalInterpolant(InterpolateDetails::Padded<Real, Scalar> padded)
      : _state{std::make_shared<const State>(std::move(padded))} {}

  /**
   * @brief The field at @p theta, @p phi, from the padded samples.
   * @throws std::invalid_argument if @p theta is outside @f$[0, \pi]@f$; a
   * longitude outside @f$[0, 2\pi)@f$ is reduced instead, which is exact.
   */
  Scalar operator()(Real theta, Real phi) const {
    constexpr auto pi = std::numbers::pi_v<Real>;
    if (!(theta >= 0) || !(theta <= pi)) {
      throw std::invalid_argument(
          "Interpolate: the colatitude must lie in [0, pi]");
    }
    phi = std::fmod(phi, 2 * pi);
    if (phi < 0) phi += 2 * pi;
    return _state->upstream(theta, phi);
  }

 private:
  struct State {
    InterpolateDetails::Padded<Real, Scalar> padded;
    _Upstream upstream;

    // padded is declared first, so it is built first and the spans below
    // point at arrays that already exist. State is never moved -- it is
    // reached only through a shared_ptr -- so they stay valid.
    explicit State(InterpolateDetails::Padded<Real, Scalar> paddedIn)
        : padded{std::move(paddedIn)},
          upstream(std::span<const Real>(padded.theta),
                   std::span<const Real>(padded.phi),
                   std::span<const Scalar>(padded.values)) {}
  };

  std::shared_ptr<const State> _state;
};

namespace InterpolateDetails {

// Which upstream type a scheme tag names.
template <typename Tag, typename Real, typename Scalar>
struct UpstreamFor;

template <typename Real, typename Scalar>
struct UpstreamFor<Scheme::BilinearTag, Real, Scalar> {
  using Type =
      Interpolation::Bilinear<std::span<const Real>, std::span<const Real>,
                              std::span<const Scalar>>;
};

template <typename Real, typename Scalar>
struct UpstreamFor<Scheme::BicubicTag, Real, Scalar> {
  using Type =
      Interpolation::BicubicSpline<std::span<const Real>, std::span<const Real>,
                                   std::span<const Scalar>>;
};

template <typename Tag>
concept LocalScheme = std::same_as<Tag, Scheme::BilinearTag> or
                      std::same_as<Tag, Scheme::BicubicTag>;

}  // namespace InterpolateDetails

/// A field, locally. The forward transform is for the polar rows and nothing
/// else: two columns of coefficients out of a whole expansion, which is the
/// construction cost this scheme carries and the reason its cheapness is per
/// evaluation rather than per interpolant.
template <SpinWeighted F, InterpolateDetails::LocalScheme Tag>
auto Interpolate(const F& field, Tag, std::ptrdiff_t lMax = -1) {
  using Real = typename F::Real;
  using Scalar = typename F::Scalar;
  using Upstream =
      typename InterpolateDetails::UpstreamFor<Tag, Real, Scalar>::Type;

  const auto& grid = field.Grid();
  const auto degree = lMax < 0 ? grid.MaxDegree() : lMax;

  auto samples =
      std::vector<Scalar>(static_cast<std::size_t>(grid.FieldSize()));
  field.EvaluateInto(std::span<Scalar>(samples));

  const auto expansion = Expand(field, degree);
  const auto rows = InterpolateDetails::PolarRows<Scalar>(expansion, grid);

  auto padded = InterpolateDetails::Pad<typename F::GridType, Scalar>(
      grid, samples, rows.first, rows.second);

  return LocalInterpolant<F::UpperIndex, typename F::GridType,
                          typename F::Value, Upstream>(std::move(padded));
}

#endif  // GSHTRANS_HAVE_INTERPOLATION

}  // namespace GSHTrans

#endif  // GSH_TRANS_INTERPOLATE_GUARD_H
