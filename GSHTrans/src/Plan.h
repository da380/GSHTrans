#ifndef GSH_TRANS_PLAN_GUARD_H
#define GSH_TRANS_PLAN_GUARD_H

#include <FFTWpp/All>
#include <GaussQuad/All>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <iterator>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <ranges>
#include <variant>
#include <vector>

#include "Concepts.h"
#include "Indexing.h"
#include "Legendre.h"
#include "Wigner.h"

namespace GSHTrans {

template <RealFloatingPoint Real = double, TransformType Type = C2C,
          Normalisation Norm = Ortho>
class Plan {
  // Local type aliases.
  using Complex = std::complex<Real>;
  using WignerType = Wigner<Real, typename Type::Orders, Norm>;
  using WignerPointer = std::unique_ptr<WignerType>;
  using WignerVector = std::vector<WignerPointer>;

 public:
  Plan() = default;

  Plan(int lMax, int nMax, FFTWpp::PlanFlag flag = FFTWpp::Measure)
      : _lMax{lMax},
        _nMax{nMax},
        _quad{GaussQuad::GaussLegendreQuadrature1D<Real>(_lMax + 1)} {
    // Check the inputs.
    assert(MaxDegree() >= 0);
    assert(MaxUpperIndex() >= 0 && MaxUpperIndex() <= MaxDegree());

    // Transform the quadrature points.
    _quad.Transform([](auto x) { return std::acos(-x); },
                    [](auto x) { return 1; });

    // Compute the Wigner d-functions.
    _d = WignerVector(MaxUpperIndex() - MinUpperIndex() + 1);
    for (auto n : UpperIndices()) {
      _d[n - MinUpperIndex()] = std::make_unique<WignerType>(
          MaxDegree(), MaxDegree(), n, _quad.Points());
    }

    // Generate wisdom for FFTs.
    if constexpr (std::same_as<Type, C2C>) {
      _inLayout = FFTWpp::DataLayout(1, std::vector{NumberOfLongitudes()}, 1,
                                     std::vector{NumberOfLongitudes()}, 1, 1);
      _outLayout = _inLayout;
      FFTWpp::GenerateWisdom<Complex, Complex, true>(_inLayout, _outLayout,
                                                     flag);
    } else {
      _inLayout = FFTWpp::DataLayout(1, std::vector{NumberOfLongitudes()}, 1,
                                     std::vector{NumberOfLongitudes()}, 1, 1);
      _outLayout =
          FFTWpp::DataLayout(1, std::vector{NumberOfLongitudes() / 2 + 1}, 1,
                             std::vector{NumberOfLongitudes() / 2 + 1}, 1, 1);
      FFTWpp::GenerateWisdom<Real, Complex, true>(_inLayout, _outLayout, flag);
    }
  }

  // Return the truncation degree.
  auto MaxDegree() const { return _lMax; }

  // Return view to the degrees.
  auto Degrees() const { return GSHTrans::Degrees(_lMax); }

  // Return the maximum upper index.
  auto MaxUpperIndex() const { return _nMax; }

  auto MinUpperIndex() const {
    if constexpr (std::same_as<Type, C2C>) {
      return -_nMax;
    } else {
      return 0;
    }
  };

  // Return view to the upper indices
  auto UpperIndices() const { return GSHTrans::UpperIndices<Type>(_nMax); }

  // Return the longitude spacing.
  auto DeltaLongitude() const {
    return 2 * std::numbers::pi_v<Real> / static_cast<Real>(2 * _lMax);
  }

  // Return numbers of points in each direction.
  auto NumberOfCoLatitudes() const { return _lMax + 1; }
  auto NumberOfLongitudes() const { return 2 * _lMax; }

  // Return view to the co-latitudes.
  auto CoLatitudes() const { return std::ranges::views::all(_quad.Points()); }

  // Return view to the longitudes.
  auto Longitudes() const {
    auto dPhi = DeltaLongitude();
    return std::ranges::views::iota(0, 2 * _lMax) |
           std::ranges::views::transform([dPhi](auto i) { return i * dPhi; });
  }

  // View over colatitude index.
  auto IndexCoLatitudes() const {
    return std::ranges::views::iota(0, _lMax + 1);
  }

  // View over longitude index
  auto IndexLongitudes() const {
    return std::ranges::views::iota(0, 2 * _lMax);
  }

  // Integrate a function over the unit sphere.
  template <typename Function>
  requires requires(Function f, Real theta, Real phi, Real w) {
    std::invocable<Function, Real, Real>;
    {f(theta, phi) * w};
  }
  auto Integrate(Function f) {
    using FunctionValue = decltype(f(0, 0));
    auto thetaIntegrand = [this, &f](auto theta) {
      auto dPhi = DeltaLongitude();
      auto phi = Longitudes();
      return dPhi * std::accumulate(phi.begin(), phi.end(), FunctionValue{0},
                                    [&theta, &f](auto acc, auto phi) {
                                      return acc + f(theta, phi);
                                    });
    };
    return _quad.Integrate(thetaIntegrand);
  };

  // Forward  transformation.
  template <std::ranges::random_access_range RangeIn,
            std::ranges::random_access_range RangeOut>
  void Execute(RangeIn& in, RangeOut& out, int n, int lMax) {
    const auto nPhi = NumberOfLongitudes();

    // Initialise tempory vector to store FFT results.
    auto tmp = std::vector<Complex>(nPhi);
    auto tmpView = FFTWpp::DataView(tmp.begin(), tmp.end(), _outLayout);

    // Set scale factor for nomalising the FFTs
    const auto scaleFactor =
        2 * std::numbers::pi_v<Real> / static_cast<Real>(nPhi);

    // Iterator to the start of the in-data.
    auto inStart = in.begin();
    for (auto w : _quad.Weights()) {
      // Scale the quadrature weight
      w *= scaleFactor;

      // Iterator to the end of the in-data for this colatitude.
      auto inFinish = std::next(inStart, nPhi);

      // Make views to the FFT data.
      auto inView = FFTWpp::DataView(inStart, inFinish, _inLayout);

      // Form the FFT plan.
      auto FFTPlan =
          FFTWpp::Plan(inView, tmpView, FFTWpp::WisdomOnly, FFTWpp::Forward);

      // Perform the FFT
      FFTPlan.Execute();

      // Iterator to the start of the coefficients.
      auto coefficientIterator = out.begin();

      // Iterator to the start to the Wigner functions.
      auto wignerIterator = _d[n - MinUpperIndex()]->begin();

      // Loop over the spherical harmonic coefficients.
      /*
      for (auto [l, m] : SphericalHarmonicIndices(lMax)) {
        auto i = m < 0 ? nPhi + m : m;
        *coefficientIterator++ += *wignerIterator++ * tmp[i] * w;
      }
      */

      // Set the start for the next section of the in-data.
      inStart = inFinish;
    }
  }

  template <std::ranges::random_access_range RangeIn,
            std::ranges::random_access_range RangeOut>
  void Execute(RangeIn& in, RangeOut& out,
               int n) requires std::same_as<Type, C2C> {
    Execute(in, out, n, MaxDegree());
  }

 private:
  // Store maximum orders and upper indices.
  int _lMax;
  int _nMax;

  // Store the data layouts
  FFTWpp::DataLayout _inLayout;
  FFTWpp::DataLayout _outLayout;

  // Store quadrature points and weights.
  GaussQuad::Quadrature1D<Real> _quad;

  // Vector of pointers to the Wigner values for difference upper indices.
  WignerVector _d;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_PLAN_GUARD_H
