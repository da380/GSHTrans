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
#include "Legendre.h"
#include "Wigner.h"

namespace GSHTrans {

template <std::floating_point Real, TransformType Type = C2C,
          Normalisation Norm = Ortho>
class Plan {
  // Local type aliases.
  using Complex = std::complex<Real>;
  using WignerAll = Wigner<Real, All, Norm>;
  using WignerNonNegative = Wigner<Real, NonNegative, Norm>;
  using WignerAllPointer = std::unique_ptr<WignerAll>;
  using WignerNonNegativePointer = std::unique_ptr<WignerNonNegative>;
  using WignerAllVector = std::vector<WignerAllPointer>;
  using WignerNonNegativeVector = std::vector<WignerNonNegativePointer>;

 public:
  Plan() = default;

  Plan(int lMax, int nMax, FFTWpp::PlanFlag flag = FFTWpp::Measure)
      : _lMax{lMax},
        _nMax{nMax},
        _quad{GaussQuad::GaussLegendreQuadrature1D<Real>(_lMax + 1)} {
    // Check the inputs.
    assert(_lMax >= 0);
    assert(_nMax >= 0 && _nMax <= _lMax);

    if constexpr (std::same_as<Type, R2R>) {
      assert(MaxUpperIndex() == 0);
    }

    // Transform the quadrature points.
    _quad.Transform([](auto x) { return std::acos(-x); },
                    [](auto x) { return 1; });

    // Compute the Wigner d-functions.
    if constexpr (std::same_as<Type, C2C>) {
      _d = WignerAllVector(2 * _nMax + 1);
      for (auto n = -_nMax; n <= _nMax; n++) {
        std::get<WignerAllVector>(_d)[n + _nMax] = std::make_unique<WignerAll>(
            MaxDegree(), MaxDegree(), n, _quad.Points());
      }
    }

    if constexpr (std::same_as<Type, R2C> or std::same_as<Type, R2R>) {
      _d = WignerNonNegativeVector(_nMax + 1);
      for (auto n = 0; n <= _nMax; n++) {
        std::get<WignerNonNegativeVector>(_d)[n] =
            std::make_unique<WignerNonNegative>(MaxDegree(), MaxDegree(), n,
                                                _quad.Points());
      }
    }

    // Generate wisdom for FFTW.
    if constexpr (std::same_as<Type, C2C>) {
      _inLayout = FFTWpp::DataLayout(1, std::vector{NumberOfLongitudes()}, 1,
                                     std::vector{NumberOfLongitudes()}, 1, 1);
      _outLayout = _inLayout;
      FFTWpp::GenerateWisdom<Complex, Complex, true>(_inLayout, _outLayout,
                                                     flag);
    }

    if constexpr (std::same_as<Type, R2C> or std::same_as<Type, R2R>) {
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

  // Return the maximum upper index.
  auto MaxUpperIndex() const { return _nMax; }

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
    auto dphi = DeltaLongitude();
    return std::ranges::views::iota(0, 2 * _lMax) |
           std::ranges::views::transform([dphi](auto i) { return i * dphi; });
  }

  // View over colatitude index.
  auto IndexCoLatitudes() const {
    return std::ranges::views::iota(0, _lMax + 1);
  }

  // View over longitude index
  auto IndexLongitudes() const {
    return std::ranges::views::iota(0, 2 * _lMax);
  }

  // Integrate function.
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

    // Forward C2C transformation.
  template <std::ranges::random_access_range RangeIn,
            std::ranges::random_access_range RangeOut>
  void Execute(int n, RangeIn& in,
               RangeOut& out) requires std::same_as<Type, C2C> {
    auto nPhi = NumberOfLongitudes();

    auto tmp = std::vector<Complex>(nPhi);
    auto tmpView = FFTWpp::DataView(tmp.begin(), tmp.end(), _outLayout);

    auto start = in.begin();
    for (auto ith : IndexCoLatitudes()) {
      // Set the end of the current longitudes.
      auto finish = std::next(start, nPhi);

      // Make views to the FFT data.
      auto inView = FFTWpp::DataView(start, finish, _inLayout);

      // Form the FFT plan.
      auto FFTPlan =
          FFTWpp::Plan(inView, tmpView, FFTWpp::WisdomOnly, FFTWpp::Forward);

      // Perform the FFT
      FFTPlan.Execute();

      // Loop over spherical harmonic degree.
      for (int l = 0; l < _lMax; l++) {
        // Loop over negative orders.
        for (int m = -l; m < 0; m++) {
        }

        // Loop over positive orders.
        for (int m = 0; m <= l; m++) {
        }
      }

      // Set the start for the next longitudes
      auto start = std::next(finish);
    }
  }

 private:
  // Store maximum orders and upper indices.
  int _lMax;
  int _nMax;

  // View to the longitudes.

  // Store the data layouts
  FFTWpp::DataLayout _inLayout;
  FFTWpp::DataLayout _outLayout;

  // Store quadrature points and weights.
  GaussQuad::Quadrature1D<Real> _quad;

  // Vector of pointers to the Wigner values for difference upper indices.
  std::variant<WignerAllVector, WignerNonNegativeVector> _d;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_PLAN_GUARD_H
