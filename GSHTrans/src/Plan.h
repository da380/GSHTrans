#ifndef GSH_TRANS_PLAN_GUARD_H
#define GSH_TRANS_PLAN_GUARD_H

#include <FFTWpp/All>
#include <GaussQuad/All>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <execution>
#include <iterator>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include "Concepts.h"
#include "Legendre.h"
#include "Wigner.h"

namespace GSHTrans {

template <std::floating_point Real>
class PlanComplex {
  using Complex = std::complex<Real>;
  using wigner = Wigner<Real, All>;

 public:
  PlanComplex() = default;

  PlanComplex(int lMax, int nMax, FFTWpp::PlanFlag flag = FFTWpp::Measure)
      : _lMax{lMax},
        _nMax{nMax},
        _quad{GaussQuad::GaussLegendreQuadrature1D<Real>(_lMax + 1)} {
    // Check the inputs.
    assert(_lMax >= 0);
    assert(_nMax >= 0 && _nMax <= _lMax);

    // Transform the quadrature points.
    _quad.Transform([](auto x) { return std::acos(x); },
                    [](auto x) { return 1; });

    // Compute the Wigner d-functions.
    for (auto n = -_nMax; n <= _nMax; n++) {
      _d.push_back(std::make_unique<wigner>(_lMax, _lMax, n, _quad.Points()));
    }

    // Generate wisdom for FFTW
    _layout = FFTWpp::DataLayout(1, {2 * _lMax}, _lMax + 1, {2 * _lMax}, 1,
                                 2 * _lMax);
    FFTWpp::GenerateWisdom<Complex, Complex, true>(_layout, _layout, flag);
  }

  // Return the longitude spacing.
  auto DeltaLongitude() const {
    return 2 * std::numbers::pi_v<Real> / static_cast<Real>(2 * _lMax - 1);
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

  // Return view to the longitudes.
  auto LongitudesWithoutRepeat() const {
    return Longitudes() | std::ranges::views::take(NumberOfLongitudes() - 1);
  }

  // Integrate function.
  template <typename Function>
  auto Integrate(Function f) {
    auto integral = Real{0};
    auto dphi = DeltaLongitude();
    auto phi = LongitudesWithoutRepeat();
    auto w = _quad.Weights().begin();
    for (auto theta : CoLatitudes()) {
      auto sum = std::accumulate(std::begin(phi), std::end(phi), Real{0},
                                 [&theta, &dphi, &f](auto acc, auto phi) {
                                   return acc + f(theta, phi) * dphi;
                                 });
      integral += sum * (*w++);
    }
    return integral;
  }

   private:
    // Store maximum orders and upper indices.
    int _lMax;
    int _nMax;

    // View to the longitudes.

    // Store the data layout
    FFTWpp::DataLayout _layout;

    // Store quadrature points and weights.
    GaussQuad::Quadrature1D<Real> _quad;

    // Store the Wigner values
    std::vector<std::unique_ptr<wigner>> _d;
  };

}  // namespace GSHTrans

#endif  // GSH_TRANS_PLAN_GUARD_H
