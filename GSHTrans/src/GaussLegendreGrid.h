#ifndef GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
#define GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H

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

template <RealFloatingPoint Real, TransformType Type, Normalisation Norm>
class GaussLegendreGrid {
  using Int = std::ptrdiff_t;
  using Complex = std::complex<Real>;
  using MRange = Type::IndexRange;
  using NRange = Type::IndexRange;
  using WignerType = Wigner<Real, MRange, NRange, Norm>;
  using QuadType = GaussQuad::Quadrature1D<Real>;

 public:
  GaussLegendreGrid() = default;

  GaussLegendreGrid(int lMax, int nMax, FFTWpp::PlanFlag flag = FFTWpp::Measure)
      : _lMax{lMax}, _nMax{nMax} {
    assert(MaxDegree() >= 0);
    assert(MaxUpperIndex() >= 0 && MaxUpperIndex() <= MaxDegree());
    _quad = GaussQuad::GaussLegendreQuadrature1D<Real>(_lMax + 1);
    _quad.Transform([](auto x) { return std::acos(-x); },
                    [](auto x) { return 1; });
    _wigner = std::make_shared<WignerType>(_lMax, _lMax, _nMax, _quad.Points());
    _inLayout = FFTWpp::DataLayout(1, std::vector{NumberOfLongitudes()}, 1,
                                   std::vector{NumberOfLongitudes()}, 1, 1);
    _outLayout = FFTWpp::DataLayout(1, std::vector{FFTWorkDimension()}, 1,
                                    std::vector{FFTWorkDimension()}, 1, 1);
    if constexpr (std::same_as<Type, C2C>) {
      FFTWpp::GenerateWisdom<Complex, Complex, true>(_inLayout, _outLayout,
                                                     flag);
    } else {
      FFTWpp::GenerateWisdom<Real, Complex, true>(_inLayout, _outLayout, flag);
    }
  }

  GaussLegendreGrid(const GaussLegendreGrid&) = default;

  GaussLegendreGrid(GaussLegendreGrid&&) = default;

  GaussLegendreGrid& operator=(const GaussLegendreGrid&) = default;

  GaussLegendreGrid& operator=(GaussLegendreGrid&&) = default;

  auto MaxDegree() const { return _lMax; }
  auto MaxUpperIndex() const { return _nMax; }

  auto CoLatitudes() const { return std::ranges::views::all(_quad.Points()); }
  auto NumberOfCoLatitudes() const { return CoLatitudes().size(); }
  auto CoLatitudeIndices() const { return _wigner->AngleIndices(); }

  auto LongitudeSpacing() const {
    return 2 * std::numbers::pi_v<Real> / static_cast<Real>(2 * _lMax);
  }
  auto LongitudeIndices() const {
    return std::ranges::views::iota(0, 2 * _lMax);
  }
  auto Longitudes() const {
    auto dPhi = LongitudeSpacing();
    return LongitudeIndices() |
           std::ranges::views::transform([dPhi](auto i) { return i * dPhi; });
  }
  auto NumberOfLongitudes() const { return Longitudes().size(); }

  template <typename Function>
  requires requires(Function f, Real theta, Real phi, Real w) {
    std::invocable<Function, Real, Real>;
    {f(theta, phi) * w};
  }
  auto Integrate(Function f) {
    using FunctionValue = decltype(f(0, 0));
    auto thetaIntegrand = [this, &f](auto theta) {
      return LongitudeSpacing() *
             std::accumulate(Longitudes().begin(), Longitudes().end(),
                             FunctionValue{0},
                             [&theta, &f](auto acc, auto phi) {
                               return acc + f(theta, phi);
                             });
    };
    return _quad.Integrate(thetaIntegrand);
  };

  template <RealOrComplexFloatingPointIterator InIterator,
            ComplexFloatingPointIterator OutIterator>
  void ForwardTransformation(Int lMax, Int n, InIterator in,
                             OutIterator out) const {
    const auto nPhi = NumberOfLongitudes();
    const auto scaleFactor = static_cast<Real>(2) * std::numbers::pi_v<Real> /
                             static_cast<Real>(nPhi);
    auto work = FFTWork();
    auto outView = FFTWpp::DataView(work.begin(), work.end(), _outLayout);

#pragma omp parallel for
    for (auto iTheta : CoLatitudeIndices()) {
      auto offset = iTheta * nPhi;
      auto inStart = std::next(in, offset);
      auto inFinish = std::next(inStart, nPhi);
      auto inView = FFTWpp::DataView(inStart, inFinish, _inLayout);
      auto plan =
          FFTWpp::Plan(inView, outView, FFTWpp::Forward, FFTWpp::WisdomOnly);
      plan.Execute();

      const auto d = _wigner->operator()(n)(iTheta);
      const auto w = _quad.X(iTheta) * scaleFactor;

      auto outIter = out;
      auto wigIter = d.cbegin();
      const auto degrees =
          d.Degrees() |
          std::ranges::views::filter([lMax](auto l) { return l <= lMax; });
      for (auto l : degrees) {
        auto dl = d(l);
        if constexpr (std::same_as<Type, C2C>) {
          auto workIter = std::prev(work.end(), l);
          for (auto m : dl.NegativeOrders()) {
            *outIter++ += *wigIter++ * *work++ * w;
          }
        }
        {
          auto workIter = work.begin();
          for (auto m : dl.NonNegativeOrders()) {
            *outIter++ += *wigIter++ * *work++ * w;
          }
        }
      }
    }
  }

  template <ComplexFloatingPointIterator InIterator,
            RealOrComplexFloatingPointIterator OutIterator>
  void InverseTransformation(Int lMax, Int n, InIterator in, OutIterator out) {
    const auto nPhi = NumberOfLongitudes();
    auto work = FFTWork();
    auto inView = FFTWpp::DataView(work.begin(), work.end(), _outLayout);
    for (auto iTheta : CoLatitudeIndices()) {
      for (auto& x : work) x = 0;
      auto d = _wigner->operator()(n)(iTheta);
      auto inIter = in;
      auto wigIter = d.begin();
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });
      for (auto l : degrees) {
        auto dl = d(l);
        if constexpr (std::same_as<Type, C2C>) {
          auto workIter = std::prev(work.end(), l);
          for (auto m : dl.NegativeOrders()) {
            *workIter++ += *inIter++ * *wigIter++;
          }
        }
        {
          auto workIter = work.begin();
          for (auto m : dl.NonNegativeOrders()) {
            *workIter++ += *inIter++ * *wigIter++;
          }
        }
      }
      auto offset = iTheta * nPhi;
      auto outStart = std::next(out, offset);
      auto outFinish = std::next(outStart, nPhi);
      auto outView = FFTWpp::DataView(outStart, outFinish, _inLayout);
      auto plan =
          FFTWpp::Plan(outView, outView, FFTWpp::Backward, FFTWpp::WisdomOnly);
      plan.Execute();
    }
  }

 private:
  int _lMax;
  int _nMax;

  QuadType _quad;

  std::shared_ptr<WignerType> _wigner;

  FFTWpp::DataLayout _inLayout;
  FFTWpp::DataLayout _outLayout;

  auto FFTWorkDimension() const {
    if constexpr (std::same_as<Type, C2C>) {
      return 2 * _lMax;
    } else {
      return _lMax + 1;
    }
  }

  auto FFTWork() const requires std::same_as<Type, C2C> {
    return FFTWpp::vector<Complex>(2 * _lMax);
  }

  auto FFTWork() const requires std::same_as<Type, R2C> {
    return FFTWpp::vector<Complex>(_lMax + 1);
  }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
