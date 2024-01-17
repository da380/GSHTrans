#ifndef GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
#define GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H

#include <FFTWpp/All>
#include <GaussQuad/All>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <iterator>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include "Concepts.h"
#include "Indexing.h"
#include "Wigner.h"

namespace GSHTrans {

template <RealFloatingPoint Real, TransformType Type>
class GaussLegendreGrid {
 public:
  using real_type = Real;
  using scalar_type = typename Type::Scalar<Real>;

 private:
  using Int = std::ptrdiff_t;
  using Complex = std::complex<Real>;
  using MRange = Type::IndexRange;
  using NRange = Type::IndexRange;
  using WignerType = Wigner<Real, MRange, NRange, Ortho>;
  using QuadType = GaussQuad::Quadrature1D<Real>;

 public:
  // Constructors.
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

  // Grid information.
  Int MaxDegree() const { return _lMax; }
  Int MaxUpperIndex() const { return _nMax; }

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

  // Integration a function over the grid.
  template <typename Function>
  requires ScalarFunction2D<Function, Real>
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

  // Integration of discretised function over the grid.
  template <RealOrComplexFloatingPointRange Range>
  auto Integrate(Range&& range) const {
    using Scalar = std::ranges::range_value_t<Range>;
    auto nPhi = NumberOfLongitudes();
    auto dPhi = LongitudeSpacing();
    auto summand = [this,nPhi,dPhi,range](auto iTheta) {
      auto start = std::next(range.begin(),iTheta*nPhi);
      auto finish = std::next(start,nPhi);
      return dPhi*std::accumulate(start,finish,Scalar{0});
    };
    return std::inner_product(CoLatitudeIndices().begin(), CoLatitudeIndices().end(),
			      _quad.Weights().begin(),Scalar{0}, std::plus<>(),
			      [this,nPhi,dPhi,range,summand](auto i, auto w){ return summand(i) * w ;});
  }

  // Fast transformations.
  template <RealOrComplexFloatingPointIterator InIterator,
            ComplexFloatingPointIterator OutIterator>
  void ForwardTransformation(Int lMax, Int n, InIterator in,
                             OutIterator out) const {
    assert(std::abs(n) <= _nMax);

    const auto nPhi = NumberOfLongitudes();
    const auto scaleFactor = static_cast<Real>(2) * std::numbers::pi_v<Real> /
                             static_cast<Real>(nPhi);
    auto work = FFTWork();
    auto outView = FFTWpp::DataView(work.begin(), work.end(), _outLayout);

    for (auto iTheta : CoLatitudeIndices()) {
      auto offset = iTheta * nPhi;
      auto inStart = std::next(in, offset);
      auto inFinish = std::next(inStart, nPhi);
      auto inView = FFTWpp::DataView(inStart, inFinish, _inLayout);
      auto plan =
          FFTWpp::Plan(inView, outView, FFTWpp::WisdomOnly, FFTWpp::Forward);
      plan.Execute();

      auto d = (_wigner->operator()(n))(iTheta);
      auto w = _quad.W(iTheta) * scaleFactor;

      auto outIter = out;
      auto wigIter = d.cbegin();
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });

      for (auto l : degrees) {
        auto dl = d(l);

        if constexpr (std::same_as<Type, C2C>) {
          auto workIter = std::prev(work.end(), l);
          for (auto m : dl.NegativeOrders()) {
            *outIter++ += *wigIter++ * *workIter++ * w;
          }
        }

        {
          auto workIter = work.begin();
          for (auto m : dl.NonNegativeOrders()) {
            *outIter++ += *wigIter++ * *workIter++ * w;
          }
        }
      }
    }
  }

  template <RealOrComplexFloatingPointIterator InIterator,
            ComplexFloatingPointIterator OutIterator>
  void ForwardTransformation(Int n, InIterator in, OutIterator out) const {
    ForwardTransformation(_lMax, n, in, out);
  }

  template <ComplexFloatingPointIterator InIterator,
            RealOrComplexFloatingPointIterator OutIterator>
  void InverseTransformation(Int lMax, Int n, InIterator in,
                             OutIterator out) const {
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
          FFTWpp::Plan(inView, outView, FFTWpp::WisdomOnly, FFTWpp::Backward);
      plan.Execute();
    }
  }

  // Return vector to store function values.
  template <RealOrComplexFloatingPoint Scalar>
  auto FunctionVector() const {
    return FFTWpp::vector<Scalar>(NumberOfCoLatitudes() * NumberOfLongitudes());
  }

  // Interpolation of function onto the grid.
  template <typename Function>
  requires ScalarFunction2D<Function, Real>
  auto Interpolate(Function f) const {
    using FunctionValue = decltype(f(0, 0));
    auto vec = FunctionVector<FunctionValue>();
    auto iter = vec.begin();
    for (auto theta : CoLatitudes()) {
      for (auto phi : Longitudes()) {
        *iter++ = f(theta, phi);
      }
    }
    return vec;
  }

  // Return vector to store coefficients.
  auto CoefficientVector(Int n) const {
    auto size = GSHIndices<MRange>(_lMax, _lMax, n).size();
    return FFTWpp::vector<Complex>(size);
  }

 private:
  Int _lMax;
  Int _nMax;

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
