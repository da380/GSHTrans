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
#include <random>
#include <ranges>
#include <vector>

#include "Concepts.h"
#include "Indexing.h"
#include "Wigner.h"

namespace GSHTrans {

template <RealFloatingPoint Real, OrderIndexRange MRange, IndexRange NRange>
class GaussLegendreGrid {
 private:
  using Int = std::ptrdiff_t;
  using Complex = std::complex<Real>;

  using WignerType = Wigner<Real, MRange, NRange, Ortho>;
  using QuadType = GaussQuad::Quadrature1D<Real>;

 public:
  using real_type = Real;
  using complex_type = Complex;
  using MRange_type = MRange;
  using NRange_type = NRange;

  // Constructors.
  GaussLegendreGrid() = default;

  GaussLegendreGrid(int lMax, int nMax, FFTWpp::PlanFlag flag = FFTWpp::Measure)
      : _lMax{lMax}, _nMax{nMax} {
    // Check the inputs.
    assert(MaxDegree() >= 0);
    assert(MaxUpperIndex() <= MaxDegree());
    assert(std::abs(MinUpperIndex()) <= MaxDegree());

    // Get the quadrature points.
    _quad = GaussQuad::GaussLegendreQuadrature1D<Real>(_lMax + 1);
    _quad.Transform([](auto x) { return std::acos(-x); },
                    [](auto x) { return 1; });

    //  Get the Winger values.
    _wigner = std::make_shared<WignerType>(_lMax, _lMax, _nMax, _quad.Points());

    if (_lMax > 0) {
      // Generate wisdom for FFTs.
      auto in = FFTWpp::DataLayout(1, std::vector{NumberOfLongitudes()}, 1,
                                   std::vector{NumberOfLongitudes()}, 1, 1);
      {
        // Real to complex case.
        auto out = FFTWpp::DataLayout(1, std::vector{_lMax + 1}, 1,
                                      std::vector{_lMax + 1}, 1, 1);
        FFTWpp::GenerateWisdom<Real, Complex, true>(in, out, flag);
      }
      {
        // Complex to complex case.
        auto out = FFTWpp::DataLayout(1, std::vector{2 * _lMax}, 1,
                                      std::vector{2 * _lMax}, 1, 1);
        FFTWpp::GenerateWisdom<Complex, Complex, true>(in, out, flag);
      }
    }
  }

  GaussLegendreGrid(const GaussLegendreGrid&) = default;

  GaussLegendreGrid(GaussLegendreGrid&&) = default;

  GaussLegendreGrid& operator=(const GaussLegendreGrid&) = default;

  GaussLegendreGrid& operator=(GaussLegendreGrid&&) = default;

  // Grid information.
  auto MaxDegree() const { return _lMax; }
  auto MaxUpperIndex() const { return _nMax; }
  auto MinUpperIndex() const {
    if constexpr (std::same_as<NRange, All>) {
      return -_nMax;
    }
    if constexpr (std::same_as<NRange, NonNegative>) {
      return 0;
    }
    if constexpr (std::same_as<NRange, Single>) {
      return _nMax;
    }
  }
  auto UpperIndices() const {
    return std::ranges::views::iota(MinUpperIndex(), MaxUpperIndex() + 1);
  }

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

  auto Points() const {
    return std::ranges::views::cartesian_product(CoLatitudes(), Longitudes());
  }

  auto PointIndices() const {
    return std::ranges::views::cartesian_product(CoLatitudeIndices(),
                                                 LongitudeIndices());
  }

  auto ComponentSize() const {
    return NumberOfCoLatitudes() * NumberOfLongitudes();
  }

  auto RealCoefficientSize(Int lMax, Int n) const {
    return GSHIndices<NonNegative>(lMax, lMax, n).size();
  }

  auto ComplexCoefficientSize(Int lMax, Int n) const {
    return GSHIndices<All>(lMax, lMax, n).size();
  }

  auto RealCoefficientSize(Int n) const {
    return GSHIndices<NonNegative>(_lMax, _lMax, n).size();
  }

  auto ComplexCoefficientSize(Int n) const {
    return GSHIndices<All>(_lMax, _lMax, n).size();
  }

  // Generate random coefficient values within a given range.
  template <RealOrComplexFloatingPointRange Range,
            typename Distribution = decltype(std::normal_distribution<Real>())>
  void RandomComplexCoefficient(
      Int lMax, Int n, Range& range,
      Distribution dist = std::normal_distribution<Real>()) const {
    assert(range.size() == ComplexCoefficientSize(lMax, n));
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::ranges::generate(range, [&gen, &dist]() {
      return Complex{dist(gen), dist(gen)};
    });
    auto view = GSHView<Complex, All>(lMax, lMax, n, range.begin());
    if (lMax == _lMax) view(lMax)(lMax) = 0;
  }

  // Generate random coefficient values within a given range.
  template <RealOrComplexFloatingPointRange Range,
            typename Distribution = decltype(std::normal_distribution<Real>())>
  void RandomRealCoefficient(
      Int lMax, Int n, Range& range,
      Distribution dist = std::normal_distribution<Real>()) const {
    assert(range.size() == RealCoefficientSize(lMax, n));
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::ranges::generate(range, [&gen, &dist]() {
      return Complex{dist(gen), dist(gen)};
    });
    auto view = GSHView<Complex, NonNegative>(lMax, lMax, n, range.begin());
    for (auto l : view.Degrees()) {
      view(l)(0).imag(0);
    }
    if (lMax == _lMax) view(lMax)(lMax).imag(0);
  }

  // Integration a function over the grid.
  template <typename Function>
  requires ScalarFunction2D<Function, Real, Real> or
           ScalarFunction2D<Function, Real, Complex>
  auto Integrate(Function f) const {
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
    auto summand = [this, nPhi, dPhi, range](auto iTheta) {
      auto start = std::next(range.begin(), iTheta * nPhi);
      auto finish = std::next(start, nPhi);
      return dPhi * std::accumulate(start, finish, Scalar{0});
    };
    return std::inner_product(CoLatitudeIndices().begin(),
                              CoLatitudeIndices().end(),
                              _quad.Weights().begin(), Scalar{0}, std::plus<>(),
                              [this, nPhi, dPhi, range, summand](
                                  auto i, auto w) { return summand(i) * w; });
  }

  //-----------------------------------------------------//
  //                Forward transformation               //
  //-----------------------------------------------------//
  template <RealOrComplexFloatingPointRange InRange,
            ComplexFloatingPointRange OutRange>
  requires(std::same_as<MRange, All> and ComplexFloatingPointRange<InRange>) or
          RealFloatingPointRange<InRange>
  void ForwardTransformation(Int lMax, Int n, InRange&& in,
                             OutRange& out) const {
    // Get scalar type for field.
    using Scalar = std::ranges::range_value_t<InRange>;

    // Check upper index is possible.
    assert(
        std::ranges::any_of(UpperIndices(), [n](auto np) { return n == np; }));

    // Check dimensions of ranges.
    assert(in.size() == NumberOfCoLatitudes() * NumberOfLongitudes());
    if constexpr (RealFloatingPoint<Scalar>) {
      assert(out.size() == GSHIndices<NonNegative>(lMax, lMax, n).size());
    } else {
      assert(out.size() == GSHIndices<All>(lMax, lMax, n).size());
    }

    // Deal with lMax = 0
    if (lMax == 0) {
      out[0] = in[0] * std::numbers::inv_sqrtpi_v<Real> / static_cast<Real>(2);
      return;
    }

    // Pre compute some constants.
    const auto nPhi = NumberOfLongitudes();
    const auto scaleFactor = static_cast<Real>(2) * std::numbers::pi_v<Real> /
                             static_cast<Real>(nPhi);

    // Allocate work array for FFTs and form its FFTWpp view.
    auto workSize = WorkSize<Scalar>();
    auto work = FFTWpp::vector<Complex>(workSize);
    auto inLayout = FFTWpp::DataLayout(1, std::vector{NumberOfLongitudes()}, 1,
                                       std::vector{NumberOfLongitudes()}, 1, 1);
    auto outLayout = FFTWpp::DataLayout(1, std::vector{workSize}, 1,
                                        std::vector{workSize}, 1, 1);
    auto outView = FFTWpp::DataView(work.begin(), work.end(), outLayout);

    // Loop over the colatitudes.
    for (auto iTheta : CoLatitudeIndices()) {
      // Form FFT of current data slice.
      auto offset = iTheta * nPhi;
      auto inStart = std::next(in.begin(), offset);
      auto inFinish = std::next(inStart, nPhi);
      auto inView = FFTWpp::DataView(inStart, inFinish, inLayout);
      auto plan =
          FFTWpp::Plan(inView, outView, FFTWpp::WisdomOnly, FFTWpp::Forward);
      plan.Execute();

      // Get the Wigner values and quadrature weight.
      auto d = _wigner->operator()(n)(iTheta);
      auto w = _quad.W(iTheta) * scaleFactor;

      // Loop over the spherical harmonic coefficients
      auto outIter = out.begin();
      auto wigIter = d.cbegin();
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });

      for (auto l : degrees) {
        auto dl = d(l);

        if constexpr (ComplexFloatingPoint<Scalar>) {
          auto workIter = std::prev(work.end(), l);
          for (auto m : dl.NegativeOrders()) {
            *outIter++ += *wigIter++ * *workIter++ * w;
          }
          workIter = work.begin();
          for (auto m : dl.NonNegativeOrders()) {
            *outIter++ += *wigIter++ * *workIter++ * w;
          }
        } else {
          auto workIter = work.begin();
          if constexpr (std::same_as<MRange, All>) {
            std::advance(wigIter, l);
          }
          for (auto m : dl.NonNegativeOrders()) {
            *outIter++ += *wigIter++ * *workIter++ * w;
          }
        }
      }

      if constexpr (ComplexFloatingPoint<Scalar>) {
        // Zero the (_lMax,_lMax) coefficient.
        if (lMax == _lMax) {
          auto view = std::ranges::views::reverse(out);
          view[0] = 0;
        }
      }
    }
  }

  //------------------------------------------------//
  //            Inverse  transformation             //
  //------------------------------------------------//
  template <ComplexFloatingPointRange InRange,
            RealOrComplexFloatingPointRange OutRange>
  requires(std::same_as<MRange, All> and ComplexFloatingPointRange<OutRange>) or
          RealFloatingPointRange<OutRange>
  void InverseTransformation(Int lMax, Int n, InRange&& in,
                             OutRange& out) const {
    // Get scalar type for field.
    using Scalar = std::ranges::range_value_t<OutRange>;

    // Check upper index is possible.
    assert(
        std::ranges::any_of(UpperIndices(), [n](auto np) { return n == np; }));

    // Check dimensions of ranges.
    if constexpr (RealFloatingPoint<Scalar>) {
      assert(in.size() == GSHIndices<NonNegative>(lMax, lMax, n).size());
    } else {
      assert(in.size() == GSHIndices<All>(lMax, lMax, n).size());
    }
    assert(out.size() == NumberOfCoLatitudes() * NumberOfLongitudes());

    // Deal with lMax = 0
    if (lMax == 0) {
      out[0] = in[0] * static_cast<Real>(2) / std::numbers::inv_sqrtpi_v<Real>;
      return;
    }

    // Precompute constants
    const auto nPhi = NumberOfLongitudes();

    // Set up for FFTs.
    auto workSize = WorkSize<Scalar>();
    auto work = FFTWpp::vector<Complex>(workSize);
    auto inLayout = FFTWpp::DataLayout(1, std::vector{workSize}, 1,
                                       std::vector{workSize}, 1, 1);
    auto outLayout =
        FFTWpp::DataLayout(1, std::vector{NumberOfLongitudes()}, 1,
                           std::vector{NumberOfLongitudes()}, 1, 1);
    auto inView = FFTWpp::DataView(work.begin(), work.end(), inLayout);

    // Loop over the colatitudes.
    for (auto iTheta : CoLatitudeIndices()) {
      std::ranges::for_each(work, [](auto& x) { return x = 0; });

      // Get the Wigner values.
      auto d = (*_wigner)(n)(iTheta);

      // Loop over the coefficients.
      auto inIter = in.begin();
      auto wigIter = d.begin();
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });
      for (auto l : degrees) {
        auto dl = d(l);
        if constexpr (ComplexFloatingPoint<Scalar>) {
          auto workIter = std::prev(work.end(), l);
          for (auto m : dl.NegativeOrders()) {
            *workIter++ += *inIter++ * *wigIter++;
          }
          workIter = work.begin();
          for (auto m : dl.NonNegativeOrders()) {
            *workIter++ += *inIter++ * *wigIter++;
          }
        } else {
          auto workIter = work.begin();
          if constexpr (std::same_as<MRange, All>) {
            std::advance(wigIter, l);
          }
          for (auto m : dl.NonNegativeOrders()) {
            *workIter++ += *inIter++ * *wigIter++;
          }
        }
      }

      // Perform FFT to recover field at the colatitude.
      auto offset = iTheta * nPhi;
      auto outStart = std::next(out.begin(), offset);
      auto outFinish = std::next(outStart, nPhi);
      auto outView = FFTWpp::DataView(outStart, outFinish, outLayout);
      auto plan =
          FFTWpp::Plan(inView, outView, FFTWpp::WisdomOnly, FFTWpp::Backward);
      plan.Execute();
    }
  }

 private:
  Int _lMax;
  Int _nMax;
  QuadType _quad;
  std::shared_ptr<WignerType> _wigner;

  template <RealOrComplexFloatingPoint Scalar>
  auto WorkSize() const {
    if constexpr (RealFloatingPoint<Scalar>) {
      return _lMax + 1;
    } else {
      return 2 * _lMax;
    }
  }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
