#ifndef GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
#define GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H

#include <FFTWpp/All>
#include <GaussQuad/All>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <memory>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include "Concepts.h"
#include "GridBase.h"
#include "Indexing.h"
#include "Wigner.h"

namespace GSHTrans {

template <RealFloatingPoint Real, OrderIndexRange MRange, IndexRange NRange>
class GaussLegendreGrid
    : public GridBase<GaussLegendreGrid<Real, MRange, NRange>> {
 private:
  using Int = std::ptrdiff_t;
  using Complex = std::complex<Real>;

  using WignerType = Wigner<Real, Ortho, MRange, NRange, Multiple, ColumnMajor>;
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
    assert(std::abs(this->MinUpperIndex()) <= MaxDegree());

    // Get the quadrature points.
    _quadPointer = std::make_shared<QuadType>(
        GaussQuad::LegendrePolynomial<Real>{}.GaussQuadrature(_lMax + 1));

    _quadPointer->Transform([](auto x) { return std::acos(-x); },
                            [](auto x) { return 1; });

    //  Get the Winger values.
    _wignerPointer = std::make_shared<WignerType>(_lMax, _lMax, _nMax,
                                                  _quadPointer->Points());

    if (_lMax > 0) {
      // Generate wisdom for FFTs.
      auto in =
          FFTWpp::DataLayout(1, std::vector{this->NumberOfLongitudes()}, 1,
                             std::vector{this->NumberOfLongitudes()}, 1, 1);
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

  auto MaxDegree() const { return _lMax; }
  auto MaxUpperIndex() const { return _nMax; }

  auto CoLatitudes() const {
    return std::ranges::views::all(_quadPointer->Points());
  }
  auto CoLatitudeWeights() const {
    return std::ranges::views::all(_quadPointer->Weights());
  }

  auto LongitudeSpacing() const {
    return 2 * std::numbers::pi_v<Real> / static_cast<Real>(2 * _lMax);
  }
  auto Longitudes() const {
    auto dPhi = LongitudeSpacing();
    return std::ranges::views::iota(0, 2 * _lMax) |
           std::ranges::views::transform([dPhi](auto i) { return i * dPhi; });
  }
  auto LongitudeWeights() const {
    auto dPhi = LongitudeSpacing();
    return std::ranges::views::repeat(dPhi, 2 * _lMax);
  }

  //-----------------------------------------------------//
  //                Forward transformation               //
  //-----------------------------------------------------//
  template <std::ranges::range InRange, std::ranges::range OutRange>
  requires requires() {
    requires(std::same_as<MRange, All> and
             ComplexFloatingPoint<std::ranges::range_value_t<InRange>>) or
                RealFloatingPoint<std::ranges::range_value_t<InRange>>;
    requires std::ranges::input_range<InRange>;
    requires std::same_as<RemoveComplex<std::ranges::range_value_t<InRange>>,
                          Real>;
    requires ComplexFloatingPoint<std::ranges::range_value_t<OutRange>>;
    requires std::ranges::output_range<OutRange,
                                       std::ranges::range_value_t<OutRange>>;
    requires std::same_as<RemoveComplex<std::ranges::range_value_t<OutRange>>,
                          Real>;
  }
  void ForwardTransformation(Int lMax, Int n, InRange&& in, OutRange& out) {
    // Get scalar type for field.
    using Scalar = std::ranges::range_value_t<InRange>;

    // Check upper index is possible.
    assert(std::ranges::contains(this->UpperIndices(), n));

    // Check dimensions of ranges.
    assert(in.size() == this->ComponentSize());
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
    const auto nPhi = this->NumberOfLongitudes();
    const auto scaleFactor = static_cast<Real>(2) * std::numbers::pi_v<Real> /
                             static_cast<Real>(nPhi);

    // Allocate work array for FFTs set up the FFT plan.
    auto inWork = FFTWpp::vector<Scalar>(nPhi);
    auto inLayout =
        FFTWpp::DataLayout(1, std::vector{nPhi}, 1, std::vector{nPhi}, 1, 1);
    auto inView = FFTWpp::DataView(inWork.begin(), inWork.end(), inLayout);
    auto outWorkSize = WorkSize<Scalar>();
    auto outWork = FFTWpp::vector<Complex>(outWorkSize);
    auto outLayout = FFTWpp::DataLayout(1, std::vector{outWorkSize}, 1,
                                        std::vector{outWorkSize}, 1, 1);
    auto outView = FFTWpp::DataView(outWork.begin(), outWork.end(), outLayout);
    auto plan =
        FFTWpp::Plan(inView, outView, FFTWpp::WisdomOnly, FFTWpp::Forward);

    // Loop over the colatitudes.
    for (auto iTheta : this->CoLatitudeIndices()) {
      // FFT the current data slice.
      auto offset = iTheta * nPhi;
      auto inStart = std::next(in.begin(), offset);
      auto inFinish = std::next(inStart, nPhi);
      if constexpr (std::ranges::output_range<InRange, Scalar>) {
        auto inView = FFTWpp::DataView(inStart, inFinish, inLayout);
        plan.Execute(inView, outView);
      } else {
        std::copy(inStart, inFinish, inWork.begin());
        plan.Execute();
      }

      // Get the Wigner values and quadrature weight.
      auto d = _wignerPointer->operator()(n, iTheta);
      auto w = _quadPointer->W(iTheta) * scaleFactor;

      // Loop over the spherical harmonic coefficients
      auto outIter = out.begin();
      auto wigIter = d.cbegin();
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });

      for (auto l : degrees) {
        auto dl = d(l);

        if constexpr (ComplexFloatingPoint<Scalar>) {
          auto workIter = std::prev(outWork.end(), l);
          for (auto m : dl.NegativeOrders()) {
            *outIter++ += *wigIter++ * *workIter++ * w;
          }
          workIter = outWork.begin();
          for (auto m : dl.NonNegativeOrders()) {
            *outIter++ += *wigIter++ * *workIter++ * w;
          }
        } else {
          auto workIter = outWork.begin();
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
  template <std::ranges::range InRange, std::ranges::range OutRange>
  requires requires() {
    requires ComplexFloatingPoint<std::ranges::range_value_t<InRange>>;
    requires std::ranges::input_range<InRange>;
    requires std::same_as<RemoveComplex<std::ranges::range_value_t<InRange>>,
                          Real>;
    requires(std::same_as<MRange, All> and
             ComplexFloatingPoint<std::ranges::range_value_t<OutRange>>) or
                RealFloatingPoint<std::ranges::range_value_t<OutRange>>;
    requires std::ranges::output_range<OutRange,
                                       std::ranges::range_value_t<OutRange>>;
    requires std::same_as<RemoveComplex<std::ranges::range_value_t<OutRange>>,
                          Real>;
  }
  void InverseTransformation(Int lMax, Int n, InRange&& in, OutRange& out) {
    // Get scalar type for field.
    using Scalar = std::ranges::range_value_t<OutRange>;

    // Check upper index is possible.
    assert(std::ranges::contains(this->UpperIndices(), n));

    // Check dimensions of ranges.
    if constexpr (RealFloatingPoint<Scalar>) {
      assert(in.size() == GSHIndices<NonNegative>(lMax, lMax, n).size());
    } else {
      assert(in.size() == GSHIndices<All>(lMax, lMax, n).size());
    }
    assert(out.size() == this->ComponentSize());

    // Deal with lMax = 0
    if (lMax == 0) {
      if constexpr (RealFloatingPoint<Scalar>) {
        out[0] = std::real(in[0]) * static_cast<Real>(2) /
                 std::numbers::inv_sqrtpi_v<Real>;
      } else {
        out[0] =
            in[0] * static_cast<Real>(2) / std::numbers::inv_sqrtpi_v<Real>;
      }
      return;
    }

    // Precompute constants
    const auto nPhi = this->NumberOfLongitudes();

    // Set up for FFTs.
    auto inWorkSize = WorkSize<Scalar>();
    auto inWork = FFTWpp::vector<Complex>(inWorkSize);
    auto inLayout = FFTWpp::DataLayout(1, std::vector{inWorkSize}, 1,
                                       std::vector{inWorkSize}, 1, 1);
    auto inView = FFTWpp::DataView(inWork.begin(), inWork.end(), inLayout);
    auto outWork = FFTWpp::vector<Scalar>(nPhi);
    auto outLayout =
        FFTWpp::DataLayout(1, std::vector{nPhi}, 1, std::vector{nPhi}, 1, 1);
    auto outView = FFTWpp::DataView(outWork.begin(), outWork.end(), outLayout);
    auto plan =
        FFTWpp::Plan(inView, outView, FFTWpp::WisdomOnly, FFTWpp::Backward);

    // Loop over the colatitudes.
    for (auto iTheta : this->CoLatitudeIndices()) {
      std::ranges::for_each(inWork, [](auto& x) { return x = 0; });

      // Get the Wigner values.
      auto d = _wignerPointer->operator()(n, iTheta);

      // Loop over the coefficients.
      auto inIter = in.begin();
      auto wigIter = d.begin();
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });
      for (auto l : degrees) {
        auto dl = d(l);
        if constexpr (ComplexFloatingPoint<Scalar>) {
          auto workIter = std::prev(inWork.end(), l);
          for (auto m : dl.NegativeOrders()) {
            *workIter++ += *inIter++ * *wigIter++;
          }
          workIter = inWork.begin();
          for (auto m : dl.NonNegativeOrders()) {
            *workIter++ += *inIter++ * *wigIter++;
          }
        } else {
          auto workIter = inWork.begin();
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
      plan.Execute(inView, outView);
    }
  }

 private:
  Int _lMax;
  Int _nMax;

  std::shared_ptr<QuadType> _quadPointer;
  std::shared_ptr<WignerType> _wignerPointer;

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
