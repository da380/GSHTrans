#ifndef GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
#define GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H

#include <FFTWpp/Core>
#include <FFTWpp/Ranges>
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

template <RealFloatingPoint _Real, OrderIndexRange _MRange, IndexRange _NRange>
class GaussLegendreGrid
    : public GridBase<GaussLegendreGrid<_Real, _MRange, _NRange>> {
 public:
  // Public type aliases.
  using Int = std::ptrdiff_t;
  using Real = _Real;
  using Complex = std::complex<Real>;
  using MRange = _MRange;
  using NRange = _NRange;

 private:
  // Private type aliases.
  using WignerType =
      Wigner<Real, Ortho, _MRange, _NRange, Multiple, ColumnMajor>;
  using QuadType = GaussQuad::Quadrature1D<Real>;

 public:
  // Constructors.
  GaussLegendreGrid() = default;

  GaussLegendreGrid(int lMax, int nMax, FFTWpp::Flag flag = FFTWpp::Measure)
      : _lMax{lMax}, _nMax{nMax}, _flag{flag} {
    // Check the inputs.
    assert(MaxDegree() >= 0);
    assert(MaxUpperIndex() <= MaxDegree());
    assert(std::abs(this->MinUpperIndex()) <= MaxDegree());
    assert(_flag != FFTWpp::WisdomOnly);

    // Get the quadrature points.
    _quadPointer = std::make_shared<QuadType>(
        GaussQuad::LegendrePolynomial<Real>{}.GaussQuadrature(_lMax + 1));

    _quadPointer->Transform([](auto x) { return std::acos(-x); },
                            [](auto x) -> Real { return 1; });

    //  Get the Winger values.
    _wignerPointer = std::make_shared<WignerType>(_lMax, _lMax, _nMax,
                                                  _quadPointer->Points());

    if (_lMax > 0 && _flag != FFTWpp::Estimate) {
      // Generate wisdom for FFTs.
      auto nPhi = this->NumberOfLongitudes();
      auto in = FFTWpp::Ranges::Layout(nPhi);
      {
        // Real to complex case.
        auto out = FFTWpp::Ranges::Layout(nPhi / 2 + 1);
        FFTWpp::GenerateWisdom<Real, Complex>(in, out, flag);
      }
      {
        // Complex to complex case.
        auto out = FFTWpp::Ranges::Layout(nPhi);
        FFTWpp::GenerateWisdom<Complex, Complex>(in, out, flag);
      }
      _flag = FFTWpp::WisdomOnly;
    } else {
      _flag = FFTWpp::Estimate;
    }
  }

  GaussLegendreGrid(const GaussLegendreGrid&) = default;

  GaussLegendreGrid(GaussLegendreGrid&&) = default;

  GaussLegendreGrid& operator=(const GaussLegendreGrid&) = default;

  GaussLegendreGrid& operator=(GaussLegendreGrid&&) = default;

  //------------------------------------------------//
  //    Methods needed to inherit from GridBase     //
  //------------------------------------------------//
  auto MaxDegree() const { return _lMax; }
  auto MaxUpperIndex() const { return _nMax; }

  auto CoLatitudes() const {
    return std::ranges::views::all(_quadPointer->Points());
  }
  auto CoLatitudeWeights() const {
    return std::ranges::views::all(_quadPointer->Weights());
  }

  auto Longitudes() const {
    auto dPhi = 2 * std::numbers::pi_v<Real> / static_cast<Real>(2 * _lMax);
    return std::ranges::views::iota(0, 2 * _lMax) |
           std::ranges::views::transform([dPhi](auto i) { return i * dPhi; });
  }
  auto LongitudeWeights() const {
    auto dPhi = 2 * std::numbers::pi_v<Real> / static_cast<Real>(2 * _lMax);
    ;
    return std::ranges::views::repeat(dPhi, 2 * _lMax);
  }

  //-----------------------------------------------------//
  //                Forward transformation               //
  //-----------------------------------------------------//
  template <std::ranges::range InRange, std::ranges::range OutRange>
  requires requires() {
    requires(std::same_as<_MRange, All> and
             ComplexFloatingPoint<std::ranges::range_value_t<InRange>>) or
                RealFloatingPoint<std::ranges::range_value_t<InRange>>;
    requires std::same_as<RemoveComplex<std::ranges::range_value_t<InRange>>,
                          Real>;
    requires std::ranges::input_range<InRange>;
    requires std::same_as<std::ranges::range_value_t<OutRange>, Complex>;
    requires std::ranges::output_range<OutRange,
                                       std::ranges::range_value_t<OutRange>>;
  }
  void ForwardTransformation(Int lMax, Int n, InRange&& in,
                             OutRange& out) const {
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

    // Make the FFT plan.
    auto [inSize, outSize] = FFTWpp::DataSize<Scalar, Complex>(nPhi);
    auto inWork = FFTWpp::vector<Scalar>(inSize);
    auto outWork = FFTWpp::vector<Complex>(outSize);
    auto inView = FFTWpp::Ranges::View(inWork);
    auto outView = FFTWpp::Ranges::View(outWork);
    auto planFunction = [this](auto in, auto out) {
      if constexpr (ComplexFloatingPoint<Scalar>) {
        return FFTWpp::Ranges::Plan(in, out, _flag, FFTWpp::Forward);
      } else {
        return FFTWpp::Ranges::Plan(in, out, _flag);
      }
    };
    auto plan = planFunction(inView, outView);

    // Loop over the colatitudes.
    for (auto iTheta : this->CoLatitudeIndices()) {
      // FFT the current data slice.
      auto offset = iTheta * nPhi;
      auto inStart = std::next(in.begin(), offset);
      auto inFinish = std::next(inStart, nPhi);
      if constexpr (std::ranges::output_range<InRange, Scalar>) {
        auto inView = std::ranges::subrange(inStart, inFinish);
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
          if constexpr (std::same_as<_MRange, All>) {
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
    requires(std::same_as<_MRange, All> and
             ComplexFloatingPoint<std::ranges::range_value_t<OutRange>>) or
                RealFloatingPoint<std::ranges::range_value_t<OutRange>>;
    requires std::ranges::input_range<InRange>;
    requires std::same_as<std::ranges::range_value_t<InRange>, Complex>;
    requires std::ranges::output_range<OutRange,
                                       std::ranges::range_value_t<OutRange>>;
    requires std::same_as<RemoveComplex<std::ranges::range_value_t<OutRange>>,
                          Real>;
  }
  void InverseTransformation(Int lMax, Int n, InRange&& in,
                             OutRange& out) const {
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

    // Make the FFT plan.
    auto [inSize, outSize] = FFTWpp::DataSize<Complex, Scalar>(nPhi);
    auto inWork = FFTWpp::vector<Complex>(inSize);
    auto outWork = FFTWpp::vector<Scalar>(outSize);
    auto inView = FFTWpp::Ranges::View(inWork);
    auto outView = FFTWpp::Ranges::View(outWork);
    auto planFunction = [this](auto in, auto out) {
      if constexpr (ComplexFloatingPoint<Scalar>) {
        return FFTWpp::Ranges::Plan(in, out, _flag, FFTWpp::Backward);
      } else {
        return FFTWpp::Ranges::Plan(in, out, _flag);
      }
    };
    auto plan = planFunction(inView, outView);

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
          if constexpr (std::same_as<_MRange, All>) {
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
      auto outView = std::ranges::subrange(outStart, outFinish);
      plan.Execute(inView, outView);
    }
  }

 private:
  Int _lMax;
  Int _nMax;
  FFTWpp::Flag _flag;

  std::shared_ptr<QuadType> _quadPointer;
  std::shared_ptr<WignerType> _wignerPointer;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
