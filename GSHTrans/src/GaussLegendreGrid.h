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
 public:
  using real_type = Real;
  using complex_type = std::complex<Real>;

 private:
  using Int = std::ptrdiff_t;
  using Complex = complex_type;

  using WignerType = Wigner<Real, MRange, NRange, Ortho>;
  using QuadType = GaussQuad::Quadrature1D<Real>;

 public:
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

    // Set data layouts for possible FFTs.
    _inLayout = FFTWpp::DataLayout(1, std::vector{NumberOfLongitudes()}, 1,
                                   std::vector{NumberOfLongitudes()}, 1, 1);
    _outLayoutNonNegative = FFTWpp::DataLayout(1, std::vector{_lMax + 1}, 1,
                                               std::vector{_lMax + 1}, 1, 1);
    if constexpr (std::same_as<MRange, All>) {
      _outLayoutAll = FFTWpp::DataLayout(1, std::vector{2 * _lMax}, 1,
                                         std::vector{2 * _lMax}, 1, 1);
    }

    // Generate FFT wisdom.
    FFTWpp::GenerateWisdom<Real, Complex, true>(_inLayout,
                                                _outLayoutNonNegative, flag);
    if constexpr (std::same_as<MRange, All>) {
      FFTWpp::GenerateWisdom<Complex, Complex, true>(_inLayout, _outLayoutAll,
                                                     flag);
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

  // Forward transformation.
  template <RealOrComplexFloatingPointIterator InIterator,
            ComplexFloatingPointIterator OutIterator>
  void ForwardTransformation(Int lMax, Int n, InIterator in,
                             OutIterator out) const {
    // Determine whether this is a transformation of real data.
    constexpr auto realTransform = RealFloatingPointIterator<InIterator>;
    constexpr auto complexTransform = ComplexFloatingPointIterator<InIterator>;

    // Check appropriate orders are present.
    if constexpr (complexTransform) {
      constexpr auto all = std::same_as<MRange, All>;
      assert(all);
    }

    // Check upper index is possible.
    assert(
        std::ranges::any_of(UpperIndices(), [n](auto np) { return n == np; }));

    // Pre compute some constants.
    const auto nPhi = NumberOfLongitudes();
    const auto scaleFactor = static_cast<Real>(2) * std::numbers::pi_v<Real> /
                             static_cast<Real>(nPhi);

    // Allocate work array for FFTs and form its FFTWpp view.
    auto [work, outView] = FFTWork<realTransform>();

    // Loop over the colatitudes.
    for (auto iTheta : CoLatitudeIndices()) {
      // Form FFT of current data slice.
      auto offset = iTheta * nPhi;
      auto inStart = std::next(in, offset);
      auto inFinish = std::next(inStart, nPhi);
      auto inView = FFTWpp::DataView(inStart, inFinish, _inLayout);
      auto plan =
          FFTWpp::Plan(inView, outView, FFTWpp::WisdomOnly, FFTWpp::Forward);
      plan.Execute();

      // Get the Wigner values and quadrature weight.
      auto d = _wigner->operator()(n)(iTheta);
      auto w = _quad.W(iTheta) * scaleFactor;

      // Loop over the spherical harmonic coefficients
      auto outIter = out;
      auto wigIter = d.cbegin();
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });
      for (auto l : degrees) {
        auto dl = d(l);

        if constexpr (complexTransform) {
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

    // Zero the (_lMax,_lMax) coefficient if needed.
    if constexpr (complexTransform) {
      if (lMax = _lMax) {
        auto view = GSHView<Complex, MRange>(_lMax, _lMax, n, out);
        view(_lMax)(_lMax) = 0;
      }
    }
  }

  template <RealOrComplexFloatingPointIterator InIterator,
            ComplexFloatingPointIterator OutIterator>
  void ForwardTransformation(Int n, InIterator in, OutIterator out) const {
    ForwardTransformation(_lMax, n, in, out);
  }

  // Inverse  transformation.
  template <ComplexFloatingPointIterator InIterator,
            RealOrComplexFloatingPointIterator OutIterator>
  void InverseTransformation(Int lMax, Int n, InIterator in,
                             OutIterator out) const {
    // Determine whether this is a transformation of real data.
    constexpr auto realTransform = RealFloatingPointIterator<OutIterator>;
    constexpr auto complexTransform = ComplexFloatingPointIterator<OutIterator>;

    // Check appropriate orders are present.
    if constexpr (complexTransform) {
      constexpr auto all = std::same_as<MRange, All>;
      assert(all);
    }

    // Check upper index is possible.
    assert(
        std::ranges::any_of(UpperIndices(), [n](auto np) { return n == np; }));

    // Precompute constants
    const auto nPhi = NumberOfLongitudes();

    // Allocate work array for FFTs and form its FFTWpp view.
    auto [work, inView] = FFTWork<realTransform>();

    // Loop over the colatitudes.
    for (auto iTheta : CoLatitudeIndices()) {
      std::ranges::for_each(work, [](auto x) { return 0; });

      // Get the Wigner values.
      auto d = (*_wigner)(n)(iTheta);

      // Loop over the coefficients.
      auto inIter = in;
      auto wigIter = d.begin();
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });
      for (auto l : degrees) {
        auto dl = d(l);
        if constexpr (complexTransform) {
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

      for (auto val : work) std::cout << val << "    ";
      std::cout << std::endl;

      // Perform FFT to recover field at the colatitude.
      auto offset = iTheta * nPhi;
      auto outStart = std::next(out, offset);
      auto outFinish = std::next(outStart, nPhi);
      auto outView = FFTWpp::DataView(outStart, outFinish, _inLayout);
      auto plan =
          FFTWpp::Plan(inView, outView, FFTWpp::WisdomOnly, FFTWpp::Backward);
      plan.Execute();
    }
  }

  template <ComplexFloatingPointIterator InIterator,
            RealOrComplexFloatingPointIterator OutIterator>
  void InverseTransformation(Int n, InIterator in, OutIterator out) const {
    InverseTransformation(_lMax, n, in, out);
  }

 private:
  Int _lMax;
  Int _nMax;

  QuadType _quad;

  std::shared_ptr<WignerType> _wigner;

  FFTWpp::DataLayout _inLayout;
  FFTWpp::DataLayout _outLayoutAll;
  FFTWpp::DataLayout _outLayoutNonNegative;

  template <bool RealTransform>
  auto FFTWork() const {
    using Vector = FFTWpp::vector<Complex>;
    using Iterator = Vector::iterator;
    auto work = Vector();
    auto view = FFTWpp::DataView<Iterator>();
    if constexpr (RealTransform) {
      work.resize(_lMax + 1);
      view = FFTWpp::DataView(work.begin(), work.end(), _outLayoutNonNegative);
    } else {
      work.resize(2 * _lMax);
      view = FFTWpp::DataView(work.begin(), work.end(), _outLayoutAll);
    }
    return std::pair(work, view);
  }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
