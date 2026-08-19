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
#include <stdexcept>
#include <string>
#include <vector>

#include "Concepts.h"
#include "GridBase.h"
#include "Indexing.h"
#include "Utility.h"
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

  // Constructors.
  GaussLegendreGrid() = default;

  GaussLegendreGrid(Int lMax, Int nMax, FFTWpp::Flag flag = FFTWpp::Measure)
      : _lMax{lMax}, _nMax{nMax}, _flag{flag} {
    // Check the inputs.
    assert(MaxDegree() >= 0);
    assert(MaxUpperIndex() <= MaxDegree());
    assert(std::abs(this->MinUpperIndex()) <= MaxDegree());
    assert(_flag != FFTWpp::WisdomOnly);

    // An MRange = NonNegative grid stores only m >= 0, so it cannot serve a
    // complex-valued transform at all, and its real-valued transforms exist
    // only at upper index zero. Such a grid with nMax != 0 could serve no
    // transform whatever, so it is a configuration error rather than a
    // wasteful but usable choice. See core-plan.md step A and [C1]: this is
    // the "real scalar grid" reading of MRange.
    if constexpr (std::same_as<_MRange, NonNegative>) {
      if (_nMax != 0) {
        throw std::invalid_argument(
            "A grid storing only non-negative orders serves real-valued "
            "transforms at upper index zero, so its maximum upper index must "
            "be zero");
      }
    }

    // Get the quadrature points.
    _quad = GaussQuad::LegendrePolynomial<Real>{}.GaussQuadrature(_lMax + 1);
    _quad.Transform([](auto x) { return std::acos(-x); },
                    [](auto x) -> Real { return 1; });

    //  Get the Winger values.
    _wigner = Wigner<Real, _MRange, _NRange, Multiple, ColumnMajor>(
        _lMax, _lMax, _nMax, _quad.Points());

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

  // A grid for working with fields of maximum degree lBand, with quadrature
  // headroom for degree oversampling * lBand.
  //
  // The distinction this expresses is between the *band* of the fields being
  // worked with and the *resolution* of the grid, which the library could not
  // previously say: the transforms already take lMax per call, so an
  // oversampled grid with truncated transforms works, but there was no way to
  // ask for one. Headroom is wanted whenever a quantity of higher degree than
  // the fields themselves is formed on the grid -- a product of two band-L
  // fields has band 2L, and integrating |f|^2 for band-L f integrates a
  // degree-2L quantity. The 3/2 rule is oversampling = 1.5; oversampling = 2
  // is exact for a single product.
  static auto ForBand(Int lBand, Int nMax, Real oversampling = 1,
                      FFTWpp::Flag flag = FFTWpp::Measure) {
    if (lBand < 0) {
      throw std::invalid_argument("Band must be non-negative");
    }
    if (!(oversampling >= 1)) {
      throw std::invalid_argument("Oversampling factor must be at least one");
    }
    const auto lGrid = static_cast<Int>(
        std::ceil(oversampling * static_cast<Real>(lBand)));
    return GaussLegendreGrid(lGrid, nMax, flag);
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

  auto CoLatitudes() const { return std::ranges::views::all(_quad.Points()); }
  auto CoLatitudeWeights() const {
    return std::ranges::views::all(_quad.Weights());
  }

  auto Longitudes() const {
    const auto nPhi = NPhi();
    const auto dPhi = 2 * std::numbers::pi_v<Real> / static_cast<Real>(nPhi);
    return std::ranges::views::iota(Int{0}, nPhi) |
           std::ranges::views::transform([dPhi](auto i) { return i * dPhi; });
  }
  auto LongitudeWeights() const {
    const auto nPhi = NPhi();
    const auto dPhi = 2 * std::numbers::pi_v<Real> / static_cast<Real>(nPhi);
    return std::ranges::views::repeat(dPhi, nPhi);
  }

  //-----------------------------------------------------//
  //          Forward transformation for ranges          //
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

    ValidateTransformRequest<Scalar>(lMax, n);
    CheckSize(std::ranges::size(in), this->FieldSize(), "field");
    CheckSize(std::ranges::size(out), CoefficientSizeFor<Scalar>(lMax, n),
              "coefficient");

    // The colatitude loop below is the quadrature sum, so it accumulates into
    // out. Nothing else initialises it, and every caller happened to arrive
    // with a zeroed buffer, which made "out arrives zeroed" an unstated and
    // unchecked precondition: transforming twice into the same buffer doubled
    // the answer (core-plan.md F1). Zero it here and the routine owns its
    // output. Accumulation, if ever wanted, is a separate named entry point.
    std::ranges::fill(out, Complex{});

    // A one-point grid needs no FFT.
    if (_lMax == 0) {
      out[0] =
          in[0] * static_cast<Real>(2) / std::numbers::inv_sqrtpi_v<Real>;
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
      // FFT the current data slice, through the plan's own buffer.
      PackRow(std::next(in.begin(), iTheta * nPhi), nPhi, inWork);
      plan.Execute();

      // Get the Wigner values and quadrature weight.
      auto d = _wigner[n, iTheta];
      auto w = _quad.W(iTheta) * scaleFactor;

      // Loop over the spherical harmonic coefficients
      auto outIter = out.begin();
      auto wigIter = d.begin();
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });

      for (auto l : degrees) {
        auto dl = d[l];

        if constexpr (ComplexFloatingPoint<Scalar>) {
          auto workIter = std::prev(outWork.end(), dl.MaxOrder());
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
            std::advance(wigIter, dl.MaxOrder());
          }
          for (auto m : dl.NonNegativeOrders()) {
            *outIter++ += *wigIter++ * *workIter++ * w;
          }
        }
      }

    }
  }

  //------------------------------------------------//
  //       Inverse  transformation for ranges       //
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

    ValidateTransformRequest<Scalar>(lMax, n);
    CheckSize(std::ranges::size(in), CoefficientSizeFor<Scalar>(lMax, n),
              "coefficient");
    CheckSize(std::ranges::size(out), this->FieldSize(), "field");

    // A one-point grid needs no FFT.
    if (_lMax == 0) {
      if constexpr (RealFloatingPoint<Scalar>) {
        out[0] = std::real(in[0]) * std::numbers::inv_sqrtpi_v<Real> /
                 static_cast<Real>(2);
      } else {
        out[0] = in[0] * std::numbers::inv_sqrtpi_v<Real> /
                 static_cast<Real>(2);
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
      auto d = _wigner[n, iTheta];

      // Loop over the coefficients.
      auto inIter = in.begin();
      auto wigIter = d.begin();
      auto degrees = d.Degrees() | std::ranges::views::filter(
                                       [lMax](auto l) { return l <= lMax; });
      for (auto l : degrees) {
        auto dl = d[l];
        if constexpr (ComplexFloatingPoint<Scalar>) {
          auto workIter = std::prev(inWork.end(), dl.MaxOrder());
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
            std::advance(wigIter, dl.MaxOrder());
          }
          for (auto m : dl.NonNegativeOrders()) {
            *workIter++ += *inIter++ * *wigIter++;
          }
        }
      }

      // Perform FFT to recover field at the colatitude, through the plan's
      // own buffer, then hand the row to the caller.
      plan.Execute();
      UnpackRow(outWork, nPhi, std::next(out.begin(), iTheta * nPhi));
    }
  }

 private:
  // Move one colatitude row between the caller's field and the plan's own
  // buffer.
  //
  // These are the only points at which caller storage is touched during a
  // transform. The FFT plans are created on the grid's own fftw_malloc'd
  // buffers and are executed on those buffers alone; FFTW's new-array execute
  // is valid only for storage with the same alignment characteristics as the
  // planning buffers, and neither FFTW nor FFTWpp checks (core-plan.md F3).
  // The forward transform used to take that path on any writable input and
  // the inverse took it unconditionally on the caller's output, which made
  // alignment an obligation propagating outward into every stride the field
  // and tensor layers might choose. Copying instead costs one pass over a row
  // against an FFT of the same row, and it is what makes the field plan's
  // promise -- that a slice target needs no alignment beyond Scalar's --
  // true rather than aspirational.
  //
  // They are also the seam step F widens: the (stride, dist) of a batch
  // ([C9]) enters here and nowhere else.
  template <typename Iterator, typename Buffer>
  static void PackRow(Iterator first, Int nPhi, Buffer& work) {
    std::copy_n(first, nPhi, work.begin());
  }

  template <typename Buffer, typename Iterator>
  static void UnpackRow(const Buffer& work, Int nPhi, Iterator first) {
    std::copy_n(work.begin(), nPhi, first);
  }

  // The coefficient count for a field of the given scalar type: reduced
  // m >= 0 storage for a real field, all orders for a complex one.
  // Named distinctly from GridBase::CoefficientSize, which it would otherwise
  // hide.
  template <RealOrComplexFloatingPoint Scalar>
  auto CoefficientSizeFor(Int lMax, Int n) const {
    if constexpr (RealFloatingPoint<Scalar>) {
      return GSHIndices<NonNegative>(lMax, lMax, n).Size();
    } else {
      return GSHIndices<All>(lMax, lMax, n).Size();
    }
  }

  // Size mismatches were assert-only, so under NDEBUG a short output range was
  // a silent heap overflow (core-plan.md F5). Checked in all build modes.
  static void CheckSize(std::size_t given, std::integral auto expected,
                        const char* what) {
    if (given != static_cast<std::size_t>(expected)) {
      throw std::invalid_argument(
          std::string("Transform ") + what + " range has size " +
          std::to_string(given) + ", but this request needs " +
          std::to_string(expected));
    }
  }

  // The longitude quadrature is the trapezoid rule on nPhi equally spaced
  // points, which is exact for exp(i (m - m') phi) only when |m - m'| < nPhi.
  // Resolving orders |m| <= lMax therefore needs nPhi >= 2 * lMax + 1, not
  // 2 * lMax: at 2 * lMax the orders m = +lMax and m = -lMax are the same
  // discrete mode and cannot be separated, which is why the transform used to
  // zero the (lMax, lMax) coefficient rather than compute it (core-plan.md
  // F2). The smallest fast FFT length at or above the bound is used, so that
  // the fix does not land on a length with a large prime factor.
  auto NPhi() const { return FastFFTSize(2 * _lMax + 1); }

  template <RealOrComplexFloatingPoint Scalar>
  void ValidateTransformRequest(Int lMax, Int n) const {
    if (lMax < 0 || lMax > MaxDegree()) {
      throw std::invalid_argument(
          "Transform degree must be between zero and the grid maximum degree");
    }
    if (std::abs(n) > lMax ||
        !std::ranges::contains(this->UpperIndices(), n)) {
      throw std::invalid_argument(
          "Transform upper index is not supported at the requested degree");
    }

    // A spin-weighted field of nonzero upper index cannot be real-valued:
    // real-valuedness is not preserved by the frame rotation
    // e_{+-} -> e^{-+ i psi} e_{+-}, so it is not a property any component of
    // any tensor can have at N != 0 (theory note section 7, item 5). The
    // reduced m >= 0 coefficient storage that a real transform uses assumes
    // the self-relation f^N_{l,-m} = (-1)^{m-N} conj(f^N_{lm}), which holds
    // only when f is its own conjugate, i.e. only at N = 0. n is a runtime
    // argument, so this is a throw rather than a static_assert
    // (core-plan.md step A, [C2]).
    if constexpr (RealFloatingPoint<Scalar>) {
      if (n != 0) {
        throw std::invalid_argument(
            "Real-valued fields exist only at upper index zero");
      }
    }
  }

  Int _lMax;
  Int _nMax;
  FFTWpp::Flag _flag;

  GaussQuad::Quadrature1D<Real> _quad;
  Wigner<Real, _MRange, _NRange, Multiple, ColumnMajor> _wigner;

  // std::shared_ptr<QuadType> _quadPointer;
  // std::shared_ptr<WignerType> _wignerPointer;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
