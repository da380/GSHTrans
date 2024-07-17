#ifndef GSH_TRANS_WIGNER_NEW_GUARD_H
#define GSH_TRANS_WIGNER_NEW_GUARD_H

#include <omp.h>

#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <numbers>
#include <ranges>
#include <vector>

#include "Concepts.h"
#include "Indexing.h"
#include "Utility.h"

namespace GSHTrans {

namespace WignerNewDetails {

template <RealFloatingPoint Real>
class Arguments {
 public:
  Arguments() = default;

  constexpr Arguments(Real theta) {
    constexpr auto half = static_cast<Real>(1) / static_cast<Real>(2);
    _logSinHalf = std::sin(half * theta);
    _logCosHalf = std::cos(half * theta);
    _atLeft = _logSinHalf < std::numeric_limits<Real>::min();
    _atRight = _logCosHalf < std::numeric_limits<Real>::min();
    _logSinHalf = _atLeft ? static_cast<Real>(0) : std::log(_logSinHalf);
    _logCosHalf = _atRight ? static_cast<Real>(0) : std::log(_logCosHalf);
  }

  constexpr auto AtLeft() const { return _atLeft; }
  constexpr auto AtRight() const { return _atRight; }

  constexpr auto LogSinHalf() const { return _logSinHalf; }
  constexpr auto LogCosHalf() const { return _logCosHalf; }

 private:
  Real _logSinHalf;
  Real _logCosHalf;
  bool _atLeft;
  bool _atRight;
};

template <std::integral Int, RealFloatingPoint Real>
constexpr auto WignerMinOrder(Int l, Int n, const Arguments<Real> &arg) {
  // Check the inputs.
  assert(l >= 0);
  assert(std::abs(n) <= l);

  // Deal with l == 0 case
  if (l == 0) return static_cast<Real>(1);

  // Deal with special case at the left boundary.
  if (arg.AtLeft()) {
    return n == -l ? static_cast<Real>(1) : static_cast<Real>(0);
  }

  // Deal with special case at the right boundary.
  if (arg.AtRight()) {
    return n == l ? static_cast<Real>(1) : static_cast<Real>(0);
  }

  // Deal with the general case.
  constexpr auto half = static_cast<Real>(1) / static_cast<Real>(2);
  auto Fl = static_cast<Real>(l);
  auto Fn = static_cast<Real>(n);
  using std::exp;
  using std::lgamma;
  return exp(
      half * (lgamma(2 * Fl + 1) - lgamma(Fl - Fn + 1) - lgamma(Fl + Fn + 1)) +
      (Fl + Fn) * arg.LogSinHalf() + (Fl - Fn) * arg.LogCosHalf());
}

template <std::integral Int, RealFloatingPoint Real>
constexpr auto WignerMaxOrder(Int l, Int n, Arguments<Real> &arg) {
  return MinusOneToPower(n + l) * WignerMinOrder(l, -n, arg);
}

template <std::integral Int, RealFloatingPoint Real>
constexpr auto WignerMinUpperIndex(Int l, Int m, Arguments<Real> &arg) {
  return WignerMaxOrder(l, -m, arg);
}

template <std::integral Int, RealFloatingPoint Real>
constexpr auto WignerMaxUpperIndex(Int l, Int m, Arguments<Real> &arg) {
  return WignerMinOrder(l, -m, arg);
}

}  // namespace WignerNewDetails

template <RealFloatingPoint _Real, Normalisation _Norm = Ortho,
          OrderIndexRange _MRange = All, IndexRange _NRange = Single,
          AngleIndexRange _AngleRange = Single,
          WignerStorage _Storage = ColumnMajor>
class WignerNew {
 public:
  using Int = std::ptrdiff_t;
  using Real = _Real;
  using Norm = _Norm;
  using MRange = _MRange;
  using NRange = _NRange;
  using AngleRange = _AngleRange;
  using Storage = _Storage;

  // General Constructor.
  template <std::ranges::range Range>
  requires RealFloatingPoint<std::ranges::range_value_t<Range>>
  WignerNew(Int lMax, Int mMax, Int nMax, Range &&theta)
      : _lMax{lMax}, _mMax{mMax}, _nMax{nMax}, _nTheta(theta.size()) {
    // Compute the offsets and allocate memory.
    {
      _offset.reserve(NumberOfUpperIndices() * NumberOfAngles());
      auto size = std::size_t{0};
      for (auto [n, iTheta] : Indices()) {
        _offset.push_back(size);
        size += GSHIndices<MRange>(MaxDegree(), _mMax, n).size();
      }
      _data = std::vector<Real>(size);
    }
  }

  // Return degree information
  constexpr auto MinDegree(Int n) const {
    assert(std::abs(n) <= _nMax);
    return std::abs(n);
  }
  constexpr auto MaxDegree() const { return _lMax; }

  constexpr auto Degrees(Int n) const {
    return std::ranges::views::iota(MinDegree(n), MaxDegree() + 1);
  }

  // Return order information.
  constexpr auto MinOrder(Int l) const {
    assert(l >= 0 && l <= MaxDegree());
    if constexpr (std::same_as<MRange, All>) {
      return -std::min(l, _mMax);
    } else {
      return 0;
    }
  }
  constexpr auto MaxOrder(Int l) const { return std::max(l, _mMax); }

  constexpr auto Orders(Int l) const {
    return std::ranges::views::iota(MinOrder(l), MaxOrder(l) + 1);
  }

  // Return upper index information.
  constexpr auto MinUpperIndex() const {
    if constexpr (std::same_as<_NRange, All>) {
      return -_nMax;
    } else if constexpr (std::same_as<_NRange, NonNegative>) {
      return Int{0};
    } else {
      return _nMax;
    }
  }

  constexpr auto MaxUpperIndex() const { return _nMax; }

  constexpr auto UpperIndices() const {
    return std::ranges::views::iota(MinUpperIndex(), MaxUpperIndex() + 1);
  }

  constexpr auto NumberOfUpperIndices() const {
    return MaxUpperIndex() - MinUpperIndex() + 1;
  }

  // Return angle information.
  constexpr auto NumberOfAngles() const { return _nTheta; }

  auto AngleIndices() const {
    return std::ranges::views::iota(Int{0}, _nTheta);
  }

  constexpr auto Indices() const {
    if constexpr (std::same_as<Storage, ColumnMajor>) {
      return std::ranges::views::cartesian_product(UpperIndices(),
                                                   AngleIndices());
    } else {
      return std::ranges::views::cartesian_product(AngleIndices(),
                                                   UpperIndices()) |
             std::ranges::views::transform(
                 [](auto p) { return std::pair(p.second, p.first); });
    }
  }

 private:
  Int _lMax;    // Maximum degree.
  Int _mMax;    // Maximum order.
  Int _nMax;    // Maximum upper index.
  Int _nTheta;  // Number of colatitudes.

  // Vector storing the values.
  std::vector<Real> _data;

  // Vector storing data offsets.
  std::vector<std::size_t> _offset;

  // Return index of the start of values for (n,iTheta).
  // storage.
  constexpr auto Offset(Int n, Int iTheta) const {
    if constexpr (std::same_as<Storage, ColumnMajor>) {
      return _offset[NumberOfAngles() * (n - MinUpperIndex()) + iTheta];
    } else {
      return _offset[NumberOfUpperIndices() * iTheta + (n - MinUpperIndex())];
    }
  }

  // Iterators to the start of Data for (n,iTheta).
  constexpr auto BeginForUpperIndexAndAngle(Int n, Int iTheta) const {
    return std::next(_data.begin(), Offset(n, iTheta));
  }

  // Pre-compute some numerical values.
  auto PreCompute() const {
    auto size = MaxDegree() + std::max(MaxOrder(), MaxUpperIndex()) + 1;
    std::vector<Real> sqrtInt, sqrtIntInv;
    sqrtInt.reserve(size);
    sqrtIntInv.reserve(size);
    std::generate_n(std::back_inserter(sqrtInt), size, [m = Int{0}]() mutable {
      return std::sqrt(static_cast<Real>(m++));
    });
    std::transform(sqrtInt.begin(), sqrtInt.end(),
                   std::back_inserter(sqrtIntInv),
                   [](auto x) { return x > 0 ? 1 / x : 0; });
    return std::tuple(std::make_shared(sqrtInt), std::make_shared(sqrtIntInv));
  }

  // Compute values for given (m,iTheta).
  constexpr void Compute(Real theta, Int n, Int iTheta, const auto preCompute) {
    // Set some values.
    auto &sqrtInt = *std::get<0>(preCompute);
    auto &sqrtIntInv = *std::get<1>(preCompute);
    auto arg = WignerDetails::Arguments(theta);
    const auto cos = std::cos(theta);
    const auto nAbs = std::abs(n);

    index = GSHIndex<MRange>(_lMax, _mMax, n);

    /*

      // Compute values for l == |n|
      {
        const auto l = nAbs;
        auto iter = DataIterator(n, iTheta);
        if (n >= 0) {
          for (auto m : Orders(l)) {
            *iter++ = WignerDetails::WignerMaxUpperIndex(l, m, arg);
          }
        } else {
          for (auto m : Orders(l)) {
            *iter++ = WignerDetails::WignerMinUpperIndex(l, m, arg);
          }
        }
      }
  */
  }
};
}  // namespace GSHTrans

#endif