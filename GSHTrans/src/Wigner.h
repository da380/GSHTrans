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
#include "Views.h"

namespace GSHTrans {

namespace WignerDetails {

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

}  // namespace WignerDetails

template <RealFloatingPoint _Real, Normalisation _Norm = Ortho,
          OrderIndexRange _MRange = All, IndexRange _NRange = Single,
          AngleIndexRange _AngleRange = Single,
          WignerStorage _Storage = ColumnMajor>
class Wigner {
 public:
  using Int = std::ptrdiff_t;
  using Real = _Real;
  using Norm = _Norm;
  using MRange = _MRange;
  using NRange = _NRange;
  using AngleRange = _AngleRange;
  using Storage = _Storage;

  // Default constructor.
  Wigner() = default;

  // General constructor.
  template <std::ranges::range Range>
  requires RealFloatingPoint<std::ranges::range_value_t<Range>>
  Wigner(Int lMax, Int mMax, Int nMax, Range &&theta)
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
    ComputeAll(theta);
  }

  // Constructor using a single angle.
  Wigner(Int lMax, Int mMax, Int nMax, Real theta)
  requires std::same_as<AngleRange, Single>
      : Wigner(lMax, mMax, nMax, std::vector{theta}) {}

  // Recompute for new angles.
  template <std::ranges::range Range>
  requires RealFloatingPoint<std::ranges::range_value_t<Range>>
  void ReCompute(Range &&theta) {
    assert(theta.size() == NumberOfAngles());
    ComputeAll(theta);
  }

  void ReCompute(Real theta)
  requires std::same_as<AngleRange, Single>
  {
    ComputeAll(std::vector{theta});
  }

  // Return degree information
  auto MinDegree(Int n) const {
    assert(std::abs(n) <= _nMax);
    return std::abs(n);
  }
  auto MaxDegree() const { return _lMax; }

  auto Degrees(Int n) const {
    return std::ranges::views::iota(MinDegree(n), MaxDegree() + 1);
  }

  // Return order information.
  auto MinOrder(Int l) const {
    assert(l >= 0 && l <= MaxDegree());
    if constexpr (std::same_as<MRange, All>) {
      return -std::min(l, _mMax);
    } else {
      return 0;
    }
  }
  auto MaxOrder() const { return _mMax; }
  auto MaxOrder(Int l) const { return std::max(l, _mMax); }

  auto Orders(Int l) const {
    return std::ranges::views::iota(MinOrder(l), MaxOrder(l) + 1);
  }

  // Return upper index information.
  auto MinUpperIndex() const {
    if constexpr (std::same_as<_NRange, All>) {
      return -_nMax;
    } else if constexpr (std::same_as<_NRange, NonNegative>) {
      return Int{0};
    } else {
      return _nMax;
    }
  }

  auto MaxUpperIndex() const { return _nMax; }

  auto UpperIndices() const {
    return std::ranges::views::iota(MinUpperIndex(), MaxUpperIndex() + 1);
  }

  auto NumberOfUpperIndices() const {
    return MaxUpperIndex() - MinUpperIndex() + 1;
  }

  // Return angle information.
  auto NumberOfAngles() const { return _nTheta; }

  auto AngleIndices() const {
    return std::ranges::views::iota(Int{0}, _nTheta);
  }

  // Return (n,iTheta) values in order.
  auto Indices() const {
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

  // Return view to data for (n,iTheta).
  auto operator[](Int n, Int iTheta) const {
    return ConstGSHView<Real, MRange>(MaxDegree(), MaxOrder(), n,
                                      &_data[Offset(n, iTheta)]);
  }

  // Return view to data for n when AngleRange = Single.
  auto operator[](Int n) const
  requires std::same_as<AngleRange, Single> && (!std::same_as<NRange, Single>)
  {
    return operator[](n, 0);
  }

  // Return view to data for iTheta when NRange = Single.
  auto operator[](Int iTheta) const
  requires std::same_as<NRange, Single> && (!std::same_as<AngleRange, Single>)
  {
    return operator[](0, iTheta);
  }

  // Return subview to data for l when NRange = Single and AngleRange = Single
  auto operator[](Int l) const
  requires std::same_as<NRange, Single> && std::same_as<AngleRange, Single>
  {
    return operator[](0, 0)[l];
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

  // Offset for values for (n,iTheta).
  auto Offset(Int n, Int iTheta) const {
    if constexpr (std::same_as<Storage, ColumnMajor>) {
      return _offset[NumberOfAngles() * (n - MinUpperIndex()) + iTheta];
    } else {
      return _offset[NumberOfUpperIndices() * iTheta + (n - MinUpperIndex())];
    }
  }

  // Pre-compute some numerical values that are repeatedly used.
  auto PreCompute() const {
    using Vector = std::vector<Real>;
    auto size = MaxDegree() + std::max(MaxOrder(), MaxUpperIndex()) + 1;
    Vector sqrtInt, sqrtIntInv;
    sqrtInt.reserve(size);
    sqrtIntInv.reserve(size);
    std::generate_n(std::back_inserter(sqrtInt), size, [m = Int{0}]() mutable {
      return std::sqrt(static_cast<Real>(m++));
    });
    std::transform(sqrtInt.begin(), sqrtInt.end(),
                   std::back_inserter(sqrtIntInv),
                   [](auto x) { return x > 0 ? 1 / x : 0; });
    return std::tuple(std::make_shared<Vector>(sqrtInt),
                      std::make_shared<Vector>(sqrtIntInv));
  }

  template <std::ranges::range Range>
  requires RealFloatingPoint<std::ranges::range_value_t<Range>>
  void ComputeAll(Range &&thetaRange) {
    auto preComputed = PreCompute();

#pragma omp parallel for
    for (auto [n, iTheta] : Indices()) {
      Compute(n, iTheta, thetaRange[iTheta], preComputed);
    }
  }

  // Compute values for given (m,iTheta).
  constexpr void Compute(Int n, Int iTheta, Real theta,
                         const auto preComputed) {
    // Set some values.
    auto &sqrtInt = *std::get<0>(preComputed);
    auto &sqrtIntInv = *std::get<1>(preComputed);
    auto arg = WignerDetails::Arguments(theta);
    const auto cos = std::cos(theta);
    const auto nAbs = std::abs(n);

    // Set view to the data for (n,iTheta).
    auto d = GSHView<Real, MRange>(_lMax, _mMax, n, &_data[Offset(n, iTheta)]);

    // Set the values for l == |n|
    {
      const auto l = nAbs;
      auto m = d[l].MinOrder();
      auto iter = d[l].begin();
      auto finish = d[l].end();
      if (n >= 0) {
        while (iter != finish) {
          *iter++ = WignerDetails::WignerMaxUpperIndex(l, m++, arg);
        }
      } else {
        while (iter != finish) {
          *iter++ = WignerDetails::WignerMinUpperIndex(l, m++, arg);
        }
      }
    }

    // Set the values for l == n+1 if needed.
    if (nAbs < _lMax) {
      const auto l = nAbs + 1;
      const auto mMin = d[l].MinOrder();
      const auto mMax = d[l].MaxOrder();
      auto m = mMin;

      // Set iterators
      auto iterMinusOne = d[l - 1].begin();
      auto finishMinusOne = d[l - 1].end();
      auto iter = d[l].begin();

      // Add in value at m == -l if needed.
      if constexpr (std::same_as<_MRange, All>) {
        if (l <= mMax) {
          *iter++ = WignerDetails::WignerMinOrder(l, n, arg);
          m++;
        }
      }

      // Add in interior orders using one-term recursion.
      {
        const auto alpha = (2 * l - 1) * l * cos * sqrtIntInv[l + nAbs];
        const auto beta = (n < 0 ? -1 : 1) * (2 * l - 1) * sqrtIntInv[l + nAbs];

        while (iterMinusOne != finishMinusOne) {
          const auto f1 =
              (alpha - beta * m) * sqrtIntInv[l - m] * sqrtIntInv[l + m];
          *iter++ = f1 * *iterMinusOne++;
          m++;
        }
      }

      // Add in value at m == l if needed
      if (l <= mMax) {
        *iter++ = WignerDetails::WignerMaxOrder(l, n, arg);
      }
    }

    // Do the remaining degrees
    for (auto l = nAbs + 2; l <= _lMax; l++) {
      const auto mMin = d[l].MinOrder();
      const auto mMax = d[l].MaxOrder();
      auto m = mMin;

      // Set iterators.
      auto iterMinusTwo = d[l - 2].begin();
      auto finishMinusTwo = d[l - 2].end();
      auto iterMinusOne = d[l - 1].begin();
      auto iter = d[l].begin();

      // Add in lower boundary terms if still growing.
      if constexpr (std::same_as<_MRange, All>) {
        if (l <= mMax) {
          {
            // Add in the m == -l term.
            *iter++ = WignerDetails::WignerMinOrder(l, n, arg);
            m++;
          }
          {
            // Now do the m == -l+1 term using one-point recursion.
            const auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) *
                            sqrtIntInv[l - n] * sqrtIntInv[l + n] *
                            sqrtIntInv[l - m] * sqrtIntInv[l + m] /
                            static_cast<Real>(l - 1);
            *iter++ = f1 * (*iterMinusOne++);
            m++;
          }
        }

        // Add in the lower boundary term at the critical degree.
        if (l == mMax + 1) {
          const auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) *
                          sqrtIntInv[l - n] * sqrtIntInv[l + n] *
                          sqrtIntInv[l - m] * sqrtIntInv[l + m] /
                          static_cast<Real>(l - 1);
          *iter++ = f1 * (*iterMinusOne++);
          m++;
        }
      }

      // Apply two-term recusion for the interior orders.
      {
        const auto alpha =
            (2 * l - 1) * l * cos * sqrtIntInv[l - n] * sqrtIntInv[l + n];
        const auto beta = (2 * l - 1) * n * sqrtIntInv[l - n] *
                          sqrtIntInv[l + n] / static_cast<Real>(l - 1);
        const auto gamma = l * sqrtInt[l - 1 - n] * sqrtInt[l - 1 + n] *
                           sqrtIntInv[l - n] * sqrtIntInv[l + n] /
                           static_cast<Real>(l - 1);

        while (iterMinusTwo != finishMinusTwo) {
          const auto denom = sqrtIntInv[l - m] * sqrtIntInv[l + m];
          const auto f1 = (alpha - beta * m) * denom;
          const auto f2 =
              gamma * sqrtInt[l - 1 - m] * sqrtInt[l - 1 + m] * denom;
          *iter++ = f1 * *iterMinusOne++ - f2 * *iterMinusTwo++;
          m++;
        }
      }

      // Add in the upper boundary terms if still growing.
      if (l <= mMax) {
        // Add in m == l - 1 term using one-point recursion.
        {
          const auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) *
                          sqrtIntInv[l - n] * sqrtIntInv[l + n] *
                          sqrtIntInv[l - m] * sqrtIntInv[l + m] /
                          static_cast<Real>(l - 1);
          *iter++ = f1 * (*iterMinusOne++);
          m++;
        }
        // Now do m == l.
        *iter++ = WignerDetails::WignerMaxOrder(l, n, arg);
      }

      // Add in the upper boundary term at the critical degree.
      if (l == mMax + 1) {
        // Update the iterators.

        const auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) *
                        sqrtIntInv[l - n] * sqrtIntInv[l + n] *
                        sqrtIntInv[l - m] * sqrtIntInv[l + m] /
                        static_cast<Real>(l - 1);
        *iter++ = f1 * (*iterMinusOne++);
      }
    }

    // Normalise the values if needed.
    if constexpr (std::same_as<Norm, Ortho>) {
      const auto factor =
          std::numbers::inv_sqrtpi_v<Real> / static_cast<Real>(2);
      for (auto l : d.Degrees()) {
        auto start = d[l].begin();
        auto finish = d[l].end();
        std::transform(start, finish, start, [l, factor](auto p) {
          return factor * std::sqrt(static_cast<Real>(2 * l + 1)) * p;
        });
      }
    }
  }
};
}  // namespace GSHTrans

#endif