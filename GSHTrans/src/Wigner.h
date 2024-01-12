#ifndef GSH_TRANS_WIGNER_GUARD_H
#define GSH_TRANS_WIGNER_GUARD_H

#include <omp.h>

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include "Concepts.h"

namespace GSHTrans {

template <RealFloatingPoint Real, OrderIndexRange MRange, IndexRange NRange,
          Normalisation Norm = Ortho>
class Wigner {
  using Vector = std::vector<Real>;
  using Int = Vector::difference_type;
  using Mat = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
  using RowVec = Eigen::Matrix<Real, 1, Eigen::Dynamic>;

 public:
  // Set member types.
  using value_type = Real;
  using iterator = Vector::iterator;
  using const_iterator = Vector::const_iterator;
  using difference_type = Vector::difference_type;
  using size_type = Vector::size_type;

  Wigner() = default;

  // Constructors for a Upper index range.
  template <RealFloatingPointRange RealRange>
  Wigner(Int lMax, Int mMax, Int nMax, RealRange &&thetaRange)
      : _lMax{lMax}, _mMax{mMax}, _nMax{nMax}, _nTheta(thetaRange.size()) {
    AllocateStorage();
    auto preCompute = PreCompute();
    for (auto n : UpperIndices()) {
      auto iTheta = Int{0};
      for (auto theta : thetaRange) {
        auto d = operator()(n)(iTheta++);
        Compute(n, theta, d, preCompute);
      }
    }
  }

  Wigner(Int lMax, Int mMax, Int nMax, Real theta)
      : Wigner(lMax, mMax, nMax, std::vector{theta}) {}

  // Return basic information.
  auto MaxDegree() const { return _lMax; }
  auto MaxOrder() const { return _mMax; }

  auto Degrees() const { return std::ranges::views::iota(0, _lMax + 1); }
  auto NumberOfAngles() const { return _nTheta; }
  auto size() const { return _data.size(); }
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  auto MinUpperIndex() const {
    if constexpr (std::same_as<NRange, All>) {
      return -_nMax;
    } else if constexpr (std::same_as<NRange, NonNegative>) {
      return Int{0};
    } else {
      return _nMax;
    }
  }

  auto MaxUpperIndex() const { return _nMax; }

  auto UpperIndices() const {
    return std::ranges::views::iota(MinUpperIndex(), MaxUpperIndex() + 1);
  }

  // Return view to angle indices.
  auto AngleIndices() const {
    return std::ranges::views::iota(Int{0}, _nTheta);
  }

  // Return view to all angles for given upper index.
  auto operator()(Int n) {
    auto upperIndices = UpperIndices() | std::ranges::views::filter(
                                             [n](auto np) { return np < n; });
    auto offset = std::accumulate(
        upperIndices.begin(), upperIndices.end(), Int{0},
        [this](auto acc, auto np) {
          return acc +
                 GSHIndices<MRange>(_lMax, _mMax, np).size() * NumberOfAngles();
        });
    auto start = std::next(begin(), offset);
    return GSHViewAngleRange<Real, MRange>(_lMax, _mMax, n, _nTheta, start);
  }

  // Application operator in the case of a single upper index.
  auto operator()() requires std::same_as<NRange, Single> {
    return operator()(_nMax);
  }

  // Extract Wigner matrix at given angle and degree indices.
  auto Matrix(Int iTheta, Int l) requires requires() {
    std::same_as<MRange, All>;
    std::same_as<NRange, All>;
  }
  {
    auto mat = Mat(2 * l + 1, 2 * l + 1);
    auto upperIndices =
        UpperIndices() |
        std::ranges::views::filter([l](auto n) { return std::abs(n) <= l; });
    for (auto n : upperIndices) {
      auto d = operator()(n)(iTheta);
      auto i = n + l;
      for (auto m : d(l).Orders()) {
        auto j = m + l;
        mat(i, j) = d(l)(m);
      }
    }
    return mat;
  }

 private:
  Int _lMax;                // Maximum degree.
  Int _mMax;                // Maximum order.
  Int _nMax;                // Maximum upper index.
  Int _nTheta;              // Number of colatitudes.

  // Vector storing the values.
  Vector _data;

  // Compute the necessary storage capacity.
  void AllocateStorage() {
    auto upperIndices = UpperIndices();
    auto size = std::accumulate(
        upperIndices.begin(), upperIndices.end(), Int{0},
        [this](auto acc, auto n) {
          return acc +
                 GSHIndices<MRange>(_lMax, _mMax, n).size() * NumberOfAngles();
        });
    _data.resize(size);
  }

  // Pre-compute some numerical terms used repeatedly within
  // calculation of the Wigner values for each (theta,n) pair.
  auto PreCompute() const {
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

  constexpr auto MinusOneToPower(Int m) { return m % 2 ? -1 : 1; }

  auto WignerMinOrder(Int l, Int n, Real logSinHalf, Real logCosHalf,
                      bool atLeft, bool atRight) {
    // Check the inputs.
    assert(l >= 0);
    assert(std::abs(n) <= l);

    // Deal with l == 0 case
    if (l == 0) return static_cast<Real>(1);

    // Deal with special case at the left boundary.
    if (atLeft) {
      return n == -l ? static_cast<Real>(1) : static_cast<Real>(0);
    }

    // Deal with special case at the right boundary.
    if (atRight) {
      return n == l ? static_cast<Real>(1) : static_cast<Real>(0);
    }

    // Deal with the general case.
    auto Fl = static_cast<Real>(l);
    auto Fn = static_cast<Real>(n);
    using std::exp;
    using std::lgamma;
    return exp(
        static_cast<Real>(0.5) *
            (lgamma(2 * Fl + 1) - lgamma(Fl - Fn + 1) - lgamma(Fl + Fn + 1)) +
        (Fl + Fn) * logSinHalf + (Fl - Fn) * logCosHalf);
  }

  auto WignerMaxOrder(Int l, Int n, Real logSinHalf, Real logCosHalf,
                      bool atLeft, bool atRight) {
    return MinusOneToPower(n + l) *
           WignerMinOrder(l, -n, logSinHalf, logCosHalf, atLeft, atRight);
  }

  auto WignerMinUpperIndex(Int l, Int m, Real logSinHalf, Real logCosHalf,
                           bool atLeft, bool atRight) {
    return WignerMaxOrder(l, -m, logSinHalf, logCosHalf, atLeft, atRight);
  }

  auto WignerMaxUpperIndex(Int l, Int m, Real logSinHalf, Real logCosHalf,
                           bool atLeft, bool atRight) {
    return WignerMinOrder(l, -m, logSinHalf, logCosHalf, atLeft, atRight);
  }

  void Compute(Int n, Real theta, GSHView<Real, MRange> d,
               const auto preCompute) {
    // Get references to the pre-computed values.
    auto &sqrtInt = *std::get<0>(preCompute);
    auto &sqrtIntInv = *std::get<1>(preCompute);

    // Pre-compute and store some terms.
    const auto nAbs = std::abs(n);
    const auto lMax = d.MaxDegree();
    const auto cos = std::cos(theta);
    const auto SinHalf = std::sin(0.5 * theta);
    const auto CosHalf = std::cos(0.5 * theta);
    const auto atLeft = SinHalf < std::numeric_limits<Real>::min();
    const auto atRight = CosHalf < std::numeric_limits<Real>::min();
    const auto logSinHalf = atLeft ? static_cast<Real>(0) : std::log(SinHalf);
    const auto logCosHalf = atRight ? static_cast<Real>(0) : std::log(CosHalf);

    // Set the values for l == |n|
    {
      const auto l = nAbs;
      auto m = d(l).MinOrder();
      auto iter = d(l).begin();
      auto finish = d(l).end();
      if (n >= 0) {
        while (iter != finish) {
          *iter++ = WignerMaxUpperIndex(l, m++, logSinHalf, logCosHalf, atLeft,
                                        atRight);
        }
      } else {
        while (iter != finish) {
          *iter++ = WignerMinUpperIndex(l, m++, logSinHalf, logCosHalf, atLeft,
                                        atRight);
        }
      }
    }

    // Set the values for l == n+1 if needed.
    if (nAbs < lMax) {
      const auto l = nAbs + 1;
      const auto mMin = d(l).MinOrder();
      const auto mMax = d(l).MaxOrder();
      auto m = mMin;

      // Set iterators
      auto iterMinusOne = d(l - 1).begin();
      auto finishMinusOne = d(l - 1).end();
      auto iter = d(l).begin();

      // Add in value at m == -l if needed.
      if constexpr (std::same_as<MRange, All>) {
        if (l <= mMax) {
          *iter++ =
              WignerMinOrder(l, n, logSinHalf, logCosHalf, atLeft, atRight);
          m++;
        }
      }

      // Add in interior orders using one-term recursion.
      {
        const auto alpha = (2 * l - 1) * l * cos * sqrtIntInv[l + nAbs];
        const auto beta = (n < 0 ? -1 : 1) * (2 * l - 1) * sqrtIntInv[l + nAbs];

        while (iterMinusOne != finishMinusOne) {
          auto f1 = (alpha - beta * m) * sqrtIntInv[l - m] * sqrtIntInv[l + m];
          *iter++ = f1 * *iterMinusOne++;
          m++;
        }
      }

      // Add in value at m == l if needed
      if (l <= mMax) {
        *iter++ = WignerMaxOrder(l, n, logSinHalf, logCosHalf, atLeft, atRight);
      }
    }

    // Do the remaining degrees
    for (auto l = nAbs + 2; l <= lMax; l++) {
      const auto mMin = d(l).MinOrder();
      const auto mMax = d(l).MaxOrder();
      auto m = mMin;

      // Set iterators.
      auto iterMinusTwo = d(l - 2).begin();
      auto finishMinusTwo = d(l - 2).end();
      auto iterMinusOne = d(l - 1).begin();
      auto iter = d(l).begin();

      // Add in lower boundary terms if still growing.
      if constexpr (std::same_as<MRange, All>) {
        if (l <= mMax) {
          {
            // Add in the m == -l term.
            *iter++ =
                WignerMinOrder(l, n, logSinHalf, logCosHalf, atLeft, atRight);
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

        // Add in the lower boundary term at the critical degree
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
          auto denom = sqrtIntInv[l - m] * sqrtIntInv[l + m];
          auto f1 = (alpha - beta * m) * denom;
          auto f2 = gamma * sqrtInt[l - 1 - m] * sqrtInt[l - 1 + m] * denom;
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
        *iter++ = WignerMaxOrder(l, n, logSinHalf, logCosHalf, atLeft, atRight);
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
        auto start = d(l).begin();
        auto finish = d(l).end();

        std::transform(start, finish, start, [l, factor](auto p) {
          return factor * std::sqrt(static_cast<Real>(2 * l + 1)) * p;
        });
      }
    }
  }
};  // namespace GSHTrans

}  // namespace GSHTrans

#endif  // GSH_TRANS_WIGNER_GUARD_H
