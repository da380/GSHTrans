#ifndef GSH_TRANS_WIGNER_GUARD_H
#define GSH_TRANS_WIGNER_GUARD_H

#include <omp.h>

#include <cassert>
#include <cmath>
#include <concepts>
#include <limits>
#include <numbers>
#include <ranges>
#include <vector>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

template <RealFloatingPoint Real, OrderIndexRange MRange, IndexRange NRange,
          Normalisation Norm>
class Wigner {
  using Int = std::ptrdiff_t;
  using Vector = std::vector<Real>;

  template <std::ranges::view V>
  class SubView : public std::ranges::view_interface<SubView<V>>,
                  public GSHSubIndices<MRange> {
    using Indices = GSHSubIndices<MRange>;
    using std::ranges::view_interface<SubView<V>>::size;

   public:
    SubView(Int l, Int mMax, V view)
        : Indices(l, mMax), _view{std::move(view)} {}

    auto begin() { return _view.begin(); }
    auto end() { return _view.end(); }

    auto operator()(Int m) const {
      auto i = this->Index(m);
      return this->operator[](i);
    }

    auto &operator()(Int m) {
      auto i = this->Index(m);
      return this->operator[](i);
    }

   private:
    V _view;
  };

  template <std::ranges::view V>
  class View : public std::ranges::view_interface<View<V>>,
               public GSHIndices<MRange> {
    using Indices = GSHIndices<MRange>;
    using std::ranges::view_interface<View<V>>::size;

   public:
    View(Int lMax, Int mMax, Int n, V view)
        : Indices(lMax, mMax, n), _view{view} {}

    auto begin() { return _view.begin(); }
    auto end() { return _view.end(); }

    auto operator()(Int l) const {
      auto view = MakeSubRange(l);
      return SubView(l, this->MaxOrder(), view | std::ranges::views::as_const);
    }

    auto operator()(Int l) {
      auto view = MakeSubRange(l);
      return SubView(l, this->MaxOrder(), view);
    }

   private:
    V _view;

    auto MakeSubRange(Int l) {
      auto [offset, indices] = this->Index(l);
      auto size = indices.size();
      auto start = std::next(begin(), offset);
      auto end = std::next(start, size);
      return std::ranges::subrange(start, end);
    }
  };

  class Arguments {
   public:
    Arguments() = default;

    Arguments(Real theta) {
      constexpr auto half = static_cast<Real>(1) / static_cast<Real>(2);
      _logSinHalf = std::sin(half * theta);
      _logCosHalf = std::cos(half * theta);
      _atLeft = _logSinHalf < std::numeric_limits<Real>::min();
      _atRight = _logCosHalf < std::numeric_limits<Real>::min();
      _logSinHalf = _atLeft ? static_cast<Real>(0) : std::log(_logSinHalf);
      _logCosHalf = _atRight ? static_cast<Real>(0) : std::log(_logCosHalf);
    }

    auto AtLeft() const { return _atLeft; }
    auto AtRight() const { return _atRight; }

    auto LogSinHalf() const { return _logSinHalf; }
    auto LogCosHalf() const { return _logCosHalf; }

   private:
    Real _logSinHalf;
    Real _logCosHalf;
    bool _atLeft;
    bool _atRight;
  };

 public:
  using value_type = Real;

  Wigner() = default;

  // Constructors for a Upper index range with angles.
  template <RealFloatingPointRange RealRange>
  Wigner(Int lMax, Int mMax, Int nMax, RealRange &&thetaRange)
      : _lMax{lMax}, _mMax{mMax}, _nMax{nMax}, _nTheta(thetaRange.size()) {
    AllocateStorage();
    ComputeValues(thetaRange);
  }

  Wigner(Int lMax, Int mMax, Int nMax, Real theta)
      : Wigner(lMax, mMax, nMax, Vector{theta}) {}

  // Return basic information.
  auto MaxDegree() const { return _lMax; }
  auto MaxOrder() const { return _mMax; }
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

  auto operator()(Int n, Int iTheta) {
    auto upperIndices = UpperIndices() | std::ranges::views::filter(
                                             [n](auto np) { return np < n; });
    auto size = GSHIndices<MRange>(_lMax, _mMax, n).size();
    auto offset =
        size * iTheta +
        std::ranges::fold_left(upperIndices, Int{0}, [this](auto acc, auto np) {
          return acc + GSHIndices<MRange>(_lMax, _mMax, np).size();
        }) * NumberOfAngles();
    auto start = std::next(begin(), offset);
    auto finish = std::next(start, size);
    auto view = std::ranges::subrange(start, finish);
    return View(_lMax, _mMax, n, view);
  }

 private:
  Int _lMax;    // Maximum degree.
  Int _mMax;    // Maximum order.
  Int _nMax;    // Maximum upper index.
  Int _nTheta;  // Number of colatitudes.

  // Vector storing the values.
  Vector _data;

  // Compute the necessary storage capacity.
  void AllocateStorage() {
    auto upperIndices = UpperIndices();
    auto size = std::ranges::fold_left(
                    upperIndices, Int{0},
                    [this](auto acc, auto n) {
                      return acc + GSHIndices<MRange>(_lMax, _mMax, n).size();
                    }) *
                NumberOfAngles();
    _data.resize(size);
  }

  template <RealFloatingPointRange RealRange>
  void ComputeValues(RealRange &&thetaRange) {
    auto preCompute = PreCompute();
#pragma omp parallel for collapse(2)
    for (auto n : UpperIndices()) {
      for (auto iTheta : AngleIndices()) {
        auto d = operator()(n, iTheta);
        Compute(n, thetaRange[iTheta], d, preCompute);
      }
    }
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

  void Compute(Int n, Real theta, auto d, const auto preCompute) {
    // Get references to the pre-computed values.
    auto &sqrtInt = *std::get<0>(preCompute);
    auto &sqrtIntInv = *std::get<1>(preCompute);

    // Pre-compute and store some terms.
    const auto nAbs = std::abs(n);
    const auto lMax = d.MaxDegree();

    auto arg = Arguments(theta);
    const auto cos = std::cos(theta);

    // Set the values for l == |n|
    {
      const auto l = nAbs;
      auto m = d(l).MinOrder();
      auto iter = d(l).begin();
      auto finish = d(l).end();
      if (n >= 0) {
        while (iter != finish) {
          *iter++ = WignerMaxUpperIndex(l, m++, arg);
        }
      } else {
        while (iter != finish) {
          *iter++ = WignerMinUpperIndex(l, m++, arg);
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
          *iter++ = WignerMinOrder(l, n, arg);
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
        *iter++ = WignerMaxOrder(l, n, arg);
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
            *iter++ = WignerMinOrder(l, n, arg);
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
        *iter++ = WignerMaxOrder(l, n, arg);
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

  constexpr auto MinusOneToPower(Int m) { return m % 2 ? -1 : 1; }

  auto WignerMinOrder(Int l, Int n, Arguments &arg) {
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
    return exp(half * (lgamma(2 * Fl + 1) - lgamma(Fl - Fn + 1) -
                       lgamma(Fl + Fn + 1)) +
               (Fl + Fn) * arg.LogSinHalf() + (Fl - Fn) * arg.LogCosHalf());
  }

  auto WignerMaxOrder(Int l, Int n, Arguments &arg) {
    return MinusOneToPower(n + l) * WignerMinOrder(l, -n, arg);
  }

  auto WignerMinUpperIndex(Int l, Int m, Arguments &arg) {
    return WignerMaxOrder(l, -m, arg);
  }

  auto WignerMaxUpperIndex(Int l, Int m, Arguments &arg) {
    return WignerMinOrder(l, -m, arg);
  }
};

namespace Testing {

template <RealFloatingPoint Real, OrderIndexRange MRange = All,
          Normalisation Norm = Ortho, IndexRange NRange = Single,
          AngleIndexRange AngleRange = Single,
          MatrixStorage Storage = ColumnMajor>
class Wigner {
  using Int = std::ptrdiff_t;
  using Vector = std::vector<Real>;

  template <std::ranges::view V>
  class SubView : public std::ranges::view_interface<SubView<V>>,
                  public GSHSubIndices<MRange> {
    using Indices = GSHSubIndices<MRange>;
    using std::ranges::view_interface<SubView<V>>::size;

   public:
    SubView(Int l, Int mMax, V view)
        : Indices(l, mMax), _view{std::move(view)} {}

    auto begin() { return _view.begin(); }
    auto end() { return _view.end(); }

    auto operator()(Int m) const {
      auto i = this->Index(m);
      return this->operator[](i);
    }

    auto &operator()(Int m) {
      auto i = this->Index(m);
      return this->operator[](i);
    }

   private:
    V _view;
  };

  template <std::ranges::view V>
  class View : public std::ranges::view_interface<View<V>>,
               public GSHIndices<MRange> {
    using Indices = GSHIndices<MRange>;
    using std::ranges::view_interface<View<V>>::size;

   public:
    View(Int lMax, Int mMax, Int n, V view)
        : Indices(lMax, mMax, n), _view{view} {}

    auto begin() { return _view.begin(); }
    auto end() { return _view.end(); }

    auto operator()(Int l) const {
      auto view = MakeSubRange(l);
      return SubView(l, this->MaxOrder(), view | std::ranges::views::as_const);
    }

    auto operator()(Int l) {
      auto view = MakeSubRange(l);
      return SubView(l, this->MaxOrder(), view);
    }

   private:
    V _view;

    auto MakeSubRange(Int l) {
      auto [offset, indices] = this->Index(l);
      auto size = indices.size();
      auto start = std::next(begin(), offset);
      auto end = std::next(start, size);
      return std::ranges::subrange(start, end);
    }
  };

  class Arguments {
   public:
    Arguments() = default;

    Arguments(Real theta) {
      constexpr auto half = static_cast<Real>(1) / static_cast<Real>(2);
      _logSinHalf = std::sin(half * theta);
      _logCosHalf = std::cos(half * theta);
      _atLeft = _logSinHalf < std::numeric_limits<Real>::min();
      _atRight = _logCosHalf < std::numeric_limits<Real>::min();
      _logSinHalf = _atLeft ? static_cast<Real>(0) : std::log(_logSinHalf);
      _logCosHalf = _atRight ? static_cast<Real>(0) : std::log(_logCosHalf);
    }

    auto AtLeft() const { return _atLeft; }
    auto AtRight() const { return _atRight; }

    auto LogSinHalf() const { return _logSinHalf; }
    auto LogCosHalf() const { return _logCosHalf; }

   private:
    Real _logSinHalf;
    Real _logCosHalf;
    bool _atLeft;
    bool _atRight;
  };

 public:
  Wigner() = default;

  template <RealFloatingPointRange RealRange>
  Wigner(Int lMax, Int mMax, Int nMax, RealRange &&thetaRange)
      : _lMax{lMax}, _mMax{mMax}, _nMax{nMax}, _nTheta(thetaRange.size()) {
    AllocateStorage();
    ComputeValues(thetaRange);
  }

  Wigner(Int lMax, Int mMax, Int nMax, Real theta)
  requires std::same_as<AngleRange, Single>
      : Wigner(lMax, mMax, nMax, std::vector(1, theta)) {}

  // Return basic information.
  auto MaxDegree() const { return _lMax; }
  auto MaxOrder() const { return _mMax; }
  auto NumberOfAngles() const { return _nTheta; }

  auto Degrees() const
  requires std::same_as<NRange, Single>
  {
    return GSHIndices<MRange>(_lMax, _mMax, _nMax).Degrees();
  }

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

  // Return pairs of (n,iTheta) in the storage order.
  auto Indices() const {
    if constexpr (std::same_as<Storage, ColumnMajor>) {
      return std::ranges::views::cartesian_product(UpperIndices(),
                                                   AngleIndices());
    } else {
      return std::ranges::views::cartesian_product(AngleIndices(),
                                                   UpperIndices());
    }
  }

  auto operator()(Int n, Int iTheta)
  requires std::same_as<Storage, ColumnMajor>
  {
    auto upperIndices = UpperIndices() | std::ranges::views::filter(
                                             [n](auto np) { return np < n; });
    auto size = GSHIndices<MRange>(_lMax, _mMax, n).size();
    auto offset =
        size * iTheta +
        std::ranges::fold_left(upperIndices, Int{0}, [this](auto acc, auto np) {
          return acc + GSHIndices<MRange>(_lMax, _mMax, np).size();
        }) * NumberOfAngles();
    auto start = std::next(begin(), offset);
    auto finish = std::next(start, size);
    auto view = std::ranges::subrange(start, finish);
    return View(_lMax, _mMax, n, view);
  }

  auto operator()(Int n, Int iTheta)
  requires std::same_as<Storage, RowMajor>
  {
    auto size = GSHIndices<MRange>(_lMax, _mMax, n).size();
    auto offset = std::ranges::fold_left(
                      UpperIndices(), Int{0},
                      [this](auto acc, auto n) {
                        return acc + GSHIndices<MRange>(_lMax, _mMax, n).size();
                      }) *
                      iTheta +
                  std::ranges::fold_left(
                      UpperIndices() | std::ranges::views::filter(
                                           [n](auto np) { return np < n; }),
                      Int{0}, [this](auto acc, auto n) {
                        return acc + GSHIndices<MRange>(_lMax, _mMax, n).size();
                      });
    auto start = std::next(begin(), offset);
    auto finish = std::next(start, size);
    auto view = std::ranges::subrange(start, finish);
    return View(_lMax, _mMax, n, view);
  }

  auto operator()(Int n)
  requires std::same_as<AngleRange, Single>
  {
    return operator()(n, 0);
  }

  auto operator()(Int iTheta)
  requires std::same_as<NRange, Single>
  {
    return operator()(_nMax, iTheta);
  }

  auto operator()(Int l)
  requires std::same_as<NRange, Single> and std::same_as<AngleRange, Single>
  {
    return operator()(_nMax, 0)(l);
  }

 private:
  Int _lMax;    // Maximum degree.
  Int _mMax;    // Maximum order.
  Int _nMax;    // Maximum upper index.
  Int _nTheta;  // Number of colatitudes.

  // Vector storing the values.
  Vector _data;

  // Compute the necessary storage capacity.
  void AllocateStorage() {
    auto upperIndices = UpperIndices();
    auto size = std::ranges::fold_left(
                    upperIndices, Int{0},
                    [this](auto acc, auto n) {
                      return acc + GSHIndices<MRange>(_lMax, _mMax, n).size();
                    }) *
                NumberOfAngles();
    _data.resize(size);
  }

  template <RealFloatingPointRange RealRange>
  void ComputeValues(RealRange &&thetaRange) {
    auto preCompute = PreCompute();
#pragma omp parallel for
    for (auto [n, iTheta] : Indices()) {
      auto d = operator()(n, iTheta);
      Compute(n, thetaRange[iTheta], d, preCompute);
    }
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

  void Compute(Int n, Real theta, auto d, const auto preCompute) {
    // Get references to the pre-computed values.
    auto &sqrtInt = *std::get<0>(preCompute);
    auto &sqrtIntInv = *std::get<1>(preCompute);

    // Pre-compute and store some terms.
    const auto nAbs = std::abs(n);
    const auto lMax = d.MaxDegree();

    auto arg = Arguments(theta);
    const auto cos = std::cos(theta);

    // Set the values for l == |n|
    {
      const auto l = nAbs;
      auto m = d(l).MinOrder();
      auto iter = d(l).begin();
      auto finish = d(l).end();
      if (n >= 0) {
        while (iter != finish) {
          *iter++ = WignerMaxUpperIndex(l, m++, arg);
        }
      } else {
        while (iter != finish) {
          *iter++ = WignerMinUpperIndex(l, m++, arg);
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
          *iter++ = WignerMinOrder(l, n, arg);
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
        *iter++ = WignerMaxOrder(l, n, arg);
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
            *iter++ = WignerMinOrder(l, n, arg);
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
        *iter++ = WignerMaxOrder(l, n, arg);
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

  constexpr auto MinusOneToPower(Int m) { return m % 2 ? -1 : 1; }

  auto WignerMinOrder(Int l, Int n, Arguments &arg) {
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
    return exp(half * (lgamma(2 * Fl + 1) - lgamma(Fl - Fn + 1) -
                       lgamma(Fl + Fn + 1)) +
               (Fl + Fn) * arg.LogSinHalf() + (Fl - Fn) * arg.LogCosHalf());
  }

  auto WignerMaxOrder(Int l, Int n, Arguments &arg) {
    return MinusOneToPower(n + l) * WignerMinOrder(l, -n, arg);
  }

  auto WignerMinUpperIndex(Int l, Int m, Arguments &arg) {
    return WignerMaxOrder(l, -m, arg);
  }

  auto WignerMaxUpperIndex(Int l, Int m, Arguments &arg) {
    return WignerMinOrder(l, -m, arg);
  }
};

}  // namespace Testing

}  // namespace GSHTrans

#endif  // GSH_TRANS_WIGNER_GUARD_H
