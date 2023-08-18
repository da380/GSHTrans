#ifndef GSH_TRANS_WIGNER_GUARD_H
#define GSH_TRANS_WIGNER_GUARD_H

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <execution>
#include <iterator>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <vector>

namespace GSHTrans {

// Define some tag-classes.
struct All {};
struct NonNegative {};

// Define some useful concepts.
template <typename Range>
concept IndexRange =
    std::same_as<Range, All> or std::same_as<Range, NonNegative>;

// Define enum class for normalisation options
enum class Normalisation { FourPi, Ortho };

// Declare some utility functions
constexpr int MinusOneToPower(int);
template <std::floating_point Float>
Float WignerMinOrderAtUpperIndex(int, int, Float, Float, bool, bool);
template <std::floating_point Float>
Float WignerMaxOrderAtUpperIndex(int, int, Float, Float, bool, bool);
template <std::floating_point Float>
Float WignerMaxUpperIndexAtOrder(int, int, Float, Float, bool, bool);
template <std::floating_point Float>
Float WignerMinUpperIndexAtOrder(int, int, Float, Float, bool, bool);

// Declate function that returns Wigner value for given degree, order
// and upper index.
template <std::floating_point Float>
Float Wigner(int, int, int, Float, Normalisation = Normalisation::Ortho);

// Declare function that returns Wigner matrix for given degree,
// with rows index by n and columns by m.
template <std::floating_point Float>
Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> WignerMatrix(
    int, Float, Normalisation = Normalisation::Ortho);

/////////////////////////////////////////////////////////////////////////
//                      WignerArrayN class                             //
/////////////////////////////////////////////////////////////////////////

template <std::floating_point Float, IndexRange Range = All>
class WignerArrayN {
 public:
  // Define member types.
  using value_type = Float;
  using iterator = std::vector<Float>::iterator;
  using const_iterator = std::vector<Float>::const_iterator;
  using difference_type = std::vector<Float>::difference_type;

  // Constructor
  WignerArrayN(int, int, int, Float, Normalisation norm = Normalisation::Ortho);

  // Geters for basic data.
  int MaxDegree() const { return lMax; }
  int MaxOrder() const { return mMax; }
  int UpperIndex() const { return n; }

  // Returns lowest order at given degree.
  int StartingOrder(int l) requires std::same_as<Range, All> {
    assert(n <= l && l <= lMax);
    return -std::min(l, mMax);
  }
  int StartingOrder(int l) requires std::same_as<Range, NonNegative> {
    assert(n <= l && l <= lMax);
    return 0;
  }

  // Iterators that point to the start of the data.
  iterator begin() { return data.begin(); }
  const_iterator cbegin() const { return data.cbegin(); }

  // Iterators that point to the end of the data.
  iterator end() { return data.end(); }
  const_iterator cend() const { return data.cend(); }

  // Iterators that point to the start of degree l.
  iterator begin(int l) { return std::next(begin(), Count(l - 1)); }
  const_iterator cbegin(int l) const {
    return std::next(cbegin(), Count(l - 1));
  }

  // Iterators that point to the end of degree l.
  iterator end(int l) { return std::next(begin(), Count(l)); }
  const_iterator cend(int l) const { return std::next(cbegin(), Count(l)); }

  // Returns value for given degree and order when all orders are stored.
  Float operator()(int l, int m) const requires std::same_as<Range, All> {
    assert(n <= l && l <= lMax);
    auto mMaxAbs = std::min(l, mMax);
    assert(std::abs(m) <= mMaxAbs);
    return *std::next(cbegin(l), mMaxAbs + m);
  }

  // Returns value for given degree and order, m >= 0, when only non-negative
  // orders are stored.
  Float operator()(int l,
                   int m) const requires std::same_as<Range, NonNegative> {
    assert(n <= l && l <= lMax);
    assert(0 <= m && m <= std::min(l, mMax));
    return *std::next(cbegin(l), m);
  }

  void FlipUpperIndex() requires std::same_as<Range, All> {
    if (n == 0) return;
  }

 private:
  // Set the execution policy.
  static constexpr auto policy = std::execution::seq;

  int lMax;  // Maximum degree.
  int mMax;  // Maximum order.
  int n;     // Upper index.

  std::vector<Float> data;  // Vector containing values.

  // Returns the number of values at a given degree when
  // all orders are stored.
  constexpr difference_type Count(
      int l) const requires std::same_as<Range, All> {
    auto nabs = std::abs(n);
    if (l < nabs) return 0;
    if (l <= mMax) return (l + 1) * (l + 1) - nabs * nabs;
    return (mMax + 1) * (mMax + 1) - nabs * nabs + (l - mMax) * (2 * mMax + 1);
  }

  // Returns the number of values at a given degree when
  // only non-negative orders are stored.
  constexpr difference_type Count(
      int l) const requires std::same_as<Range, NonNegative> {
    auto nabs = std::abs(n);
    if (l < nabs) return 0;
    if (l <= mMax) return ((l + 1) * (l + 2)) / 2 - (nabs * (nabs + 1)) / 2;
    return ((mMax + 1) * (mMax + 2)) / 2 - (nabs * (nabs + 1)) / 2 +
           (l - mMax) * (mMax + 1);
  }

  // Returns total number of values.
  constexpr difference_type const Count() { return Count(lMax); }
};

template <std::floating_point Float, IndexRange Range>
WignerArrayN<Float, Range>::WignerArrayN(int lMax, int mMax, int n, Float theta,
                                         Normalisation norm)
    : lMax{lMax}, mMax{mMax}, n{n} {
  // Check the maximum degree is non-negative.
  assert(lMax >= 0);

  // Check the maximum order is in range.
  assert(mMax >= 0 && mMax <= lMax);

  // Initialise the data vector.
  data = std::vector<Float>(Count());

  // Pre-compute and store trigonometric terms.
  auto cos = std::cos(theta);
  auto logSinHalf = std::sin(0.5 * theta);
  auto logCosHalf = std::cos(0.5 * theta);
  auto atLeft = logSinHalf < std::numeric_limits<Float>::min();
  auto atRight = logCosHalf < std::numeric_limits<Float>::min();
  logSinHalf = atLeft ? static_cast<Float>(0) : std::log(logSinHalf);
  logCosHalf = atRight ? static_cast<Float>(0) : std::log(logCosHalf);

  // Pre-compute and store square roots and their inverses up to lMax + mMax.
  std::vector<Float> sqInt(lMax + mMax + 1);
  std::transform(
      policy, sqInt.begin(), sqInt.end(), sqInt.begin(),
      [&](auto &x) { return std::sqrt(static_cast<Float>(&x - &sqInt[0])); });
  std::vector<Float> sqIntInv(lMax + mMax + 1);
  std::transform(
      policy, sqInt.begin(), sqInt.end(), sqIntInv.begin(), [](auto x) {
        return x > static_cast<Float>(0) ? 1 / x : static_cast<Float>(0);
      });

  // Set the values for l == |n|
  const int nabs = std::abs(n);
  {
    auto l = nabs;
    auto mStart = StartingOrder(l);
    auto start = begin(l);
    auto finish = end(l);
    if (n >= 0) {
      std::transform(policy, start, finish, start, [&](auto &p) {
        auto m = mStart + static_cast<int>(&p - &*start);
        return WignerMaxUpperIndexAtOrder(l, m, logSinHalf, logCosHalf, atLeft,
                                          atRight);
      });
    } else {
      std::transform(policy, start, finish, start, [&](auto &p) {
        auto m = mStart + static_cast<int>(&p - &*start);
        return WignerMinUpperIndexAtOrder(l, m, logSinHalf, logCosHalf, atLeft,
                                          atRight);
      });
    }
  }

  // Set the values for l == n+1 if needed.
  if (nabs < lMax) {
    auto l = nabs + 1;
    auto mStart = StartingOrder(l);

    // Set iterators
    auto startMinusOne = begin(l - 1);
    auto finishMinusOne = end(l - 1);
    auto start = begin(l);

    // Add in value at m == -l if needed.
    if constexpr (std::same_as<Range, All>) {
      if (l <= mMax) {
        *start++ = WignerMinOrderAtUpperIndex(l, n, logSinHalf, logCosHalf,
                                              atLeft, atRight);
        // Update the starting order for recursion
        mStart += 1;
      }
    }

    // Add in interior orders using one-term recursion.
    {
      auto alpha = (2 * l - 1) * l * cos * sqIntInv[l + nabs];
      auto beta = (2 * l - 1) * sqIntInv[l + nabs];
      if (n < 0) beta *= -1;
      std::transform(
          policy, startMinusOne, finishMinusOne, start, [&](auto &minus1) {
            auto m = mStart + static_cast<int>(&minus1 - &*startMinusOne);
            auto f1 = (alpha - beta * m) * sqIntInv[l - m] * sqIntInv[l + m];
            return f1 * minus1;
          });
    }

    // Add in value at m == l if needed
    if (l <= mMax) {
      auto mStep = 2 * l - 1;
      if constexpr (std::same_as<Range, NonNegative>) {
        mStep = l;
      }
      *std::next(start, mStep) = WignerMaxOrderAtUpperIndex(
          l, n, logSinHalf, logCosHalf, atLeft, atRight);
    }
  }

  // Now do the remaining degrees.
  for (int l = nabs + 2; l <= lMax; l++) {
    // Starting order within two-term recursion
    auto mStart = StartingOrder(l);

    // Set iterators
    auto startMinusTwo = begin(l - 2);
    auto finishMinusTwo = end(l - 2);
    auto startMinusOne = begin(l - 1);
    auto start = begin(l);

    // Add in lower boundary terms if still growing.
    if constexpr (std::same_as<Range, All>) {
      if (l <= mMax) {
        // Add in the m == -l term.
        *start++ = WignerMinOrderAtUpperIndex(l, n, logSinHalf, logCosHalf,
                                              atLeft, atRight);
        // Now do the m == -l+1 term using one-point recursion.
        {
          auto m = -l + 1;
          auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) *
                    sqIntInv[l - n] * sqIntInv[l + n] * sqIntInv[l - m] *
                    sqIntInv[l + m] / static_cast<Float>(l - 1);
          *start++ = f1 * (*startMinusOne++);
        }
        // Update the starting order for two-term recursion.
        mStart += 2;
      }

      // Add in the lower boundary term at the critical degree
      if (l == mMax + 1) {
        auto m = -mMax;
        auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) * sqIntInv[l - n] *
                  sqIntInv[l + n] * sqIntInv[l - m] * sqIntInv[l + m] /
                  static_cast<Float>(l - 1);
        *start++ = f1 * (*startMinusOne++);
        // Update the starting order for two-term recursion.
        mStart += 1;
      }
    }

    // Apply two-term recusion for the interior orders.
    {
      auto alpha = (2 * l - 1) * l * cos * sqIntInv[l - n] * sqIntInv[l + n];
      auto beta = (2 * l - 1) * n * sqIntInv[l - n] * sqIntInv[l + n] /
                  static_cast<Float>(l - 1);
      auto gamma = l * sqInt[l - 1 - n] * sqInt[l - 1 + n] * sqIntInv[l - n] *
                   sqIntInv[l + n] / static_cast<Float>(l - 1);
      std::transform(
          policy, startMinusTwo, finishMinusTwo, startMinusOne, start,
          [&](auto &minus2, auto &minus1) {
            auto m = mStart + static_cast<int>(&minus2 - &*startMinusTwo);
            auto denom = sqIntInv[l - m] * sqIntInv[l + m];
            auto f1 = (alpha - beta * m) * denom;
            auto f2 = gamma * sqInt[l - 1 - m] * sqInt[l - 1 + m] * denom;
            return f1 * minus1 - f2 * minus2;
          });
    }

    // Add in the upper boundary terms if still growing.
    if (l <= mMax) {
      // Update the iterator
      auto mStep = 2 * l - 3;
      if constexpr (std::same_as<Range, NonNegative>) {
        mStep = l - 1;
      }
      startMinusOne = std::next(startMinusOne, mStep);
      start = std::next(start, mStep);
      // Add in m == l - 1 term using one-point recursion.
      {
        auto m = l - 1;
        auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) * sqIntInv[l - n] *
                  sqIntInv[l + n] * sqIntInv[l - m] * sqIntInv[l + m] /
                  static_cast<Float>(l - 1);
        *start++ = f1 * (*startMinusOne++);
      }
      // Now do m == l.
      *start++ = WignerMaxOrderAtUpperIndex(l, n, logSinHalf, logCosHalf,
                                            atLeft, atRight);
    }

    // Add in the upper boundary term at the crtiical degree.
    if (l == mMax + 1) {
      // Update the iterators.
      auto mStep = 2 * mMax - 1;
      if constexpr (std::same_as<Range, NonNegative>) {
        mStep = mMax;
      }
      startMinusOne = std::next(startMinusOne, mStep);
      start = std::next(start, mStep);
      auto m = mMax;
      auto f1 = (2 * l - 1) * (l * (l - 1) * cos - m * n) * sqIntInv[l - n] *
                sqIntInv[l + n] * sqIntInv[l - m] * sqIntInv[l + m] /
                static_cast<Float>(l - 1);
      *start++ = f1 * (*startMinusOne++);
    }
  }

  if (norm == Normalisation::Ortho) {
    for (int l = 0; l <= lMax; l++) {
      auto start = begin(l);
      auto finish = end(l);
      std::transform(start, finish, start, [l](auto p) {
        return 0.5 * std::sqrt(static_cast<Float>(2 * l + 1)) *
               std::numbers::inv_sqrtpi_v<Float> * p;
      });
    }
  }
}

/////////////////////////////////////////////////////////////////////////
//                   Definition of WignerArray class                   //
/////////////////////////////////////////////////////////////////////////

template <std::floating_point Float, IndexRange MRange = All,
          IndexRange NRange = All>
class WignerArray {
 public:
  using value_type = Float;
  using iterator = WignerArrayN<Float, MRange>::iterator;
  using const_iterator = WignerArrayN<Float, MRange>::const_iterator;
  using difference_type = WignerArrayN<Float, MRange>::difference_type;

  WignerArray(int lMax, int mMax, int nMax, Float theta,
              Normalisation norm = Normalisation::Ortho)
      : lMax{lMax}, mMax{mMax}, nMax{nMax} {
    assert(0 <= lMax);
    assert(0 <= mMax && mMax <= lMax);
    assert(0 <= nMax && nMax <= lMax);
    for (int n = nMin(); n <= nMax; n++) {
      data.push_back(Array(lMax, mMax, n, theta));
    }
  }

  int MaxDegree() const { return lMax; }
  int MaxOrder() const { return mMax; }
  int MaxUpperIndex() const { return nMax; }

  // Return value for given degree, order and
  // upper index.
  Float operator()(int l, int m, int n) const {
    assert(0 <= l && l <= lMax);
    assert(mMin() <= m && m <= mMax);
    assert(nMin() <= n && n <= nMax);
    return data[Index(n)](l, m);
  }

  // Iterators to the data for given upper index.
  iterator begin(int n) { return data[Index(n)].begin(); }
  iterator cbegin(int n) const { return data[Index(n)].cbegin(); }
  iterator end(int n) { return data[Index(n)].end(); }
  iterator cend(int n) const { return data[Index(n)].cend(); }

  // Iterators to the data for given degree and upper index.
  iterator begin(int l, int n) { return data[Index(n)].begin(l); }
  iterator cbegin(int l, int n) const { return data[Index(n)].cbegin(l); }
  iterator end(int l, int n) { return data[Index(n)].end(l); }
  iterator cend(int l, int n) const { return data[Index(n)].cend(l); }

 private:
  int lMax;
  int mMax;
  int nMax;
  using Array = WignerArrayN<Float, MRange>;
  std::vector<Array> data;

  constexpr int mMin() const {
    if constexpr (std::same_as<MRange, All>) {
      return -mMax;
    }
    if constexpr (std::same_as<MRange, NonNegative>) {
      return 0;
    }
  }

  constexpr int nMin() const {
    if constexpr (std::same_as<NRange, All>) {
      return -nMax;
    }
    if constexpr (std::same_as<NRange, NonNegative>) {
      return 0;
    }
  }

  constexpr int Index(int n) const {
    if constexpr (std::same_as<NRange, All>) {
      return nMax + n;
    }
    if constexpr (std::same_as<NRange, NonNegative>) {
      return n;
    }
  }
};

/////////////////////////////////////////////////////////////////////////
//                          WignerArrayLN class                        //
/////////////////////////////////////////////////////////////////////////

template <std::floating_point Float>
class WignerArrayLN {
 public:
  // Define member types.
  using value_type = Float;
  using iterator = std::vector<Float>::iterator;
  using const_iterator = std::vector<Float>::const_iterator;
  using difference_type = std::vector<Float>::difference_type;

  // Constructor.
  WignerArrayLN(int, int, Float, Normalisation norm = Normalisation::Ortho);

  // Basic data functions.
  int Degree() const { return l; }
  int UpperIndex() const { return n; }
  Float Angle() const { return theta; }

  // Iterators to the data.
  iterator begin() { return data.begin(); }
  const_iterator cbegin() const { return data.cbegin(); }
  iterator end() { return data.end(); }
  const_iterator cend() { return data.end(); }

  // Returns value at given order.
  Float operator()(int m) const { return data[l + m]; }

 private:
  int l;                    // Degree.
  int n;                    // Upper index.
  Float theta;              // Angle.
  std::vector<Float> data;  // Stored values.
};

template <std::floating_point Float>
WignerArrayLN<Float>::WignerArrayLN(int l, int n, Float theta,
                                    Normalisation norm)
    : l{l}, n{n}, theta{theta} {
  // Check the inputs
  assert(0 <= l);
  assert(std::abs(n) <= l);

  // Deal with l = 0 separately.
  if (l == 0) {
    data.push_back(1.0);
    if (norm == Normalisation::Ortho) {
      data[0] *= 0.5 * std::numbers::inv_sqrtpi_v<Float>;
    }
    return;
  }

  // Allocate the data vector.
  data = std::vector<Float>(2 * l + 1);

  // Pre-compute some trigonometric terms.
  auto logSinHalf = std::sin(0.5 * theta);
  auto logCosHalf = std::cos(0.5 * theta);
  auto atLeft = logSinHalf < std::numeric_limits<Float>::min();
  auto atRight = logCosHalf < std::numeric_limits<Float>::min();

  // Deal with values at the end points
  if (atLeft) {
    data[l + n] = 1;
    if (norm == Normalisation::Ortho) {
      data[l + n] *= 0.5 * std::sqrt(static_cast<Float>(2 * l + 1)) *
                     std::numbers::inv_sqrtpi_v<Float>;
    }
    return;
  }

  if (atRight) {
    data[l - n] = MinusOneToPower(l + n);
    if (norm == Normalisation::Ortho) {
      data[l - n] *= 0.5 * std::sqrt(static_cast<Float>(2 * l + 1)) *
                     std::numbers::inv_sqrtpi_v<Float>;
    }
    return;
  }

  // Compute remaining trigonometric terms
  logSinHalf = std::log(logSinHalf);
  logCosHalf = std::log(logCosHalf);
  auto cosec = static_cast<Float>(1) / std::sin(theta);
  auto cot = std::cos(theta) * cosec;

  // Pre-compute and store square roots and their inverses up to 2*l+1.
  std::vector<Float> sqInt(2 * l + 1);
  std::transform(sqInt.begin(), sqInt.end(), sqInt.begin(), [&](auto &x) {
    return std::sqrt(static_cast<Float>(&x - &sqInt[0]));
  });
  std::vector<Float> sqIntInv(2 * l + 1);
  std::transform(sqInt.begin(), sqInt.end(), sqIntInv.begin(), [](auto x) {
    return x > static_cast<Float>(0) ? 1 / x : static_cast<Float>(0);
  });

  // Compute the optimal meeting point for the recursion.
  int mOpt = n * std::cos(theta);

  // Set value at minimum order directly.
  data[0] =
      WignerMinOrderAtUpperIndex(l, n, logSinHalf, logCosHalf, false, false);

  //  Upwards recusion.
  Float minusOne = 0;
  for (int m = -l; m < mOpt; m++) {
    Float current = data[l + m];
    data[l + m + 1] = 2 * (n * cosec - m * cot) * sqIntInv[l - m] *
                          sqIntInv[l + m + 1] * current -
                      sqInt[l + m] * sqInt[l - m + 1] * sqIntInv[l - m] *
                          sqIntInv[l + m + 1] * minusOne;
    minusOne = current;
  }

  // Set value at maximum order directly.
  data[2 * l] =
      WignerMaxOrderAtUpperIndex(l, n, logSinHalf, logCosHalf, false, false);

  // Downwards recusion.
  Float plusOne = 0;
  for (int m = l; m > mOpt + 1; m--) {
    Float current = data[l + m];
    data[l + m - 1] = 2 * (n * cosec - m * cot) * sqIntInv[l + m] *
                          sqIntInv[l - m + 1] * current -
                      sqInt[l - m] * sqInt[l + m + 1] * sqIntInv[l + m] *
                          sqIntInv[l - m + 1] * plusOne;
    plusOne = current;
  }

  if (norm == Normalisation::Ortho) {
    std::transform(data.begin(), data.end(), data.begin(), [l](auto p) {
      return 0.5 * std::sqrt(static_cast<Float>(2 * l + 1)) *
             std::numbers::inv_sqrtpi_v<Float> * p;
    });
  }
}

/////////////////////////////////////////////////////////////////////////
//                Definition of some utility functions                 //
/////////////////////////////////////////////////////////////////////////

constexpr int MinusOneToPower(int m) { return m % 2 ? -1 : 1; }

template <std::floating_point Float>
Float WignerMinOrderAtUpperIndex(int l, int n, Float logSinHalf,
                                 Float logCosHalf, bool atLeft, bool atRight) {
  // Check the inputs.
  assert(l >= 0);
  assert(std::abs(n) <= l);

  // Deal with l == 0 case
  if (l == 0) return static_cast<Float>(1);

  // Deal with special case at the left boundary.
  if (atLeft) {
    return n == -l ? static_cast<Float>(1) : static_cast<Float>(0);
  }

  // Deal with special case at the right boundary.
  if (atRight) {
    return n == l ? static_cast<Float>(1) : static_cast<Float>(0);
  }

  // Deal with the general case.
  auto Fl = static_cast<Float>(l);
  auto Fn = static_cast<Float>(n);
  using std::exp;
  using std::lgamma;
  return exp(
      static_cast<Float>(0.5) *
          (lgamma(2 * Fl + 1) - lgamma(Fl - Fn + 1) - lgamma(Fl + Fn + 1)) +
      (Fl + Fn) * logSinHalf + (Fl - Fn) * logCosHalf);
}

template <std::floating_point Float>
Float WignerMaxOrderAtUpperIndex(int l, int n, Float logSinHalf,
                                 Float logCosHalf, bool atLeft, bool atRight) {
  return MinusOneToPower(n + l) * WignerMinOrderAtUpperIndex(l, -n, logSinHalf,
                                                             logCosHalf, atLeft,
                                                             atRight);
}

template <std::floating_point Float>
Float WignerMinUpperIndexAtOrder(int l, int m, Float logSinHalf,
                                 Float logCosHalf, bool atLeft, bool atRight) {
  return WignerMaxOrderAtUpperIndex(l, -m, logSinHalf, logCosHalf, atLeft,
                                    atRight);
}

template <std::floating_point Float>
Float WignerMaxUpperIndexAtOrder(int l, int m, Float logSinHalf,
                                 Float logCosHalf, bool atLeft, bool atRight) {
  return WignerMinOrderAtUpperIndex(l, -m, logSinHalf, logCosHalf, atLeft,
                                    atRight);
}

// Simple function to return Wigner values.
template <std::floating_point Float>
Float Wigner(int l, int m, int n, Float theta, Normalisation norm) {
  return WignerArrayLN(l, n, theta, norm)(m);
}

// Simple function to return values at degree l as a matrix.
template <std::floating_point Float>
Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> WignerMatrix(
    int l, Float theta, Normalisation norm) {
  using Matrix = Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>;
  using Vector = Eigen::Matrix<Float, Eigen::Dynamic, 1>;
  Matrix d(2 * l + 1, 2 * l + 1);
  for (int n = -l; n <= l; n++) {
    auto dn = WignerArrayLN(l, n, theta, norm);
    for (int m = -l; m <= l; m++) {
      d(n + l, m + l) = dn(m);
    }
  }
  return d;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_WIGNER_GUARD_H
