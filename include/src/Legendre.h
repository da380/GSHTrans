#ifndef Legendre_GUARD_H
#define Legendre_GUARD_H

#include <cassert>
#include <cmath>
#include <concepts>
#include <iterator>
#include <limits>
#include <numbers>
#include <utility>
#include <vector>

namespace GSHT {

template <std::floating_point Float>
class LegendreValues {
 public:
  // Define member types.
  using value_type = Float;
  using iterator = typename std::vector<Float>::iterator;

  // Constructors.
  LegendreValues() = delete;
  LegendreValues(Float theta, int L, int M);


  // Basic data functions.
  int MaxDegree() { return L; }
  int MaxOrder() { return M; }
  Float Angle() { return theta; }

  // Functions to return iterators to the values.
  iterator begin() { return data.begin(); }
  iterator end() { return data.end(); }
  iterator begin(int l) { return std::next(begin(), CountBeforeDegree(l)); }
  iterator end(int l) { return std::next(begin(),CountToDegree(l)); }

  // Function to return a specific value.
  Float operator() (int l, int m)
  {
    auto p = *std::next(begin(l), m < 0 ? -m  : m);
    return m < 0 ? Sign(m) * p : p;
  }

  
 private:
  int L;
  int M;
  Float theta;
  std::vector<Float> data;
  std::vector<std::pair<int,int>> index;

  int CountToDegree(int l) {
    return l <= M ? l + 1 + (l * (l + 1)) / 2
                  : M + 1 + (M * (M + 1)) / 2 + (l - M) * (M + 1);
  }

  int CountBeforeDegree(int l) { return CountToDegree(l - 1); }

  Float Sign(int m) { return m % 2 ? -1.0 : 1.0; }
  
};

template <std::floating_point Float>
LegendreValues<Float>::LegendreValues(Float theta, int L, int M)
    : theta{theta}, L{L}, M{M} {
  // Check the degree is non-negative.
  assert(L >= 0);

  // Check the maximum order is in range.
  assert(M >= 0 && M <= L);

  // Allocate space to store the values.
  data.reserve(CountToDegree(L));


  int count = 0;
  for (int l = 0; l <= L; l++) {
    for (int m = 0; m <= std::min(l, M); m++) {
      data.push_back(count++);
    }
  }

  /*

  // Deal with degree 0.
  {
    constexpr Float f = 0.5 * std::numbers::inv_sqrtpi_v<Float>;
    data.push_back(f);
  }

  // pre-compute trigonometric terms.
  Float sin = std::sin(theta);
  Float cos = std::cos(theta);

  // Deal with degree 1.
  if (L > 0) {
    constexpr Float f =
        0.5 * std::numbers::sqrt3_v<Float> * std::numbers::inv_sqrtpi_v<Float>;
    constexpr Float g = -f / std::numbers::sqrt2_v<Float>;
    data.push_back(f * cos);
    data.push_back(g * sin);
  }

  // Deal with the remaining terms using recursion.
  if (L > 1) {
  }

  */
}

template <std::floating_point Float>
class LegendreIterator {
 public:
  // Set up some type aliases.
  using iterator_category = std::input_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = std::vector<Float>;
  using reference = value_type&;

  // Define the constructors.
  LegendreIterator() = delete;
  LegendreIterator(Float theta) : l{-1}, theta{theta} {
    Float min = std::numeric_limits<Float>::min();
    Float eps = std::numeric_limits<Float>::epsilon();
  }

  // Return the degree.
  int degree() const { return l; }

  // Return the value at a given order (positive or negative).
  Float operator()(int m) const {
    if (m >= 0) {
      return current[m];
    } else {
      return Sign(-m) * current[-m];
    }
  }

  // Rerturn reference to the current values.
  reference operator*() { return current; }

  // Prefix increment.
  LegendreIterator& operator++() {
    // Increment the degree.
    l++;

    // Deal with l == 0 as a special case
    if (l == 0) {
      constexpr Float f = 0.5 * std::numbers::inv_sqrtpi_v<Float>;
      current.push_back(f);
      previous.push_back(0.0);
      return *this;
    }

    // Pre-compute some terms.
    auto sin = std::sin(theta);
    auto cos = std::cos(theta);
    auto Fl = static_cast<Float>(l);

    // Increment orders m < l using recursion.
    Float f1Pre = std::sqrt(4 * Fl * Fl - 1);
    Float f2Pre = l > 1 ? 1.0 / std::sqrt(4 * (Fl - 1) * (Fl - 1) - 1) : 0.0;
    for (int m = 0; m < l; m++) {
      auto Fm = static_cast<Float>(m);
      Float f1 = f1Pre / std::sqrt(Fl * Fl - Fm * Fm);
      Float f2 = std::sqrt(((Fl - 1) * (Fl - 1) - Fm * Fm)) * f2Pre;
      previous[m] = f1 * (cos * current[m] - f2 * previous[m]);
      std::swap(previous[m], current[m]);
    }

    // Add on a new element at highest order
    Float f = Sign(l) * 0.5 * std::numbers::inv_sqrtpi_v<Float> *
              std::sqrt(2 * Fl + 1);
    if (sin > std::numeric_limits<Float>::min()) {
      f *= std::exp(0.5 * std::lgamma(2 * Fl + 1) -
                    Fl * std::numbers::ln2_v<Float> - std::lgamma(Fl + 1) +
                    Fl * std::log(sin));
    } else {
      f = 0.0;
    }
    current.push_back(f);
    previous.push_back(0.0);

    return *this;
  }

  // Postfix increment.
  LegendreIterator operator++(int) {
    LegendreIterator temp = *this;
    ++*this;
    return temp;
  }

  // Increment a given numbers of steps.
  LegendreIterator& operator+(int steps) {
    for (int i = 0; i <= steps; i++) {
      ++*this;
    }
    return *this;
  }

 private:
  // Store the current degree.
  int l;

  // Store the maximum order
  int M;

  // Store the angle and trig function values
  Float theta;

  // Vectors to store the function values
  std::vector<Float> current;
  std::vector<Float> previous;

  // Function to return (-1)^m
  Float Sign(int m) const { return m % 2 ? -1.0 : 1.0; }
};
}  // namespace GSHT

#endif  // Legendre_GUARD_H
