#ifndef Legendre_GUARD_H
#define Legendre_GUARD_H

#include <cassert>
#include <cmath>
#include <concepts>
#include <limits>
#include <numbers>
#include <utility>
#include <vector>

namespace GSHT {

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
  LegendreIterator(Float theta, int mMax = -1)
      : l{-1}, theta{theta}, mMax{mMax} {
    Float min = std::numeric_limits<Float>::min();
    Float eps = std::numeric_limits<Float>::epsilon();
    bool too_low = theta < -min;
    bool too_high = theta > std::numbers::pi_v<Float> * (1 + eps);
    assert(!too_low && !too_high);
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
    for (int m = 0; m < l; m++) {
      auto Fm = static_cast<Float>(m);
      Float f1 = std::sqrt((4 * Fl * Fl - 1) / (Fl * Fl - Fm * Fm));
      Float f2 = std::sqrt(((Fl - 1) * (Fl - 1) - Fm * Fm) /
                           (4 * (Fl - 1) * (Fl - 1) - 1));
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

  // Store the maximum order to be computed
  int mMax;

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
