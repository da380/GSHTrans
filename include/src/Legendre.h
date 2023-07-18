#ifndef LEGENDRE_GUARD_H
#define LEGENDRE_GUARD_H

#include <cassert>
#include <cmath>
#include <concepts>
#include <limits>
#include <numbers>
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

  // Basic geters.
  Float Angle() const { return theta; }
  int MaxDegree() const { return L; }
  int MaxOrder() const { return M; }

  // Functions to return iterators to the values.
  iterator begin() { return data.begin(); }
  iterator end() { return data.end(); }
  iterator beginForDegree(int l) {
    return std::next(begin(), CountBeforeDegree(l));
  }
  iterator endForDegree(int l) { return std::next(begin(), CountToDegree(l)); }

  // Function to return value at given degree and order.
  // Note that negative orders are obtained using symmetry.
  Float operator()(int l, int m) {
    auto p = *std::next(beginForDegree(l), std::abs(m));
    return m < 0 ? Sign(m) * p : p;
  }

 private:
  Float theta;
  int L;
  int M;
  std::vector<Float> data;

  int CountToDegree(int l) {
    return l <= M ? l + 1 + (l * (l + 1)) / 2
                  : M + 1 + (M * (M + 1)) / 2 + (l - M) * (M + 1);
  }

  int CountBeforeDegree(int l) { return CountToDegree(l - 1); }

  Float Sign(int m) { return m % 2 ? -1.0 : 1.0; }
};

template <std::floating_point Float>
LegendreValues<Float>::LegendreValues(const Float theta, const int L,
                                      const int M)
    : theta{theta}, L{L}, M{M} {
  // Check the degree is non-negative.
  assert(L >= 0);

  // Check the maximum order is in range.
  assert(M >= 0 && M <= L);

  // Allocate space to store the values.
  data.reserve(CountToDegree(L));

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
    if (M > 0) data.push_back(g * sin);
  }

  // Deal with degrees > 1.
  if (L > 1) {
    for (int l = 2; l <= L; l++) {
      Float Fl = static_cast<Float>(l);

      // Apply the recursion relation
      auto minus2 = beginForDegree(l - 2);
      auto minus1 = beginForDegree(l - 1);
      for (int m = 0; m <= std::min(l - 1, M); m++) {
        Float Fm = static_cast<Float>(m);
        Float f1 = std::sqrt((4 * Fl * Fl - 1) / (Fl * Fl - Fm * Fm));
        Float f2 = std::sqrt(((Fl - 1) * (Fl - 1) - Fm * Fm) /
                             (4 * (Fl - 1) * (Fl - 1) - 1));
        data.push_back(f1 * (cos * (*minus1++) - f2 * (*minus2++)));
      }

      // Add in term at maximum order if needed.
      if (l <= M) {
        if (sin > std::numeric_limits<Float>::min()) {
          Float p = Sign(l) * 0.5 * std::numbers::inv_sqrtpi_v<Float> *
                    std::sqrt(2 * Fl + 1) *
                    std::exp(0.5 * std::lgamma(2 * Fl + 1) -
                             Fl * std::numbers::ln2_v<Float> -
                             std::lgamma(Fl + 1) + Fl * std::log(sin));
          data.push_back(p);
        } else {
          data.push_back(static_cast<Float>(0));
        }
      }
    }
  }
}

}  // namespace GSHT

#endif  // LEGENDRE_GUARD_H
