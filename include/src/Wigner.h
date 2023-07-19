#ifndef WIGNER_GUARD_H
#define WIGNER_GUARD_H

#include <cassert>
#include <cmath>
#include <concepts>
#include <limits>
#include <numbers>
#include <tuple>
#include <vector>

namespace GSHT {

template <std::floating_point Float>
class WignerValues {
 public:
  // Define member types.
  using value_type = Float;
  //  using iterator = typename std::vector<Float>::iterator;
  using iterator = typename std::vector<std::tuple<int, int, int> >::iterator;

  // Constructors.
  WignerValues() = delete;
  WignerValues(int L, int M, int N, Float theta);

  // Geters.
  Float Angle() const { return theta; }
  int MaxDegree() const { return L; }
  int MaxOrder() const { return M; }
  int MaxUpperIndex() const { return N; }

  // Functions to return iterators to the data.
  iterator begin() { return index.begin(); }
  iterator end() { return index.end(); }
  iterator beginForUpperIndex(const int n) {
    return std::next(begin(), ElementsIncludingUpperIndex(n - 1));
  }
  iterator endForUpperIndex(const int n) {
    return std::next(begin(), ElementsIncludingUpperIndex(n));
  }
  iterator beginForDegreeAtUpperIndex(const int l, const int n) {
    return std::next(begin(), ElementsIncludingDegreeAtUpperIndex(l - 1, n));
  }
  iterator endForDegreeAtUpperIndex(const int l, const int n) {
    return std::next(begin(), ElementsIncludingDegreeAtUpperIndex(l, n));
  }

 private:
  Float theta;
  int L;
  int M;
  int N;
  std::vector<Float> data;
  std::vector<std::tuple<int, int, int> > index;

  constexpr int ElementsIncludingUpperIndex(int n) const {
    return (n + 1) * ((M + 1) * (M + 1) + (L - M) * (2 * M + 1)) -
           (n * (n + 1) * (2 * n + 1)) / 6;
  }

  constexpr int ElementsIncludingDegreeAtUpperIndex(int l, int n) const {
    int i = ElementsIncludingUpperIndex(n - 1);
    if (l >= n) {
      i -= n * n;
      if (l <= M) {
        i += (l + 1) * (l + 1);
      } else {
        i += (M + 1) * (M + 1) + (l - M) * (2 * M + 1);
      }
    }
    return i;
  }

  constexpr int ElementsBeforeDegreeAtUpperIndex(int l, int n) const {
    return ElementsIncludingDegreeAtUpperIndex(l - 1, n);
  }

  Float Sign(int m) { return m % 2 ? -1.0 : 1.0; }
};

template <std::floating_point Float>
WignerValues<Float>::WignerValues(const int L, const int M, const int N,
                                  const Float theta)
    : L{L}, M{M}, N{N}, theta{theta} {
  // Check the maximum degree is non-negative.
  assert(L >= 0);

  // Check the maximum order is in range.
  assert(M >= 0 && M <= L);

  // Check the upper index is in range.
  assert(N >= 0 && N <= M);

  for (int n = 0; n <= N; n++) {
    for (int l = n; l <= L; l++) {
      int mm = std::min(l, M);
      for (int m = -mm; m <= mm; m++) {
        index.push_back({n, l, m});
      }
    }
  }
}


}  // namespace GSHT

#endif  // WIGNER_GUARD_H
