#ifndef GHS_TRANS_LEGNEDRE_GUARD_H
#define GHS_TRANS_LEGNEDRE_GUARD_H

#include <concepts>
#include <execution>

#include "Wigner.h"

namespace GSHTrans {

template <std::floating_point Float, OrderRange Range>
class Legendre {
 public:
  // Define member types.
  using value_type = Float;
  using iterator = typename std::vector<Float>::iterator;
  using const_iterator = typename std::vector<Float>::const_iterator;
  using difference_type = std::vector<Float>::difference_type;

  // Constructor with default policy.
  Legendre(int L, int M, Float theta, Normalisation norm)
      : d{WignerN<Float, Range>(L, M, 0, theta, norm)} {}

  // Geters for basic data.
  Float Angle() const { return d.Angle(); }
  int MaxDegree() const { return d.MaxDegree(); }
  int MaxOrder() const { return d.MaxOrder(); }

  // Iterators that point to the start of the data.
  iterator begin() { return d.begin(); }
  const_iterator cbegin() const { return d.cbegin(); }

  // Iterators that point to the end of the data.
  iterator end() { return d.end(); }
  const_iterator cend() const { return d.cend(); }

  // Iterators that point to the start of degree l.
  iterator begin(int l) { return d.begin(l); }
  const_iterator cbegin(int l) const { return d.cbegin(l); }

  // Iterators that point to the end of degree l.
  iterator end(int l) { return d.end(l); }
  const_iterator cend(int l) const { return d.cend(l); }

  // Returns value for given degree and order when all orders are stored.
  Float operator()(int l, int m) const { return d(l, m); }

 private:
  WignerN<Float, Range> d;
};

// Simple legendre function
template <std::floating_point Float>
Float legendre(int l, int m, Float theta,
               Normalisation norm = Normalisation::Ortho) {
  return wigner(l, m, 0, theta, norm);
}

}  // namespace GSHTrans

#endif  // GHS_TRANS_LEGNEDRE_GUARD_H
