#ifndef GHS_TRANS_LEGNEDRE_GUARD_H
#define GHS_TRANS_LEGNEDRE_GUARD_H

#include <concepts>
#include <execution>

#include "Wigner.h"

namespace GSHTrans {

template <std::floating_point Float, IndexRange Range>
class LegendreArray {
 public:
  // Define member types.
  using value_type = Float;
  using iterator = typename WignerArrayN<Float, Range>::iterator;
  using const_iterator = typename WignerArrayN<Float, Range>::const_iterator;
  using difference_type = WignerArrayN<Float, Range>::difference_type;

  LegendreArray(int lMax, int mMax, Float theta,
                Normalisation norm = Normalisation::Ortho)
      : d{WignerArrayN<Float, Range>(lMax, mMax, 0, theta, norm)} {}

  // Geters for basic data.
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
  WignerArrayN<Float, Range> d;
};

// Simple legendre function
template <std::floating_point Float>
Float Legendre(int l, int m, Float theta,
               Normalisation norm = Normalisation::Ortho) {
  return Wigner(l, m, 0, theta, norm);
}

}  // namespace GSHTrans

#endif  // GHS_TRANS_LEGNEDRE_GUARD_H
