#ifndef INDEXING_GUARD_H
#define INDEXING_GUARD_H

namespace GSHT {

class Range {
 private:
  int last;
  int iter;

 public:
  Range(int first, int last) : last{last}, iter{first} {}
  explicit Range(int last) : Range(0, last) {}

  // Iterable functions
  const Range& begin() const { return *this; }
  const Range& end() const { return *this; }

  // Iterator functions
  bool operator!=(const Range&) const { return iter < last; }
  void operator++() { ++iter; }
  int operator*() const { return iter; }
};

auto UpperIndices(int nmin, int nmax) { return Range(nmin, nmax + 1); }
auto UpperIndices(int nmax) { return UpperIndices(-nmax, nmax); }

auto Degrees(int lmin, int lmax) { return Range(lmin, lmax + 1); }
auto Degrees(int lmax) { return Degrees(0, lmax); }

auto Orders(int l) { return Range(-l, l + 1); }
auto Orders(int l, int m) {
  auto mMax = std::min(l, m);
  return Range(-mMax, mMax + 1);
}



}  // namespace GSHT

#endif  // INDEXING_GUARD_H
