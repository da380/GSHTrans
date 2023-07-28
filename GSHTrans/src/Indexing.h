#ifndef INDEXING_GUARD_H
#define INDEXING_GUARD_H

namespace GSHTrans {

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

}  // namespace GSHTrans

#endif  // INDEXING_GUARD_H
