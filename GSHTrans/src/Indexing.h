#ifndef GSH_TRANS_INDEXING_GUARD_H
#define GSH_TRANS_INDEXING_GUARD_H

namespace GSHTrans {

class IntegerRange {
 private:
  int last;
  int iter;

 public:
  IntegerRange(int first, int last) : last{last}, iter{first} {}
  explicit IntegerRange(int last) : IntegerRange(0, last) {}

  // Iterable functions
  const IntegerRange& begin() const { return *this; }
  const IntegerRange& end() const { return *this; }

  // Iterator functions
  bool operator!=(const IntegerRange&) const { return iter < last; }
  void operator++() { ++iter; }
  int operator*() const { return iter; }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_INDEXING_GUARD_H
