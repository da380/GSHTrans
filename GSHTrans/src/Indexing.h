#ifndef GSH_TRANS_INDEXING_GUARD_H
#define GSH_TRANS_INDEXING_GUARD_H

namespace GSHTrans {

class IntegerRange {
 private:
  int _last;
  int _iter;

 public:
  IntegerRange(int first, int last) : _last{last}, _iter{first} {}
  explicit IntegerRange(int last) : IntegerRange(0, last) {}

  // Iterable functions
  const IntegerRange& begin() const { return *this; }
  const IntegerRange& end() const { return *this; }

  // Iterator functions
  bool operator!=(const IntegerRange&) const { return _iter < _last; }
  void operator++() { ++_iter; }
  int operator*() const { return _iter; }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_INDEXING_GUARD_H
