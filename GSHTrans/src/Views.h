#ifndef GSH_TRANS_VIEWS_GUARD_H
#define GSH_TRANS_VIEWS_GUARD_H

#include "Concepts.h"
#include "Indexing.h"
namespace GSHTrans {

/*-----------------------------------------------------------------/
                    Non constant view classes
/------------------------------------------------------------------*/
template <RealOrComplexFloatingPoint Scalar, OrderIndexRange MRange>
class GSHSubView : public GSHSubIndices<MRange> {
 public:
  using Int = std::ptrdiff_t;
  constexpr GSHSubView(Int l, Int mMax, Scalar* start)
      : GSHSubIndices<MRange>(l, mMax), _start{start} {}

  constexpr auto begin() { return _start; }
  constexpr auto end() { return std::next(begin(), this->Size()); }

  constexpr auto& operator[](Int m) { return _start[this->Index(m)]; }

 private:
  Scalar* _start;
};

template <RealOrComplexFloatingPoint Scalar, OrderIndexRange MRange>
class GSHView : public GSHIndices<MRange> {
 public:
  using Int = std::ptrdiff_t;

  constexpr GSHView(Int lMax, Int mMax, Int n, Scalar* start)
      : GSHIndices<MRange>(lMax, mMax, n), _start{start} {}

  constexpr auto begin() { return _start; }
  constexpr auto end() { return std::next(begin(), this->Size()); }

  constexpr auto operator[](Int l) {
    return GSHSubView<Scalar, MRange>(
        l, this->MaxOrder(), std::next(begin(), this->OffsetForDegree(l)));
  }

  constexpr auto& operator[](Int l, Int m) { return _start[this->Index(l, m)]; }

 private:
  Scalar* _start;
};

/*-----------------------------------------------------------------/
                        constant view classes
/------------------------------------------------------------------*/
template <RealOrComplexFloatingPoint Scalar, OrderIndexRange MRange>
class ConstGSHSubView : public GSHSubIndices<MRange> {
 public:
  using Int = std::ptrdiff_t;
  constexpr ConstGSHSubView(Int l, Int mMax, const Scalar* start)
      : GSHSubIndices<MRange>(l, mMax), _start{start} {}

  constexpr auto begin() const { return _start; }
  constexpr auto end() const { return std::next(begin(), this->Size()); }

  constexpr auto operator[](Int m) const { return _start[this->Index(m)]; }

 private:
  const Scalar* _start;
};

template <RealOrComplexFloatingPoint Scalar, OrderIndexRange MRange>
class ConstGSHView : public GSHIndices<MRange> {
 public:
  using Int = std::ptrdiff_t;

  constexpr ConstGSHView(Int lMax, Int mMax, Int n, const Scalar* start)
      : GSHIndices<MRange>(lMax, mMax, n), _start{start} {}

  constexpr auto begin() const { return _start; }
  constexpr auto end() const { return std::next(begin(), this->Size()); }

  constexpr auto operator[](Int l) const {
    return ConstGSHSubView<Scalar, MRange>(
        l, this->MaxOrder(), std::next(begin(), this->OffsetForDegree(l)));
  }

  constexpr auto operator[](Int l, Int m) { return _start[this->Index(l, m)]; }

 private:
  const Scalar* _start;
};

}  // namespace GSHTrans

#endif