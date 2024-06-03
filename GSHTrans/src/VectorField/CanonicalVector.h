#ifndef GSH_TRANS_CANONICAL_VECTOR_GUARD_H
#define GSH_TRANS_CANONICAL_VECTOR_GUARD_H

#include <array>
#include <complex>
#include <concepts>

#include "../Concepts.h"

namespace GSHTrans {

template <typename Derived>
class CanonicalVectorBase {
 public:
  using Int = std::ptrdiff_t;
  auto CanonicalIndices() const { return std::ranges::views::iota(-1, 2); }
  void CheckCanonicalIndices(Int alpha) const { assert(std::abs(alpha) <= 1); }

 private:
  auto& GetDerived() const { return static_cast<const Derived&>(*this); }
  auto& GetDerived() { return static_cast<Derived&>(*this); }
};

template <RealFloatingPoint _Real, RealOrComplexValued _Value>
class CanonicalVector
    : public CanonicalVectorBase<CanonicalVector<_Real, _Value>> {
 public:
  using Int = std::ptrdiff_t;
  using Real = _Real;
  using Value = _Value;
  using Scalar = std::conditional_t<std::same_as<Value, RealValued>, Real,
                                    std::complex<Real>>;

  CanonicalVector() = default;
  constexpr CanonicalVector(Scalar m, Scalar z, Scalar p) : _data{m, z, p} {}

  CanonicalVector(const CanonicalVector&) = default;
  CanonicalVector(CanonicalVector&&) = default;

  CanonicalVector& operator=(const CanonicalVector&) = default;
  CanonicalVector& operator=(CanonicalVector&&) = default;

  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }

  auto operator[](Int alpha) const {
    this->CheckCanonicalIndices(alpha);
    return _data[alpha + 1];
  }

  auto& operator[](Int alpha) {
    this->CheckCanonicalIndices(alpha);
    return _data[alpha + 1];
  }

  friend std::ostream& operator<<(std::ostream& os, const CanonicalVector& v) {
    for (auto alpha : v.CanonicalIndices() | std::ranges::views::take(2)) {
      os << v[alpha] << std::endl;
    }
    os << v[1];
    return os;
  }

 private:
  std::array<Scalar, 3> _data;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_VECTOR_GUARD_H