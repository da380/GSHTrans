#ifndef GSH_TRANS_CANONICAL_VECTOR_GUARD_H
#define GSH_TRANS_CANONICAL_VECTOR_GUARD_H

#include <array>
#include <complex>
#include <concepts>

#include "../Concepts.h"

namespace GSHTrans {

class CanonicalVectorBase {
 public:
  using Int = std::ptrdiff_t;
  auto CanonicalIndices() const { return std::ranges::views::iota(-1, 2); }
  void CheckCanonicalIndices(Int alpha) const { assert(std::abs(alpha) <= 1); }
};

template <RealOrComplexFloatingPoint _Scalar>
class CanonicalVector : public CanonicalVectorBase {
 public:
  using Int = std::ptrdiff_t;
  using Scalar = _Scalar;

  using Value =
      std::conditional_t<RealFloatingPoint<Scalar>, RealValued, ComplexValued>;

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

// Concept for vector-valued functions
template <typename Function, typename Real, typename Value>
concept CanonicalVectorValuedFunction =
    requires(Function f, Real theta, Real phi) {
      requires RealFloatingPoint<Real>;
      requires RealOrComplexValued<Value>;
      {
        f(theta, phi)
      } -> std::same_as<CanonicalVector<std::conditional_t<
            std::same_as<RealValued, Value>, Real, std::complex<Real>>>>;
    };

// Type aliases for real and complex vectors.
template <RealFloatingPoint Real>
using RealCanonicalVector = CanonicalVector<Real>;

template <RealFloatingPoint Real>
using ComplexCanonicalVector = CanonicalVector<std::complex<Real>>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_VECTOR_GUARD_H