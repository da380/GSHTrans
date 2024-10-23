#ifndef GSH_TRANS_CANONICAL_COMPONENT_FIELD_UNARY_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_FIELD_UNARY_GUARD_H

#include <complex>
#include <concepts>
#include <cstddef>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../Utility.h"
#include "CanonicalComponentFieldBase.h"

namespace GSHTrans {

// Forward declare the classes.

template <std::ptrdiff_t N, typename Derived>
class CanonicalComponentFieldConj;

template <std::ptrdiff_t N, typename Derived>
class CanonicalComponentFieldReal;

template <std::ptrdiff_t N, typename Derived>
class CanonicalComponentFieldImag;

template <std::ptrdiff_t N, typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class CanonicalComponentFieldUnary;

template <std::ptrdiff_t N, typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar,
                          typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar,
                           typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class CanonicalComponentFieldUnaryWithScalar;

// Set the traits
namespace Internal {

template <std::ptrdiff_t N, typename Derived>
struct Traits<CanonicalComponentFieldConj<N, Derived>> {
  using Int = typename Derived::Int;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Writeable = std::false_type;
};

template <std::ptrdiff_t N, typename Derived>
struct Traits<CanonicalComponentFieldReal<N, Derived>> {
  using Int = typename Derived::Int;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using Scalar = typename Derived::Real;
  using Value = RealValued;
  using Writeable = std::false_type;
};

template <std::ptrdiff_t N, typename Derived>
struct Traits<CanonicalComponentFieldImag<N, Derived>> {
  using Int = typename Derived::Int;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using Scalar = typename Derived::Real;
  using Value = RealValued;
  using Writeable = std::false_type;
};

template <std::ptrdiff_t N, typename Derived, typename Function>
struct Traits<CanonicalComponentFieldUnary<N, Derived, Function>> {
  using Int = typename Derived::Int;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Writeable = std::false_type;
};

template <std::ptrdiff_t N, typename Derived, typename Function>
struct Traits<CanonicalComponentFieldUnaryWithScalar<N, Derived, Function>> {
  using Int = typename Derived::Int;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Writeable = std::false_type;
};

}  // namespace Internal

// Class for the conjugate operation.
template <std::ptrdiff_t _N, typename Derived>
class CanonicalComponentFieldConj
    : public CanonicalComponentFieldBase<
          _N, CanonicalComponentFieldConj<_N, Derived>> {
 public:
  using Int =
      typename Internal::Traits<CanonicalComponentFieldConj<_N, Derived>>::Int;
  using Real =
      typename Internal::Traits<CanonicalComponentFieldConj<_N, Derived>>::Real;
  using Complex = typename Internal::Traits<
      CanonicalComponentFieldConj<_N, Derived>>::Complex;
  using Scalar = typename Internal::Traits<
      CanonicalComponentFieldConj<_N, Derived>>::Scalar;
  using Value = typename Internal::Traits<
      CanonicalComponentFieldConj<_N, Derived>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldConj<_N, Derived>>::Writeable;

  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    if constexpr (std::same_as<Value, RealValued>) {
      return _u[iTheta, iPhi];
    } else {
      return std::conj(_u[iTheta, iPhi]);
    }
  }

  // Constructors.
  CanonicalComponentFieldConj() = delete;
  CanonicalComponentFieldConj(const CanonicalComponentFieldBase<_N, Derived>& u)
      : _u{u} {}

  CanonicalComponentFieldConj(const CanonicalComponentFieldConj&) = default;
  CanonicalComponentFieldConj(CanonicalComponentFieldConj&&) = default;

  // Assignment.
  CanonicalComponentFieldConj& operator=(CanonicalComponentFieldConj&) =
      default;
  CanonicalComponentFieldConj& operator=(CanonicalComponentFieldConj&&) =
      default;

 private:
  const CanonicalComponentFieldBase<_N, Derived>& _u;
};

// Class for the real operation.
template <std::ptrdiff_t _N, typename Derived>
class CanonicalComponentFieldReal
    : public CanonicalComponentFieldBase<
          _N, CanonicalComponentFieldReal<_N, Derived>> {
 public:
  using Int =
      typename Internal::Traits<CanonicalComponentFieldReal<_N, Derived>>::Int;
  using Real =
      typename Internal::Traits<CanonicalComponentFieldReal<_N, Derived>>::Real;
  using Complex = typename Internal::Traits<
      CanonicalComponentFieldReal<_N, Derived>>::Complex;
  using Scalar = typename Internal::Traits<
      CanonicalComponentFieldReal<_N, Derived>>::Scalar;
  using Value = typename Internal::Traits<
      CanonicalComponentFieldReal<_N, Derived>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldReal<_N, Derived>>::Writeable;

  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    if constexpr (std::same_as<typename Derived::Value, RealValued>) {
      return _u[iTheta, iPhi];
    } else {
      return std::real(_u[iTheta, iPhi]);
    }
  }

  // Constructors.
  CanonicalComponentFieldReal() = delete;
  CanonicalComponentFieldReal(const CanonicalComponentFieldBase<_N, Derived>& u)
      : _u{u} {}

  CanonicalComponentFieldReal(const CanonicalComponentFieldReal&) = default;
  CanonicalComponentFieldReal(CanonicalComponentFieldReal&&) = default;

  // Assignment.
  CanonicalComponentFieldReal& operator=(CanonicalComponentFieldReal&) =
      default;
  CanonicalComponentFieldReal& operator=(CanonicalComponentFieldReal&&) =
      default;

 private:
  const CanonicalComponentFieldBase<_N, Derived>& _u;
};

// Class for the imag operation.
template <std::ptrdiff_t _N, typename Derived>
class CanonicalComponentFieldImag
    : public CanonicalComponentFieldBase<
          _N, CanonicalComponentFieldImag<_N, Derived>> {
 public:
  using Int =
      typename Internal::Traits<CanonicalComponentFieldImag<_N, Derived>>::Int;
  using Real =
      typename Internal::Traits<CanonicalComponentFieldImag<_N, Derived>>::Real;
  using Complex = typename Internal::Traits<
      CanonicalComponentFieldImag<_N, Derived>>::Complex;
  using Scalar = typename Internal::Traits<
      CanonicalComponentFieldImag<_N, Derived>>::Scalar;
  using Value = typename Internal::Traits<
      CanonicalComponentFieldImag<_N, Derived>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldImag<_N, Derived>>::Writeable;

  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    if constexpr (std::same_as<typename Derived::Value, RealValued>) {
      return static_cast<Scalar>(0);
    } else {
      return std::imag(_u[iTheta, iPhi]);
    }
  }

  // Constructors.
  CanonicalComponentFieldImag() = delete;
  CanonicalComponentFieldImag(const CanonicalComponentFieldBase<_N, Derived>& u)
      : _u{u} {}

  CanonicalComponentFieldImag(const CanonicalComponentFieldImag&) = default;
  CanonicalComponentFieldImag(CanonicalComponentFieldImag&&) = default;

  // Assignment.
  CanonicalComponentFieldImag& operator=(CanonicalComponentFieldImag&) =
      default;
  CanonicalComponentFieldImag& operator=(CanonicalComponentFieldImag&&) =
      default;

 private:
  const CanonicalComponentFieldBase<_N, Derived>& _u;
};

// Class for unary transformation.
template <std::ptrdiff_t _N, typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class CanonicalComponentFieldUnary
    : public CanonicalComponentFieldBase<
          _N, CanonicalComponentFieldUnary<_N, Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      CanonicalComponentFieldUnary<_N, Derived, Function>>::Int;
  using Real = typename Internal::Traits<
      CanonicalComponentFieldUnary<_N, Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      CanonicalComponentFieldUnary<_N, Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      CanonicalComponentFieldUnary<_N, Derived, Function>>::Scalar;
  using Value = typename Internal::Traits<
      CanonicalComponentFieldUnary<_N, Derived, Function>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldUnary<_N, Derived, Function>>::Writeable;

  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u[iTheta, iPhi]);
  }

  // Constructors.
  CanonicalComponentFieldUnary() = delete;
  CanonicalComponentFieldUnary(
      const CanonicalComponentFieldBase<_N, Derived>& u, Function&& f)
      : _u{u}, _f{f} {}

  CanonicalComponentFieldUnary(const CanonicalComponentFieldUnary&) = default;
  CanonicalComponentFieldUnary(CanonicalComponentFieldUnary&&) = default;

  // Assignment.
  CanonicalComponentFieldUnary& operator=(CanonicalComponentFieldUnary&) =
      default;
  CanonicalComponentFieldUnary& operator=(CanonicalComponentFieldUnary&&) =
      default;

 private:
  const CanonicalComponentFieldBase<_N, Derived>& _u;
  Function& _f;
};

// Class for unary transformation with a scalar parameter.
template <std::ptrdiff_t _N, typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar,
                          typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar,
                           typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class CanonicalComponentFieldUnaryWithScalar
    : public CanonicalComponentFieldBase<
          _N, CanonicalComponentFieldUnaryWithScalar<_N, Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<_N, Derived, Function>>::Int;
  using Real = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<_N, Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<_N, Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<_N, Derived, Function>>::Scalar;
  using Value = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<_N, Derived, Function>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<_N, Derived, Function>>::Writeable;

  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u[iTheta, iPhi], _s);
  }

  // Constructors.
  CanonicalComponentFieldUnaryWithScalar() = delete;
  CanonicalComponentFieldUnaryWithScalar(
      const CanonicalComponentFieldBase<_N, Derived>& u, Function&& f, Scalar s)
      : _u{u}, _f{f}, _s{s} {}

  CanonicalComponentFieldUnaryWithScalar(
      const CanonicalComponentFieldUnaryWithScalar&) = default;
  CanonicalComponentFieldUnaryWithScalar(
      CanonicalComponentFieldUnaryWithScalar&&) = default;

  // Assignment.
  CanonicalComponentFieldUnaryWithScalar& operator=(
      CanonicalComponentFieldUnaryWithScalar&) = default;
  CanonicalComponentFieldUnaryWithScalar& operator=(
      CanonicalComponentFieldUnaryWithScalar&&) = default;

 private:
  const CanonicalComponentFieldBase<_N, Derived>& _u;
  const Function& _f;
  Scalar _s;
};

}  // namespace GSHTrans

#endif