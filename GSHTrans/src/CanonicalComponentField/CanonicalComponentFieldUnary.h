#ifndef GSH_TRANS_CANONICAL_COMPONENT_FIELD_UNARY_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_FIELD_UNARY_GUARD_H

#include <complex>
#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../Utility.h"
#include "CanonicalComponentFieldBase.h"

namespace GSHTrans {

// Forward declare the classes.

template <typename Derived>
class CanonicalComponentFieldConj;

template <typename Derived>
class CanonicalComponentFieldReal;

template <typename Derived>
class CanonicalComponentFieldImag;

template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class CanonicalComponentFieldUnary;

template <typename Derived, typename Function>
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

template <typename Derived>
struct Traits<CanonicalComponentFieldConj<Derived>> {
  using Int = typename Derived::Int;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Writeable = std::false_type;
};

template <typename Derived>
struct Traits<CanonicalComponentFieldReal<Derived>> {
  using Int = typename Derived::Int;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using Scalar = typename Derived::Real;
  using Value = RealValued;
  using Writeable = std::false_type;
};

template <typename Derived>
struct Traits<CanonicalComponentFieldImag<Derived>> {
  using Int = typename Derived::Int;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using Scalar = typename Derived::Real;
  using Value = RealValued;
  using Writeable = std::false_type;
};

template <typename Derived, typename Function>
struct Traits<CanonicalComponentFieldUnary<Derived, Function>> {
  using Int = typename Derived::Int;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Writeable = std::false_type;
};

template <typename Derived, typename Function>
struct Traits<CanonicalComponentFieldUnaryWithScalar<Derived, Function>> {
  using Int = typename Derived::Int;
  using Real = typename Derived::Real;
  using Complex = typename Derived::Complex;
  using Scalar = typename Derived::Scalar;
  using Value = typename Derived::Value;
  using Writeable = std::false_type;
};

}  // namespace Internal

// Class for the conjugate operation.
template <typename Derived>
class CanonicalComponentFieldConj
    : public CanonicalComponentFieldBase<CanonicalComponentFieldConj<Derived>> {
 public:
  using Int =
      typename Internal::Traits<CanonicalComponentFieldConj<Derived>>::Int;
  using Real =
      typename Internal::Traits<CanonicalComponentFieldConj<Derived>>::Real;
  using Complex =
      typename Internal::Traits<CanonicalComponentFieldConj<Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<CanonicalComponentFieldConj<Derived>>::Scalar;
  using Value =
      typename Internal::Traits<CanonicalComponentFieldConj<Derived>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldConj<Derived>>::Writeable;

  auto UpperIndex() const { return _N; }
  auto& UpperIndex() { return _N; }

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
  CanonicalComponentFieldConj(const CanonicalComponentFieldBase<Derived>& u)
      : _u{u}, _N{u.UpperIndex()} {}

  CanonicalComponentFieldConj(const CanonicalComponentFieldConj&) = default;
  CanonicalComponentFieldConj(CanonicalComponentFieldConj&&) = default;

  // Assignment.
  CanonicalComponentFieldConj& operator=(CanonicalComponentFieldConj&) =
      default;
  CanonicalComponentFieldConj& operator=(CanonicalComponentFieldConj&&) =
      default;

 private:
  const CanonicalComponentFieldBase<Derived>& _u;
  Int _N;
};

// Class for the real operation.
template <typename Derived>
class CanonicalComponentFieldReal
    : public CanonicalComponentFieldBase<CanonicalComponentFieldReal<Derived>> {
 public:
  using Int =
      typename Internal::Traits<CanonicalComponentFieldReal<Derived>>::Int;
  using Real =
      typename Internal::Traits<CanonicalComponentFieldReal<Derived>>::Real;
  using Complex =
      typename Internal::Traits<CanonicalComponentFieldReal<Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<CanonicalComponentFieldReal<Derived>>::Scalar;
  using Value =
      typename Internal::Traits<CanonicalComponentFieldReal<Derived>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldReal<Derived>>::Writeable;

  auto UpperIndex() const { return _N; }
  auto& UpperIndex() { return _N; }

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
  CanonicalComponentFieldReal(const CanonicalComponentFieldBase<Derived>& u)
      : _u{u}, _N{u.UpperIndex()} {}

  CanonicalComponentFieldReal(const CanonicalComponentFieldReal&) = default;
  CanonicalComponentFieldReal(CanonicalComponentFieldReal&&) = default;

  // Assignment.
  CanonicalComponentFieldReal& operator=(CanonicalComponentFieldReal&) =
      default;
  CanonicalComponentFieldReal& operator=(CanonicalComponentFieldReal&&) =
      default;

 private:
  const CanonicalComponentFieldBase<Derived>& _u;
  Int _N;
};

// Class for the imag operation.
template <typename Derived>
class CanonicalComponentFieldImag
    : public CanonicalComponentFieldBase<CanonicalComponentFieldImag<Derived>> {
 public:
  using Int =
      typename Internal::Traits<CanonicalComponentFieldImag<Derived>>::Int;
  using Real =
      typename Internal::Traits<CanonicalComponentFieldImag<Derived>>::Real;
  using Complex =
      typename Internal::Traits<CanonicalComponentFieldImag<Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<CanonicalComponentFieldImag<Derived>>::Scalar;
  using Value =
      typename Internal::Traits<CanonicalComponentFieldImag<Derived>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldImag<Derived>>::Writeable;

  auto UpperIndex() const { return _N; }
  auto& UpperIndex() { return _N; }

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
  CanonicalComponentFieldImag(const CanonicalComponentFieldBase<Derived>& u)
      : _u{u}, _N{u.UpperIndex()} {}

  CanonicalComponentFieldImag(const CanonicalComponentFieldImag&) = default;
  CanonicalComponentFieldImag(CanonicalComponentFieldImag&&) = default;

  // Assignment.
  CanonicalComponentFieldImag& operator=(CanonicalComponentFieldImag&) =
      default;
  CanonicalComponentFieldImag& operator=(CanonicalComponentFieldImag&&) =
      default;

 private:
  const CanonicalComponentFieldBase<Derived>& _u;
  Int _N;
};

// Class for unary transformation.
template <typename Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename Derived::Scalar>,
      typename Derived::Scalar>;
}
class CanonicalComponentFieldUnary
    : public CanonicalComponentFieldBase<
          CanonicalComponentFieldUnary<Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      CanonicalComponentFieldUnary<Derived, Function>>::Int;
  using Real = typename Internal::Traits<
      CanonicalComponentFieldUnary<Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      CanonicalComponentFieldUnary<Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      CanonicalComponentFieldUnary<Derived, Function>>::Scalar;
  using Value = typename Internal::Traits<
      CanonicalComponentFieldUnary<Derived, Function>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldUnary<Derived, Function>>::Writeable;

  auto UpperIndex() const { return _N; }
  auto& UpperIndex() { return _N; }
  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u[iTheta, iPhi]);
  }

  // Constructors.
  CanonicalComponentFieldUnary() = delete;
  CanonicalComponentFieldUnary(const CanonicalComponentFieldBase<Derived>& u,
                               Function&& f)
      : _u{u}, _N{u.UpperIndex()}, _f{f} {}

  CanonicalComponentFieldUnary(const CanonicalComponentFieldUnary&) = default;
  CanonicalComponentFieldUnary(CanonicalComponentFieldUnary&&) = default;

  // Assignment.
  CanonicalComponentFieldUnary& operator=(CanonicalComponentFieldUnary&) =
      default;
  CanonicalComponentFieldUnary& operator=(CanonicalComponentFieldUnary&&) =
      default;

 private:
  const CanonicalComponentFieldBase<Derived>& _u;
  Int _N;
  Function& _f;
};

// Class for unary transformation with a scalar parameter.
template <typename Derived, typename Function>
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
          CanonicalComponentFieldUnaryWithScalar<Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<Derived, Function>>::Int;
  using Real = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<Derived, Function>>::Scalar;
  using Value = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<Derived, Function>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldUnaryWithScalar<Derived, Function>>::Writeable;

  auto UpperIndex() const { return _N; }
  auto& UpperIndex() { return _N; }

  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u[iTheta, iPhi], _s);
  }

  // Constructors.
  CanonicalComponentFieldUnaryWithScalar() = delete;
  CanonicalComponentFieldUnaryWithScalar(
      const CanonicalComponentFieldBase<Derived>& u, Function&& f, Scalar s)
      : _u{u}, _N{u.UpperIndex()}, _f{f}, _s{s} {}

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
  const CanonicalComponentFieldBase<Derived>& _u;
  Int _N;
  const Function& _f;
  Scalar _s;
};

}  // namespace GSHTrans

#endif