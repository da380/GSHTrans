#ifndef GSH_TRANS_CANONICAL_COMPONENT_FIELD_UNARY_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_FIELD_UNARY_GUARD_H

#include <complex>
#include <concepts>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../Utility.h"
#include "CanonicalComponentField.h"
#include "CanonicalComponentFieldBase.h"

namespace GSHTrans {

// Forward declare the classes.

template <typename _Derived>
class CanonicalComponentFieldConj;

template <typename _Derived>
class CanonicalComponentFieldReal;

template <typename _Derived>
class CanonicalComponentFieldImag;

template <typename _Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename _Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename _Derived::Scalar>,
      typename _Derived::Scalar>;
}
class CanonicalComponentFieldUnaryTransformation;

// Set the traits
namespace Internal {

template <typename _Derived>
struct Traits<CanonicalComponentFieldConj<_Derived>> {
  using Int = typename _Derived::Int;
  using Real = typename _Derived::Real;
  using Complex = typename _Derived::Complex;
  using Scalar = typename _Derived::Scalar;
  using Value = typename _Derived::Value;
  using Writeable = std::false_type;
};

template <typename _Derived>
struct Traits<CanonicalComponentFieldReal<_Derived>> {
  using Int = typename _Derived::Int;
  using Real = typename _Derived::Real;
  using Complex = typename _Derived::Complex;
  using Scalar = typename _Derived::Real;
  using Value = RealValued;
  using Writeable = std::false_type;
};

template <typename _Derived>
struct Traits<CanonicalComponentFieldImag<_Derived>> {
  using Int = typename _Derived::Int;
  using Real = typename _Derived::Real;
  using Complex = typename _Derived::Complex;
  using Scalar = typename _Derived::Real;
  using Value = RealValued;
  using Writeable = std::false_type;
};

template <typename _Derived, typename Function>
struct Traits<CanonicalComponentFieldUnaryTransformation<_Derived, Function>> {
  using Int = typename _Derived::Int;
  using Real = typename _Derived::Real;
  using Complex = typename _Derived::Complex;
  using Scalar = typename _Derived::Scalar;
  using Value = typename _Derived::Value;
  using Writeable = std::false_type;
};

}  // namespace Internal

// Class for the conjugate operation.
template <typename _Derived>
class CanonicalComponentFieldConj : public CanonicalComponentFieldBase<
                                        CanonicalComponentFieldConj<_Derived>> {
 public:
  using Int =
      typename Internal::Traits<CanonicalComponentFieldConj<_Derived>>::Int;
  using Real =
      typename Internal::Traits<CanonicalComponentFieldConj<_Derived>>::Real;
  using Complex =
      typename Internal::Traits<CanonicalComponentFieldConj<_Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<CanonicalComponentFieldConj<_Derived>>::Scalar;
  using Value =
      typename Internal::Traits<CanonicalComponentFieldConj<_Derived>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldConj<_Derived>>::Writeable;

  auto UpperIndex() const { return _u.UpperIndex(); }
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
  CanonicalComponentFieldConj(const CanonicalComponentFieldBase<_Derived>& u)
      : _u{u} {}

  CanonicalComponentFieldConj(const CanonicalComponentFieldConj&) = default;
  CanonicalComponentFieldConj(CanonicalComponentFieldConj&&) = default;

  // Assignment.
  CanonicalComponentFieldConj& operator=(CanonicalComponentFieldConj&) =
      default;
  CanonicalComponentFieldConj& operator=(CanonicalComponentFieldConj&&) =
      default;

 private:
  const CanonicalComponentFieldBase<_Derived>& _u;
};

// Class for the real operation.
template <typename _Derived>
class CanonicalComponentFieldReal : public CanonicalComponentFieldBase<
                                        CanonicalComponentFieldReal<_Derived>> {
 public:
  using Int =
      typename Internal::Traits<CanonicalComponentFieldReal<_Derived>>::Int;
  using Real =
      typename Internal::Traits<CanonicalComponentFieldReal<_Derived>>::Real;
  using Complex =
      typename Internal::Traits<CanonicalComponentFieldReal<_Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<CanonicalComponentFieldReal<_Derived>>::Scalar;
  using Value =
      typename Internal::Traits<CanonicalComponentFieldReal<_Derived>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldReal<_Derived>>::Writeable;

  auto UpperIndex() const { return _u.UpperIndex(); }
  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    if constexpr (std::same_as<typename _Derived::Value, RealValued>) {
      return _u[iTheta, iPhi];
    } else {
      return std::real(_u[iTheta, iPhi]);
    }
  }

  // Constructors.
  CanonicalComponentFieldReal() = delete;
  CanonicalComponentFieldReal(const CanonicalComponentFieldBase<_Derived>& u)
      : _u{u} {}

  CanonicalComponentFieldReal(const CanonicalComponentFieldReal&) = default;
  CanonicalComponentFieldReal(CanonicalComponentFieldReal&&) = default;

  // Assignment.
  CanonicalComponentFieldReal& operator=(CanonicalComponentFieldReal&) =
      default;
  CanonicalComponentFieldReal& operator=(CanonicalComponentFieldReal&&) =
      default;

 private:
  const CanonicalComponentFieldBase<_Derived>& _u;
};

// Class for the imag operation.
template <typename _Derived>
class CanonicalComponentFieldImag : public CanonicalComponentFieldBase<
                                        CanonicalComponentFieldImag<_Derived>> {
 public:
  using Int =
      typename Internal::Traits<CanonicalComponentFieldImag<_Derived>>::Int;
  using Real =
      typename Internal::Traits<CanonicalComponentFieldImag<_Derived>>::Real;
  using Complex =
      typename Internal::Traits<CanonicalComponentFieldImag<_Derived>>::Complex;
  using Scalar =
      typename Internal::Traits<CanonicalComponentFieldImag<_Derived>>::Scalar;
  using Value =
      typename Internal::Traits<CanonicalComponentFieldImag<_Derived>>::Value;
  using Writeable = typename Internal::Traits<
      CanonicalComponentFieldImag<_Derived>>::Writeable;

  auto UpperIndex() const { return _u.UpperIndex(); }
  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    if constexpr (std::same_as<typename _Derived::Value, RealValued>) {
      return static_cast<Scalar>(0);
    } else {
      return std::imag(_u[iTheta, iPhi]);
    }
  }

  // Constructors.
  CanonicalComponentFieldImag() = delete;
  CanonicalComponentFieldImag(const CanonicalComponentFieldBase<_Derived>& u)
      : _u{u} {}

  CanonicalComponentFieldImag(const CanonicalComponentFieldImag&) = default;
  CanonicalComponentFieldImag(CanonicalComponentFieldImag&&) = default;

  // Assignment.
  CanonicalComponentFieldImag& operator=(CanonicalComponentFieldImag&) =
      default;
  CanonicalComponentFieldImag& operator=(CanonicalComponentFieldImag&&) =
      default;

 private:
  const CanonicalComponentFieldBase<_Derived>& _u;
};

// Class for unary transformations.
template <typename _Derived, typename Function>
requires requires() {
  requires std::invocable<Function, typename _Derived::Scalar>;
  requires std::convertible_to<
      std::invoke_result_t<Function, typename _Derived::Scalar>,
      typename _Derived::Scalar>;
}
class CanonicalComponentFieldUnaryTransformation
    : public CanonicalComponentFieldBase<
          CanonicalComponentFieldUnaryTransformation<_Derived, Function>> {
 public:
  using Int = typename Internal::Traits<
      CanonicalComponentFieldUnaryTransformation<_Derived, Function>>::Int;
  using Real = typename Internal::Traits<
      CanonicalComponentFieldUnaryTransformation<_Derived, Function>>::Real;
  using Complex = typename Internal::Traits<
      CanonicalComponentFieldUnaryTransformation<_Derived, Function>>::Complex;
  using Scalar = typename Internal::Traits<
      CanonicalComponentFieldUnaryTransformation<_Derived, Function>>::Scalar;
  using Value = typename Internal::Traits<
      CanonicalComponentFieldUnaryTransformation<_Derived, Function>>::Value;
  using Writeable =
      typename Internal::Traits<CanonicalComponentFieldUnaryTransformation<
          _Derived, Function>>::Writeable;

  auto UpperIndex() const { return _u.UpperIndex(); }
  auto& Grid() const { return _u.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _f(_u[iTheta, iPhi]);
  }

  // Constructors.
  CanonicalComponentFieldUnaryTransformation() = delete;
  CanonicalComponentFieldUnaryTransformation(
      const CanonicalComponentFieldBase<_Derived>& u, Function&& f)
      : _u{u}, _f{f} {}

  CanonicalComponentFieldUnaryTransformation(
      const CanonicalComponentFieldUnaryTransformation&) = default;
  CanonicalComponentFieldUnaryTransformation(
      CanonicalComponentFieldUnaryTransformation&&) = default;

  // Assignment.
  CanonicalComponentFieldUnaryTransformation& operator=(
      CanonicalComponentFieldUnaryTransformation&) = default;
  CanonicalComponentFieldUnaryTransformation& operator=(
      CanonicalComponentFieldUnaryTransformation&&) = default;

 private:
  Function& _f;
  const CanonicalComponentFieldBase<_Derived>& _u;
};

}  // namespace GSHTrans

#endif