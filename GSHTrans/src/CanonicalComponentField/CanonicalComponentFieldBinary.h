#ifndef GSH_TRANS_CANONICAL_COMPONENT_FIELD_BINARY_GUARD_H
#define GSH_TRANS_CANONICAL_COMPONENT_FIELD_BINARY_GUARD_H

#include <complex>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../Utility.h"
#include "CanonicalComponentFieldBase.h"

namespace GSHTrans {

// Forward declare the classes.

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
requires requires() {
  requires std::same_as<typename Derived0::Scalar, typename Derived1::Scalar>;
  requires(N0 == N1);
}
class CanonicalComponentFieldAdd;

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
requires requires() {
  requires std::same_as<typename Derived0::Scalar, typename Derived1::Scalar>;
  requires(N0 == N1);
}
class CanonicalComponentFieldSubtract;

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
requires requires() {
  requires std::same_as<typename Derived0::Scalar, typename Derived1::Scalar>;
}
class CanonicalComponentFieldMultiply;

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
requires requires() {
  requires std::same_as<typename Derived0::Scalar, typename Derived1::Scalar>;
  requires(N1 == 0);
}
class CanonicalComponentFieldDivide;

// Set up the traits

namespace Internal {

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
struct Traits<CanonicalComponentFieldAdd<N0, Derived0, N1, Derived1>> {
  using Int = typename Derived0::Int;
  using Real = typename Derived0::Real;
  using Complex = typename Derived0::Complex;
  using Scalar = typename Derived0::Scalar;
  using Value = typename Derived0::Value;
  using Writeable = std::false_type;
};

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
struct Traits<CanonicalComponentFieldSubtract<N0, Derived0, N1, Derived1>> {
  using Int = typename Derived0::Int;
  using Real = typename Derived0::Real;
  using Complex = typename Derived0::Complex;
  using Scalar = typename Derived0::Scalar;
  using Value = typename Derived0::Value;
  using Writeable = std::false_type;
};

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
struct Traits<CanonicalComponentFieldMultiply<N0, Derived0, N1, Derived1>> {
  using Int = typename Derived0::Int;
  using Real = typename Derived0::Real;
  using Complex = typename Derived0::Complex;
  using Scalar = typename Derived0::Scalar;
  using Value = typename Derived0::Value;
  using Writeable = std::false_type;
};

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
struct Traits<CanonicalComponentFieldDivide<N0, Derived0, N1, Derived1>> {
  using Int = typename Derived0::Int;
  using Real = typename Derived0::Real;
  using Complex = typename Derived0::Complex;
  using Scalar = typename Derived0::Scalar;
  using Value = typename Derived0::Value;
  using Writeable = std::false_type;
};
}  // namespace Internal

// Add two fields.
template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
requires requires() {
  requires std::same_as<typename Derived0::Scalar, typename Derived1::Scalar>;
  requires(N0 == N1);
}
class CanonicalComponentFieldAdd
    : public CanonicalComponentFieldBase<
          N0, CanonicalComponentFieldAdd<N0, Derived0, N1, Derived1>> {
 public:
  using Int = Internal::Traits<
      CanonicalComponentFieldAdd<N0, Derived0, N1, Derived1>>::Int;
  using Real = Internal::Traits<
      CanonicalComponentFieldAdd<N0, Derived0, N1, Derived1>>::Real;
  using Complex = Internal::Traits<
      CanonicalComponentFieldAdd<N0, Derived0, N1, Derived1>>::Complex;
  using Scalar = Internal::Traits<
      CanonicalComponentFieldAdd<N0, Derived0, N1, Derived1>>::Scalar;
  using Value = Internal::Traits<
      CanonicalComponentFieldAdd<N0, Derived0, N1, Derived1>>::Value;
  using Writeable = Internal::Traits<
      CanonicalComponentFieldAdd<N0, Derived0, N1, Derived1>>::Writeable;

  auto& Grid() const { return _u0.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _u0[iTheta, iPhi] + _u1[iTheta, iPhi];
  }

  // Constructors.
  CanonicalComponentFieldAdd() = delete;
  CanonicalComponentFieldAdd(
      const CanonicalComponentFieldBase<N0, Derived0>& u0,
      const CanonicalComponentFieldBase<N1, Derived1>& u1)
      : _u0{u0}, _u1{u1} {}

  CanonicalComponentFieldAdd(const CanonicalComponentFieldAdd&) = default;
  CanonicalComponentFieldAdd(CanonicalComponentFieldAdd&&) = default;

  // Assignment.
  CanonicalComponentFieldAdd& operator=(CanonicalComponentFieldAdd&) = default;
  CanonicalComponentFieldAdd& operator=(CanonicalComponentFieldAdd&&) = default;

 private:
  const CanonicalComponentFieldBase<N0, Derived0>& _u0;
  const CanonicalComponentFieldBase<N1, Derived1>& _u1;
};

// Subtract two fields.
template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
requires requires() {
  requires std::same_as<typename Derived0::Scalar, typename Derived1::Scalar>;
  requires(N0 == N1);
}
class CanonicalComponentFieldSubtract
    : public CanonicalComponentFieldBase<
          N0, CanonicalComponentFieldSubtract<N0, Derived0, N1, Derived1>> {
 public:
  using Int = Internal::Traits<
      CanonicalComponentFieldSubtract<N0, Derived0, N1, Derived1>>::Int;
  using Real = Internal::Traits<
      CanonicalComponentFieldSubtract<N0, Derived0, N1, Derived1>>::Real;
  using Complex = Internal::Traits<
      CanonicalComponentFieldSubtract<N0, Derived0, N1, Derived1>>::Complex;
  using Scalar = Internal::Traits<
      CanonicalComponentFieldSubtract<N0, Derived0, N1, Derived1>>::Scalar;
  using Value = Internal::Traits<
      CanonicalComponentFieldSubtract<N0, Derived0, N1, Derived1>>::Value;
  using Writeable = Internal::Traits<
      CanonicalComponentFieldSubtract<N0, Derived0, N1, Derived1>>::Writeable;

  auto& Grid() const { return _u0.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _u0[iTheta, iPhi] - _u1[iTheta, iPhi];
  }

  // Constructors.
  CanonicalComponentFieldSubtract() = delete;
  CanonicalComponentFieldSubtract(
      const CanonicalComponentFieldBase<N0, Derived0>& u0,
      const CanonicalComponentFieldBase<N1, Derived1>& u1)
      : _u0{u0}, _u1{u1} {}

  CanonicalComponentFieldSubtract(const CanonicalComponentFieldSubtract&) =
      default;
  CanonicalComponentFieldSubtract(CanonicalComponentFieldSubtract&&) = default;

  // Assignment.
  CanonicalComponentFieldSubtract& operator=(CanonicalComponentFieldSubtract&) =
      default;
  CanonicalComponentFieldSubtract& operator=(
      CanonicalComponentFieldSubtract&&) = default;

 private:
  const CanonicalComponentFieldBase<N0, Derived0>& _u0;
  const CanonicalComponentFieldBase<N1, Derived1>& _u1;
};

// Multiply two fields.
template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
requires requires() {
  requires std::same_as<typename Derived0::Scalar, typename Derived1::Scalar>;
}
class CanonicalComponentFieldMultiply
    : public CanonicalComponentFieldBase<
          N0 + N1,
          CanonicalComponentFieldMultiply<N0, Derived0, N1, Derived1>> {
 public:
  using Int = Internal::Traits<
      CanonicalComponentFieldMultiply<N0, Derived0, N1, Derived1>>::Int;
  using Real = Internal::Traits<
      CanonicalComponentFieldMultiply<N0, Derived0, N1, Derived1>>::Real;
  using Complex = Internal::Traits<
      CanonicalComponentFieldMultiply<N0, Derived0, N1, Derived1>>::Complex;
  using Scalar = Internal::Traits<
      CanonicalComponentFieldMultiply<N0, Derived0, N1, Derived1>>::Scalar;
  using Value = Internal::Traits<
      CanonicalComponentFieldMultiply<N0, Derived0, N1, Derived1>>::Value;
  using Writeable = Internal::Traits<
      CanonicalComponentFieldMultiply<N0, Derived0, N1, Derived1>>::Writeable;

  auto& Grid() const { return _u0.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _u0[iTheta, iPhi] * _u1[iTheta, iPhi];
  }

  // Constructors.
  CanonicalComponentFieldMultiply() = delete;
  CanonicalComponentFieldMultiply(
      const CanonicalComponentFieldBase<N0, Derived0>& u0,
      const CanonicalComponentFieldBase<N1, Derived1>& u1)
      : _u0{u0}, _u1{u1} {}

  CanonicalComponentFieldMultiply(const CanonicalComponentFieldMultiply&) =
      default;
  CanonicalComponentFieldMultiply(CanonicalComponentFieldMultiply&&) = default;

  // Assignment.
  CanonicalComponentFieldMultiply& operator=(CanonicalComponentFieldMultiply&) =
      default;
  CanonicalComponentFieldMultiply& operator=(
      CanonicalComponentFieldMultiply&&) = default;

 private:
  const CanonicalComponentFieldBase<N0, Derived0>& _u0;
  const CanonicalComponentFieldBase<N1, Derived1>& _u1;
};

template <std::ptrdiff_t N0, typename Derived0, std::ptrdiff_t N1,
          typename Derived1>
requires requires() {
  requires std::same_as<typename Derived0::Scalar, typename Derived1::Scalar>;
  requires(N1 == 0);
}
class CanonicalComponentFieldDivide
    : public CanonicalComponentFieldBase<
          N0, CanonicalComponentFieldDivide<N0, Derived0, N1, Derived1>> {
 public:
  using Int = Internal::Traits<
      CanonicalComponentFieldDivide<N0, Derived0, N1, Derived1>>::Int;
  using Real = Internal::Traits<
      CanonicalComponentFieldDivide<N0, Derived0, N1, Derived1>>::Real;
  using Complex = Internal::Traits<
      CanonicalComponentFieldDivide<N0, Derived0, N1, Derived1>>::Complex;
  using Scalar = Internal::Traits<
      CanonicalComponentFieldDivide<N0, Derived0, N1, Derived1>>::Scalar;
  using Value = Internal::Traits<
      CanonicalComponentFieldDivide<N0, Derived0, N1, Derived1>>::Value;
  using Writeable = Internal::Traits<
      CanonicalComponentFieldDivide<N0, Derived0, N1, Derived1>>::Writeable;

  auto& Grid() const { return _u0.Grid(); }

  auto operator[](Int iTheta, Int iPhi) const {
    this->CheckPointIndices(iTheta, iPhi);
    return _u0[iTheta, iPhi] / _u1[iTheta, iPhi];
  }

  // Constructors.
  CanonicalComponentFieldDivide() = delete;
  CanonicalComponentFieldDivide(
      const CanonicalComponentFieldBase<N0, Derived0>& u0,
      const CanonicalComponentFieldBase<N1, Derived1>& u1)
      : _u0{u0}, _u1{u1} {}

  CanonicalComponentFieldDivide(const CanonicalComponentFieldDivide&) = default;
  CanonicalComponentFieldDivide(CanonicalComponentFieldDivide&&) = default;

  // Assignment.
  CanonicalComponentFieldDivide& operator=(CanonicalComponentFieldDivide&) =
      default;
  CanonicalComponentFieldDivide& operator=(CanonicalComponentFieldDivide&&) =
      default;

 private:
  const CanonicalComponentFieldBase<N0, Derived0>& _u0;
  const CanonicalComponentFieldBase<N1, Derived1>& _u1;
};

}  // namespace GSHTrans

#endif