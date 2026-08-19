#ifndef GSH_TRANS_SPIN_FIELD_OVERLOADS_GUARD_H
#define GSH_TRANS_SPIN_FIELD_OVERLOADS_GUARD_H

#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>

#include "../Concepts.h"
#include "SpinField.h"
#include "SpinFieldNodes.h"
#include "SpinWeighted.h"

namespace GSHTrans {

// Anything usable as an operand: a node, however it was passed. The parameter
// is deduced from a forwarding reference, so it may be a reference type.
template <typename T>
concept SpinFieldExpr = SpinWeighted<std::remove_cvref_t<T>>;

namespace SpinFieldOps {

// The functors. Each returns whatever the arithmetic gives, and the node
// derives its scalar and value kind from that -- so the reality constraint
// holds by construction rather than by a separate rule.

struct Plus {
  template <typename A, typename B>
  auto operator()(A a, B b) const {
    return a + b;
  }
};

struct Minus {
  template <typename A, typename B>
  auto operator()(A a, B b) const {
    return a - b;
  }
};

struct Times {
  template <typename A, typename B>
  auto operator()(A a, B b) const {
    return a * b;
  }
};

struct DividedBy {
  template <typename A, typename B>
  auto operator()(A a, B b) const {
    return a / b;
  }
};

struct Negated {
  template <typename A>
  auto operator()(A a) const {
    return -a;
  }
};

// Identity on a real scalar. std::conj would return a complex there and
// promote a real-valued field for no reason; the plan asks for conjugation at
// N = 0 on a real operand to be the identity, permitted and not special-cased.
struct Conjugated {
  template <typename A>
  auto operator()(A a) const {
    if constexpr (ComplexFloatingPoint<A>) {
      return std::conj(a);
    } else {
      return a;
    }
  }
};

struct Modulus {
  template <typename A>
  auto operator()(A a) const {
    return std::abs(a);
  }
};

// |f|^2 in its own right, not real(conj(f) * f), which would cost a complex
// multiply per point and a complex accumulator for a quantity known to be
// real.
struct SquaredModulus {
  template <typename A>
  auto operator()(A a) const {
    if constexpr (ComplexFloatingPoint<A>) {
      return std::norm(a);
    } else {
      return a * a;
    }
  }
};

struct RealPart {
  template <typename A>
  auto operator()(A a) const {
    return std::real(a);
  }
};

struct ImaginaryPart {
  template <typename A>
  auto operator()(A a) const {
    return std::imag(a);
  }
};

// Scalars are captured by value.
template <RealOrComplexFloatingPoint S>
struct TimesScalar {
  S s;
  template <typename A>
  auto operator()(A a) const {
    return s * a;
  }
};

template <RealOrComplexFloatingPoint S>
struct OverScalar {
  S s;
  template <typename A>
  auto operator()(A a) const {
    return a / s;
  }
};

template <RealOrComplexFloatingPoint S>
struct ScalarOver {
  S s;
  template <typename A>
  auto operator()(A a) const {
    return s / a;
  }
};

}  // namespace SpinFieldOps

// One precision per expression tree; no mixed-precision promotion.
template <typename L, typename R>
concept SamePrecisionAs =
    std::same_as<typename Node<L>::Real, typename Node<R>::Real>;

template <typename S, typename A>
concept ScalarFor_ = std::same_as<RemoveComplex<S>, typename Node<A>::Real>;

//--------------------------------------------------------------------------//
//                                  Binary                                   //
//--------------------------------------------------------------------------//

// Addition and subtraction: equal upper indices.
template <SpinFieldExpr L, SpinFieldExpr R>
requires SamePrecisionAs<L, R> and
    IndexRules::Equal::Admissible<Node<L>::UpperIndex, Node<R>::UpperIndex>
auto operator+(L&& l, R&& r) {
  return Binary<SpinFieldOps::Plus, IndexRules::Equal, L, R>(
      std::forward<L>(l), std::forward<R>(r));
}

template <SpinFieldExpr L, SpinFieldExpr R>
requires SamePrecisionAs<L, R> and
    IndexRules::Equal::Admissible<Node<L>::UpperIndex, Node<R>::UpperIndex>
auto operator-(L&& l, R&& r) {
  return Binary<SpinFieldOps::Minus, IndexRules::Equal, L, R>(
      std::forward<L>(l), std::forward<R>(r));
}

// Multiplication: upper indices add. The product of two band-limited fields
// exceeds the grid's truncation; the type system permits it, as it must, and
// the dealiasing question belongs to the grid -- see
// GaussLegendreGrid::ForBand.
template <SpinFieldExpr L, SpinFieldExpr R>
requires SamePrecisionAs<L, R> and
    IndexRules::Sum::Admissible<Node<L>::UpperIndex, Node<R>::UpperIndex>
auto operator*(L&& l, R&& r) {
  return Binary<SpinFieldOps::Times, IndexRules::Sum, L, R>(
      std::forward<L>(l), std::forward<R>(r));
}

// Division: by a scalar field only. Zeros of the divisor are the caller's
// problem and are not checked.
template <SpinFieldExpr L, SpinFieldExpr R>
requires SamePrecisionAs<L, R> and
    IndexRules::FirstOnly::Admissible<Node<L>::UpperIndex,
                                      Node<R>::UpperIndex>
auto operator/(L&& l, R&& r) {
  return Binary<SpinFieldOps::DividedBy, IndexRules::FirstOnly, L, R>(
      std::forward<L>(l), std::forward<R>(r));
}

//--------------------------------------------------------------------------//
//                                   Unary                                   //
//--------------------------------------------------------------------------//

template <SpinFieldExpr A>
auto operator-(A&& a) {
  return Unary<SpinFieldOps::Negated, IndexRules::Same, A>(std::forward<A>(a));
}

// Conjugation reverses the upper index (theory note section 5).
template <SpinFieldExpr A>
auto conj(A&& a) {
  return Unary<SpinFieldOps::Conjugated, IndexRules::Negate, A>(
      std::forward<A>(a));
}

// |f| is admissible at every upper index, and is real. Note that it is neither
// band-limited nor smooth at zeros of f, so the spectral layer must not assume
// an upper-index-zero expression is truncatable at the grid's lMax.
template <SpinFieldExpr A>
auto abs(A&& a) {
  return Unary<SpinFieldOps::Modulus, IndexRules::Zero, A>(std::forward<A>(a));
}

template <SpinFieldExpr A>
auto abs2(A&& a) {
  return Unary<SpinFieldOps::SquaredModulus, IndexRules::Zero, A>(
      std::forward<A>(a));
}

// real and imag are covariant only at upper index zero, so they exist only
// there. This was the third of the structural defects in the layer being
// replaced, which offered them at every upper index.
template <SpinFieldExpr A>
requires(Node<A>::UpperIndex == 0)
auto real(A&& a) {
  return Unary<SpinFieldOps::RealPart, IndexRules::Zero, A>(
      std::forward<A>(a));
}

template <SpinFieldExpr A>
requires(Node<A>::UpperIndex == 0)
auto imag(A&& a) {
  return Unary<SpinFieldOps::ImaginaryPart, IndexRules::Zero, A>(
      std::forward<A>(a));
}

//--------------------------------------------------------------------------//
//                             Scalar arithmetic                             //
//--------------------------------------------------------------------------//

// A real scalar preserves the value kind; a complex one promotes it. Neither
// touches the upper index.
template <SpinFieldExpr A, RealOrComplexFloatingPoint S>
requires ScalarFor_<S, A>
auto operator*(A&& a, S s) {
  return Unary<SpinFieldOps::TimesScalar<S>, IndexRules::Same, A>(
      std::forward<A>(a), SpinFieldOps::TimesScalar<S>{s});
}

template <RealOrComplexFloatingPoint S, SpinFieldExpr A>
requires ScalarFor_<S, A>
auto operator*(S s, A&& a) {
  return Unary<SpinFieldOps::TimesScalar<S>, IndexRules::Same, A>(
      std::forward<A>(a), SpinFieldOps::TimesScalar<S>{s});
}

template <SpinFieldExpr A, RealOrComplexFloatingPoint S>
requires ScalarFor_<S, A>
auto operator/(A&& a, S s) {
  return Unary<SpinFieldOps::OverScalar<S>, IndexRules::Same, A>(
      std::forward<A>(a), SpinFieldOps::OverScalar<S>{s});
}

// Dividing a scalar by a field needs the field at upper index zero, for the
// same reason division between fields does.
template <RealOrComplexFloatingPoint S, SpinFieldExpr A>
requires ScalarFor_<S, A> and (Node<A>::UpperIndex == 0)
auto operator/(S s, A&& a) {
  return Unary<SpinFieldOps::ScalarOver<S>, IndexRules::Same, A>(
      std::forward<A>(a), SpinFieldOps::ScalarOver<S>{s});
}

//--------------------------------------------------------------------------//
//                              Callable nodes                               //
//--------------------------------------------------------------------------//

// Apply an arbitrary callable pointwise. This is where pow, exp, log and
// anything else of that kind live: they are all Map at upper index zero, and
// naming each of them would be sugar over one node.
//
// Called Map rather than Transform, which would collide with
// ForwardTransformation and InverseTransformation and with "spherical harmonic
// transform" throughout this codebase's vocabulary.
//
// Restricted to upper index zero, like real and imag: a function applied to a
// component's value is a statement about that value, and only at N = 0 is the
// value frame-independent enough for the statement to mean anything.
//
// The callable is decayed and stored by value, so the node owns it and an
// expression outlives the caller's lambda. Its result determines the node's
// scalar and hence its value kind -- a trait on the return type rather than a
// guess.
//
// Invoked as f(value). Point-dependent callables, f(theta, phi, value), are
// deliberately not offered here; if wanted they are a second overload rather
// than a change to this one.
template <SpinFieldExpr A, typename F>
requires(Node<A>::UpperIndex == 0) and
    std::invocable<std::decay_t<F>, typename Node<A>::Scalar>
auto Map(A&& a, F&& f) {
  using Functor = std::decay_t<F>;
  return Unary<Functor, IndexRules::Zero, A>(std::forward<A>(a),
                                             Functor(std::forward<F>(f)));
}

//--------------------------------------------------------------------------//
//                             Materialisation                               //
//--------------------------------------------------------------------------//

// Break a lazy chain deliberately -- before feeding a product into a transform
// twice, say, or before a loop that would otherwise re-evaluate it.
template <SpinFieldExpr Expr>
auto Materialise(Expr&& expr) {
  using E = Node<Expr>;
  return SpinField<E::UpperIndex, typename E::GridType, typename E::Value>(
      expr);
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_SPIN_FIELD_OVERLOADS_GUARD_H
