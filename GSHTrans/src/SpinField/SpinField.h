#ifndef GSH_TRANS_SPIN_FIELD_GUARD_H
#define GSH_TRANS_SPIN_FIELD_GUARD_H

#include <FFTWpp/Core>
#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "../Concepts.h"
#include "SpinWeighted.h"

namespace GSHTrans {

// The owning terminal: a spin-weighted field of definite upper index N,
// holding its own samples on a grid.
//
// Named SpinField rather than CanonicalComponentField because "canonical
// component" belongs to the tensor layer, where a component is identified by a
// multi-index and not by an upper index -- for rank >= 2 the two are different
// things, and a collection labelled only by N does not determine a tensor
// (theory note section 2).
template <std::ptrdiff_t _N, AngularGrid _Grid,
          RealOrComplexValued _Value = ComplexValued>
class SpinField {
 public:
  using Int = std::ptrdiff_t;

  static constexpr Int UpperIndex = _N;
  using Value = _Value;
  using GridType = _Grid;
  using Real = typename _Grid::Real;
  using Complex = std::complex<Real>;
  using Scalar = ScalarFor<Real, Value>;

  // The reality constraint, restated here so that a mistake in a terminal's
  // template arguments reports itself rather than showing up as a failure to
  // satisfy SpinWeighted somewhere downstream.
  static_assert(std::same_as<Value, ComplexValued> or UpperIndex == 0,
                "A spin-weighted field can be real-valued only at upper index "
                "zero: real-valuedness is not preserved by the frame rotation "
                "e_{+-} -> e^{-+ i psi} e_{+-}");

  // Whether the grid admits negative upper indices at all is a property of its
  // type, so it is a static_assert and not a runtime test.
  static_assert(
      not std::same_as<typename _Grid::NRange, NonNegative> or UpperIndex >= 0,
      "This grid stores only non-negative upper indices");

  // A field is always built on a grid; there is no valid empty state.
  SpinField() = delete;

  explicit SpinField(GridType grid)
      : _grid{std::move(grid)},
        _data(CheckedSize(_grid), Scalar{}) {}

  // Sample a function of position. The callable is invoked as f(theta, phi)
  // over the grid's points, in the canonical order.
  template <typename Function>
  requires requires(Function f, Real theta, Real phi) {
    { f(theta, phi) } -> std::convertible_to<Scalar>;
  }
  SpinField(GridType grid, Function&& f) : SpinField(std::move(grid)) {
    auto iter = _data.begin();
    for (auto [theta, phi] : _grid.Points()) {
      *iter++ = static_cast<Scalar>(f(theta, phi));
    }
  }

  SpinField(const SpinField&) = default;
  SpinField(SpinField&&) = default;
  SpinField& operator=(const SpinField&) = default;
  SpinField& operator=(SpinField&&) = default;

  //------------------------------------------------------------------------//
  //                    Construction from an expression                     //
  //------------------------------------------------------------------------//

  // What may be evaluated into a field of this type: the same upper index, the
  // same precision, and a scalar that converts. The last condition is what
  // permits a real-valued expression into a complex field and makes the
  // reverse a compile error rather than a truncation.
  template <typename Expr>
  static constexpr bool Compatible =
      SpinWeighted<Node<Expr>> and Node<Expr>::UpperIndex == UpperIndex and
      std::same_as<typename Node<Expr>::Real, Real> and
      std::convertible_to<typename Node<Expr>::Scalar, Scalar>;

  // The grid comes from the expression, so a materialised field is on the
  // grid its operands were on and nowhere else.
  template <typename Expr>
  requires Compatible<Expr> and (not std::same_as<Node<Expr>, SpinField>)
  SpinField(const Expr& expr)
      : _grid{expr.Grid()}, _data(CheckedSize(_grid)) {
    expr.EvaluateInto(std::span<Scalar>(_data));
  }

  //------------------------------------------------------------------------//
  //                     Assignment from an expression                      //
  //------------------------------------------------------------------------//

  // Assignment never rebinds the destination grid: a field stays on the grid
  // it was built on, and an expression from elsewhere is an error rather than
  // a silent migration.
  //
  // Evaluation is in place and needs no temporary. Every node in this layer is
  // pointwise and index-preserving, so writing element (iTheta, iPhi) of the
  // destination happens after reading element (iTheta, iPhi) -- and only that
  // element -- of every operand, including the destination itself. `u = conj(u)
  // * v + u` is therefore safe as written. If a re-indexing node ever enters
  // this layer the argument fails, which is why the invariant is stated as
  // "pointwise *and index-preserving*" and why there is a regression test for
  // exactly this shape.
  template <typename Expr>
  requires Compatible<Expr> and (not std::same_as<Node<Expr>, SpinField>)
  SpinField& operator=(const Expr& expr) {
    CheckSameGrid(expr);
    expr.EvaluateInto(std::span<Scalar>(_data));
    return *this;
  }

  //------------------------------------------------------------------------//
  //                          Compound assignment                           //
  //------------------------------------------------------------------------//

  // Each carries exactly the constraints of "form the binary expression, then
  // assign it", so an unlawful use is an overload-resolution failure at the
  // call site rather than an error inside the body. That is what makes
  // static_assert(!requires { u *= v; }) a usable negative test.
  //
  // It also gets the index rules right without restating them: += needs equal
  // upper indices because operator+ does; *= and /= need the right-hand side
  // at zero because the result must still have the destination's index; and a
  // complex scalar cannot multiply a real field in place, because the
  // destination's value kind cannot change under assignment.
  template <typename Expr>
  static constexpr bool CanAddAssign =
      requires(const SpinField& f, const Expr& e) {
        requires Compatible<decltype(f + e)>;
      };

  template <typename Expr>
  static constexpr bool CanSubtractAssign =
      requires(const SpinField& f, const Expr& e) {
        requires Compatible<decltype(f - e)>;
      };

  template <typename Expr>
  static constexpr bool CanMultiplyAssign =
      requires(const SpinField& f, const Expr& e) {
        requires Compatible<decltype(f * e)>;
      };

  template <typename Expr>
  static constexpr bool CanDivideAssign =
      requires(const SpinField& f, const Expr& e) {
        requires Compatible<decltype(f / e)>;
      };

  // Expressions and scalars alike: whichever the right-hand side is, the
  // question is the same one.
  template <typename Expr>
  requires CanAddAssign<Expr>
  SpinField& operator+=(const Expr& expr) {
    return *this = *this + expr;
  }

  template <typename Expr>
  requires CanSubtractAssign<Expr>
  SpinField& operator-=(const Expr& expr) {
    return *this = *this - expr;
  }

  template <typename Expr>
  requires CanMultiplyAssign<Expr>
  SpinField& operator*=(const Expr& expr) {
    return *this = *this * expr;
  }

  template <typename Expr>
  requires CanDivideAssign<Expr>
  SpinField& operator/=(const Expr& expr) {
    return *this = *this / expr;
  }

  //------------------------------------------------------------------------//
  //                            The node interface                          //
  //------------------------------------------------------------------------//

  const GridType& Grid() const { return _grid; }

  // By value, as on every node: uniform value return is what lets terminals,
  // views and expressions be used interchangeably.
  Scalar operator[](Int iTheta, Int iPhi) const {
    return _data[FlatIndex(iTheta, iPhi)];
  }

  // Mutable access is a terminal's own, outside the concept.
  Scalar& operator[](Int iTheta, Int iPhi) {
    return _data[FlatIndex(iTheta, iPhi)];
  }

  // Terminals override the generic element loop with a contiguous copy.
  // Templated on the destination scalar so that a real-valued field can be
  // written into a complex destination; the constraint makes the reverse a
  // compile error rather than a silent truncation.
  template <typename S>
  requires std::convertible_to<Scalar, S>
  void EvaluateInto(std::span<S> target) const {
    CheckTargetSize(target.size());
    std::ranges::copy(_data, target.begin());
  }

  //------------------------------------------------------------------------//
  //                              Storage access                            //
  //------------------------------------------------------------------------//

  auto Size() const { return static_cast<Int>(_data.size()); }
  auto Data() const { return std::span<const Scalar>(_data); }
  auto Data() { return std::span<Scalar>(_data); }

  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }
  auto begin() const { return _data.begin(); }
  auto end() const { return _data.end(); }

 private:
  GridType _grid;
  FFTWpp::vector<Scalar> _data;

  // (iTheta, iPhi) with phi fastest, matching the transform's own layout.
  Int FlatIndex(Int iTheta, Int iPhi) const {
    assert(iTheta >= 0 && iTheta < static_cast<Int>(_grid.NumberOfCoLatitudes()));
    assert(iPhi >= 0 && iPhi < static_cast<Int>(_grid.NumberOfLongitudes()));
    return iTheta * static_cast<Int>(_grid.NumberOfLongitudes()) + iPhi;
  }

  // Checked in all build modes: whether this grid carries the field's upper
  // index is a configuration question, and a wrong answer is not something to
  // discover under NDEBUG.
  static std::size_t CheckedSize(const GridType& grid) {
    if (!std::ranges::contains(grid.UpperIndices(), UpperIndex)) {
      throw std::invalid_argument(
          "This grid does not carry upper index " +
          std::to_string(UpperIndex) + "; it carries " +
          std::to_string(*std::ranges::begin(grid.UpperIndices())) + " to " +
          std::to_string(grid.MaxUpperIndex()));
    }
    return static_cast<std::size_t>(grid.FieldSize());
  }

  template <typename Expr>
  void CheckSameGrid(const Expr& expr) const {
    if (expr.Grid().Identity() != _grid.Identity()) {
      throw std::invalid_argument(
          "Cannot assign an expression on a different grid: assignment does "
          "not rebind the destination's grid");
    }
  }

  void CheckTargetSize(std::size_t given) const {
    if (given != _data.size()) {
      throw std::invalid_argument(
          "Evaluation target has size " + std::to_string(given) +
          ", but this field has " + std::to_string(_data.size()) + " points");
    }
  }
};

template <std::ptrdiff_t N, AngularGrid Grid, RealOrComplexValued Value>
struct IsTerminalTrait<SpinField<N, Grid, Value>> : std::true_type {};

}  // namespace GSHTrans

#endif  // GSH_TRANS_SPIN_FIELD_GUARD_H
