#pragma once

#include <FFTWpp/Core>
#include <cassert>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "../Concepts.hpp"
#include "../Indexing.hpp"
#include "../Policies.hpp"
#include "../SpinField/SpinField.hpp"
#include "../SpinField/SpinFieldView.hpp"
#include "../SpinField/SpinWeighted.hpp"
#include "../Views.hpp"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                       A spin field in the spectral domain                 //
//--------------------------------------------------------------------------//

/// The coefficients f^N_{lm} of a field of upper index N, for degrees
/// |N| <= l <= lMax. See eq:expansion in docs/gshtrans-reference.tex.
///
/// This is the spectral counterpart of SpinField, and it exists for the same
/// reason: a raw coefficient buffer is untyped, so nothing stops a caller
/// reading the wrong block of one, and there is nowhere for the raising and
/// lowering operators to live. The layout is GSHIndices', which the transform
/// already writes, so this is a type over storage rather than a new storage
/// scheme.
///
/// Two things it records that a buffer cannot. The upper index is a
/// compile-time property, so the operators that change it change the type. And
/// the reduced m >= 0 storage of a real-valued field is a type distinction
/// rather than a convention: a RealValued expansion holds only the
/// non-negative orders, the rest being fixed by f_{l,-m} = (-1)^m conj(f_{lm}).
///
/// The degrees start at |N| because d^l_{mN} vanishes identically below it, so
/// there is no coefficient there to hold.
template <std::ptrdiff_t N_, AngularGrid Grid_,
          RealOrComplexValued Value_ = ComplexValued,
          typename Element_ = std::complex<typename Grid_::Real>>
class SpinExpansionBase {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.

  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = N_;
  using Value = Value_;    ///< Whether the samples are real-valued or complex.
  using GridType = Grid_;  ///< The angular grid this is defined on.
  using Real = typename Grid_::Real;   ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.

  /// Coefficients are complex whatever the field is; what a real field
  /// changes is how many of them there are.
  using Scalar = Complex;
  using MRange =
      std::conditional_t<std::same_as<Value, RealValued>, NonNegative,
                         All>;  ///< Which orders are stored.

  static_assert(std::same_as<Value, ComplexValued> or UpperIndex == 0,
                "A spin-weighted field can be real-valued only at upper index "
                "zero, so only there can its expansion use the reduced "
                "m >= 0 storage");

  SpinExpansionBase() = delete;

  /**
   * @brief A view of @p data as the coefficients of a field on @p grid.
   * @param grid The angular grid the coefficients belong to.
   * @param lMax The largest degree stored.
   * @param data The buffer, which must be the right size.
   * @throws std::invalid_argument if @p lMax is below the upper index, since
   * the harmonics there do not exist.
   */
  SpinExpansionBase(GridType grid, Int lMax, std::span<Element_> data)
      : grid_{std::move(grid)},
        indices_{Checked(lMax), lMax, UpperIndex},
        data_{data} {
    if (data_.size() != static_cast<std::size_t>(indices_.Size())) {
      throw std::invalid_argument(
          "An expansion over " + std::to_string(data_.size()) +
          " coefficients does not hold the " + std::to_string(indices_.Size()) +
          " this degree and upper index need");
    }
  }

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return grid_; }

  /** @brief The largest degree stored. */
  auto MaxDegree() const { return indices_.MaxDegree(); }
  /** @brief The smallest degree stored. */
  auto MinDegree() const { return indices_.MinDegree(); }
  /** @brief Every degree stored. */
  auto Degrees() const { return indices_.Degrees(); }
  /** @brief Every order stored at degree @p l. */
  auto Orders(Int l) const {
    return GSHSubIndices<MRange>(l, MaxDegree()).Orders();
  }
  /** @brief How many elements are stored. */
  auto Size() const { return static_cast<Int>(data_.size()); }
  /** @brief The underlying buffer. */
  auto Data() const { return data_; }

  /** @brief The coefficient at degree @p l and order @p m. */
  Complex operator[](Int l, Int m) const { return data_[Index(l, m)]; }

  /// The same, writable. Present only on a view over mutable storage.
  Element_& operator[](Int l, Int m)
  requires(not std::is_const_v<Element_>)
  {
    return data_[Index(l, m)];
  }

 private:
  // Checked on the way into the index block, which asserts the same thing: as
  // a member it is built before the constructor's body runs, so a check made
  // there would come too late in a Debug build, and the assert rather than
  // the documented exception would fire.
  static Int Checked(Int lMax) {
    if (lMax < std::abs(UpperIndex)) {
      throw std::invalid_argument(
          "An expansion's degree cannot be below its upper index, since the "
          "harmonics there do not exist");
    }
    return lMax;
  }

  GridType grid_;
  GSHIndices<MRange> indices_;
  std::span<Element_> data_;

  std::size_t Index(Int l, Int m) const {
    assert(l >= MinDegree() && l <= MaxDegree());
    return static_cast<std::size_t>(indices_.Index(l, m));
  }
};

// A non-owning expansion, which is what a tensor's component is.
template <std::ptrdiff_t N, AngularGrid Grid,
          RealOrComplexValued Value = ComplexValued>
using SpinExpansionView =
    SpinExpansionBase<N, Grid, Value, std::complex<typename Grid::Real>>;

template <std::ptrdiff_t N, AngularGrid Grid,
          RealOrComplexValued Value = ComplexValued>
using ConstSpinExpansionView =
    SpinExpansionBase<N, Grid, Value, const std::complex<typename Grid::Real>>;

/// The owning expansion: the same interface over storage it holds itself.
template <std::ptrdiff_t N_, AngularGrid Grid_,
          RealOrComplexValued Value_ = ComplexValued>
class SpinExpansion {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.

  /** @brief The upper index N of what this evaluates to. */
  static constexpr Int UpperIndex = N_;
  using Value = Value_;    ///< Whether the samples are real-valued or complex.
  using GridType = Grid_;  ///< The angular grid this is defined on.
  using Real = typename Grid_::Real;   ///< The precision.
  using Complex = std::complex<Real>;  ///< `std::complex` over the precision.
  /// A writable view over this object.
  using ViewType = SpinExpansionView<N_, Grid_, Value_>;
  /// A read-only view over this object.
  using ConstViewType = ConstSpinExpansionView<N_, Grid_, Value_>;
  using MRange =
      std::conditional_t<std::same_as<Value_, RealValued>, NonNegative,
                         All>;  ///< Which orders are stored.

  static_assert(std::same_as<Value_, ComplexValued> or UpperIndex == 0);

  SpinExpansion() = delete;

  /**
   * @brief A zero expansion on @p grid.
   * @param grid The angular grid the coefficients belong to.
   * @param lMax The largest degree stored.
   */
  SpinExpansion(GridType grid, Int lMax)
      : grid_{std::move(grid)},
        indices_{Checked(lMax), lMax, UpperIndex},
        data_(static_cast<std::size_t>(indices_.Size())) {}

  /** @brief The angular grid this is defined on. */
  const GridType& Grid() const { return grid_; }
  /** @brief The largest degree stored. */
  auto MaxDegree() const { return indices_.MaxDegree(); }
  /** @brief The smallest degree stored. */
  auto MinDegree() const { return indices_.MinDegree(); }
  /** @brief Every degree stored. */
  auto Degrees() const { return indices_.Degrees(); }
  /** @brief Every order stored at degree @p l. */
  auto Orders(Int l) const {
    return GSHSubIndices<MRange>(l, MaxDegree()).Orders();
  }
  /** @brief How many elements are stored. */
  auto Size() const { return static_cast<Int>(data_.size()); }
  /** @brief The underlying buffer. */
  auto Data() { return std::span<Complex>(data_); }
  /** @brief The underlying buffer. */
  auto Data() const { return std::span<const Complex>(data_); }

  /** @brief A writable view over this expansion's storage. */
  auto View() { return ViewType(grid_, MaxDegree(), Data()); }
  /** @brief A read-only view over this expansion's storage. */
  auto View() const { return ConstViewType(grid_, MaxDegree(), Data()); }

  /** @brief The coefficient at degree @p l and order @p m. */
  Complex operator[](Int l, Int m) const { return data_[Index(l, m)]; }
  /// The same, writable.
  Complex& operator[](Int l, Int m) { return data_[Index(l, m)]; }

 private:
  GridType grid_;
  GSHIndices<MRange> indices_;
  FFTWpp::vector<Complex> data_;

  static Int Checked(Int lMax) {
    if (lMax < std::abs(UpperIndex)) {
      throw std::invalid_argument(
          "An expansion's degree cannot be below its upper index, since the "
          "harmonics there do not exist");
    }
    return lMax;
  }

  std::size_t Index(Int l, Int m) const {
    assert(l >= MinDegree() && l <= MaxDegree());
    return static_cast<std::size_t>(indices_.Index(l, m));
  }
};

//--------------------------------------------------------------------------//
//                          Between the two domains                          //
//--------------------------------------------------------------------------//

/**
 * @brief Anything that holds the coefficients of one spin-weighted field.
 *
 * @details The owning SpinExpansion, and the views a tensor expansion hands
 * out for its components, which are the same thing over someone else's
 * storage. The free functions below and in Eth.hpp are written against this and
 * not against the owning type, so that a tensor's component can be raised,
 * lowered, evaluated or interpolated as it stands -- which is most of what one
 * wants to do with one.
 *
 * A layered expansion is deliberately not one: its coefficients are indexed
 * by radius as well, and it has free functions of its own.
 */
template <typename E>
concept SpinCoefficients =
    requires(const E& e, std::ptrdiff_t l, std::ptrdiff_t m) {
      requires std::same_as<std::remove_cv_t<decltype(E::UpperIndex)>,
                            std::ptrdiff_t>;
      typename E::Value;
      requires RealOrComplexValued<typename E::Value>;
      typename E::GridType;
      requires AngularGrid<typename E::GridType>;
      { e.Grid() } -> std::convertible_to<const typename E::GridType&>;
      { e.MinDegree() } -> std::convertible_to<std::ptrdiff_t>;
      { e.MaxDegree() } -> std::convertible_to<std::ptrdiff_t>;
      {
        e[l, m]
      } -> std::convertible_to<std::complex<typename E::GridType::Real>>;
      e.Data();
    };

// Transform a field into its expansion, and back.
//
// Free functions rather than members, because they belong to neither type: a
// transform is the grid's, and these only arrange the call. They are also the
// point at which a field's Value decides which of the transform's two paths
// runs -- a real field through the real one, with its reduced m >= 0 storage.
//
// Each comes in two forms. The returning form allocates its result. The
// writing form takes the result as an argument -- an owning object kept from
// call to call, or a view over the caller's own storage -- and allocates
// nothing; the degree is the result's. The returning form is the writing form
// applied to a new result, so the two give the same bits.

namespace SpinExpansionDetails {

// A field or view of upper index N and value kind V, whose storage can be
// written: what the writing Evaluate fills.
template <typename T, std::ptrdiff_t N, typename V>
concept WritableSpinField =
    SpinWeighted<T> && (T::UpperIndex == N) &&
    std::same_as<typename T::Value, V> && requires(T& t) {
      { t.Data() } -> std::same_as<std::span<typename T::Scalar>>;
    };

// Coefficients of upper index N and value kind V whose storage can be
// written: what the writing Expand fills.
template <typename T, std::ptrdiff_t N, typename V>
concept WritableSpinCoefficients =
    SpinCoefficients<T> && (T::UpperIndex == N) &&
    std::same_as<typename T::Value, V> && requires(T& t) {
      {
        t.Data()
      } -> std::same_as<std::span<std::complex<typename T::GridType::Real>>>;
    };

// The separation of successive samples in a field's storage: a view's stride,
// one for anything else.
template <typename T>
std::ptrdiff_t StrideOf(const T& field) {
  if constexpr (requires { field.Stride(); }) {
    return field.Stride();
  } else {
    return 1;
  }
}

// One field of `size` samples, `stride` apart, as a batch descriptor.
inline Batch SingleField(std::ptrdiff_t size, std::ptrdiff_t stride) {
  return stride == 1 ? Batch::One(size)
                     : Batch::Strided(1, stride, size * stride);
}

// Hand `run` a field's samples as the transform reads them, with the batch
// descriptor that says where they are. An owning field or a view is read in
// place, at its stride; an expression is evaluated into per-thread scratch,
// which is kept between calls so that a repeated transform allocates nothing.
template <typename F, typename Run>
void WithSamples(const F& field, Run&& run) {
  const auto size = field.Grid().FieldSize();
  if constexpr (requires { field.Data(); } &&
                (IsTerminal<F> || requires { field.Stride(); })) {
    run(field.Data(), SingleField(size, StrideOf(field)));
  } else {
    thread_local auto scratch = std::vector<typename F::Scalar>{};
    if (scratch.size() < static_cast<std::size_t>(size)) {
      scratch.resize(static_cast<std::size_t>(size));
    }
    auto samples = std::span(scratch.data(), static_cast<std::size_t>(size));
    field.EvaluateInto(samples);
    run(std::span<const typename F::Scalar>(samples), Batch::One(size));
  }
}

template <typename A, typename B>
void CheckSameGrid(const A& a, const B& b) {
  if (a.Grid().Identity() != b.Grid().Identity()) {
    throw std::invalid_argument(
        "A field and the expansion it is transformed to or from must be on "
        "the same grid");
  }
}

}  // namespace SpinExpansionDetails

/**
 * @brief Expands a field, writing into @p out.
 * @details The degree is @p out's. @p field may be an owning field, a view
 * (read in place at its stride) or a lazy expression (evaluated into
 * per-thread scratch first).
 * @param field What to expand.
 * @param out Where the coefficients go: a SpinExpansion, or a writable
 * SpinExpansionView over the caller's storage.
 * @param policy Whether the transform may thread.
 * @throws std::invalid_argument if the two are not on the same grid.
 */
template <typename FieldType, typename OutType>
requires SpinWeighted<std::remove_cvref_t<FieldType>> &&
         SpinExpansionDetails::WritableSpinCoefficients<
             std::remove_cvref_t<OutType>,
             std::remove_cvref_t<FieldType>::UpperIndex,
             typename std::remove_cvref_t<FieldType>::Value>
void Expand(const FieldType& field, OutType&& out,
            Execution policy = Execution::Sequential()) {
  using F = std::remove_cvref_t<FieldType>;
  SpinExpansionDetails::CheckSameGrid(field, out);
  auto target = out.Data();
  const auto lMax = static_cast<std::ptrdiff_t>(out.MaxDegree());
  SpinExpansionDetails::WithSamples(field, [&](auto samples, Batch batch) {
    field.Grid().ForwardTransformation(
        lMax, F::UpperIndex, samples, batch, target,
        Batch::One(static_cast<std::ptrdiff_t>(target.size())), policy);
  });
}

/**
 * @brief Expands a field to degree @p lMax.
 * @return A new SpinExpansion.
 */
template <typename FieldType>
requires SpinWeighted<std::remove_cvref_t<FieldType>>
auto Expand(const FieldType& field, std::ptrdiff_t lMax,
            Execution policy = Execution::Sequential()) {
  using F = std::remove_cvref_t<FieldType>;
  auto expansion =
      SpinExpansion<F::UpperIndex, typename F::GridType, typename F::Value>(
          field.Grid(), lMax);
  Expand(field, expansion, policy);
  return expansion;
}

/**
 * @brief Evaluates an expansion on its grid, writing into @p out.
 * @param expansion The coefficients, to their own degree.
 * @param out Where the samples go: a SpinField, or a writable SpinFieldView
 * over the caller's storage, at any stride.
 * @param policy Whether the transform may thread.
 * @throws std::invalid_argument if the two are not on the same grid.
 */
template <SpinCoefficients E, typename OutType>
requires SpinExpansionDetails::WritableSpinField<
    std::remove_cvref_t<OutType>, E::UpperIndex, typename E::Value>
void Evaluate(const E& expansion, OutType&& out,
              Execution policy = Execution::Sequential()) {
  SpinExpansionDetails::CheckSameGrid(expansion, out);
  auto target = out.Data();
  expansion.Grid().InverseTransformation(
      expansion.MaxDegree(), E::UpperIndex, expansion.Data(),
      Batch::One(static_cast<std::ptrdiff_t>(expansion.Data().size())), target,
      SpinExpansionDetails::SingleField(out.Grid().FieldSize(),
                                        SpinExpansionDetails::StrideOf(out)),
      policy);
}

/**
 * @brief Evaluates an expansion on its grid.
 * @return A new SpinField.
 */
template <SpinCoefficients E>
auto Evaluate(const E& expansion, Execution policy = Execution::Sequential()) {
  auto field =
      SpinField<E::UpperIndex, typename E::GridType, typename E::Value>(
          expansion.Grid());
  Evaluate(expansion, field, policy);
  return field;
}

}  // namespace GSHTrans
