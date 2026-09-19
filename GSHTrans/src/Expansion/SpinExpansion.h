#ifndef GSH_TRANS_SPIN_EXPANSION_GUARD_H
#define GSH_TRANS_SPIN_EXPANSION_GUARD_H

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

#include "../Concepts.h"
#include "../Indexing.h"
#include "../Policies.h"
#include "../SpinField/SpinField.h"
#include "../SpinField/SpinFieldView.h"
#include "../SpinField/SpinWeighted.h"
#include "../Views.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                       A spin field in the spectral domain                 //
//--------------------------------------------------------------------------//

/// The coefficients f^N_{lm} of a field of upper index N, for degrees
/// |N| <= l <= lMax. See eq:expansion of the theory note,
/// docs/canonical-components.tex.
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
      : grid_{std::move(grid)}, indices_{lMax, lMax, UpperIndex}, data_{data} {
    if (lMax < std::abs(UpperIndex)) {
      throw std::invalid_argument(
          "An expansion's degree cannot be below its upper index, since the "
          "harmonics there do not exist");
    }
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

// Transform a field into its expansion, and back.
//
// Free functions rather than members, because they belong to neither type: a
// transform is the grid's, and these only arrange the call. They are also the
// point at which a field's Value decides which of the transform's two paths
// runs -- a real field through the real one, with its reduced m >= 0 storage.
template <typename FieldType>
requires SpinWeighted<std::remove_cvref_t<FieldType>>
auto Expand(const FieldType& field, std::ptrdiff_t lMax,
            Execution policy = Execution::Sequential()) {
  using F = std::remove_cvref_t<FieldType>;
  auto expansion =
      SpinExpansion<F::UpperIndex, typename F::GridType, typename F::Value>(
          field.Grid(), lMax);
  auto values = std::vector<typename F::Scalar>(
      static_cast<std::size_t>(field.Grid().FieldSize()));
  field.EvaluateInto(std::span(values));
  auto out = expansion.Data();
  field.Grid().ForwardTransformation(lMax, F::UpperIndex, values, out, policy);
  return expansion;
}

template <std::ptrdiff_t N, AngularGrid Grid, RealOrComplexValued Value>
auto Evaluate(const SpinExpansion<N, Grid, Value>& expansion,
              Execution policy = Execution::Sequential()) {
  auto field = SpinField<N, Grid, Value>(expansion.Grid());
  auto out = field.Data();
  expansion.Grid().InverseTransformation(expansion.MaxDegree(), N,
                                         expansion.Data(), out, policy);
  return field;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_SPIN_EXPANSION_GUARD_H
