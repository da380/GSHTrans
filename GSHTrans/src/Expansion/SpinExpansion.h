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
#include "../Policies.h"
#include "../Indexing.h"
#include "../SpinField/SpinField.h"
#include "../SpinField/SpinFieldView.h"
#include "../SpinField/SpinWeighted.h"
#include "../Views.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                       A spin field in the spectral domain                 //
//--------------------------------------------------------------------------//

// The coefficients f^N_{lm} of a field of upper index N (theory note
// eq:expansion), for degrees |N| <= l <= lMax.
//
// This is the spectral counterpart of SpinField, and it exists for the same
// reason: a raw coefficient buffer is untyped, so nothing stops a caller
// reading the wrong block of one, and there is nowhere for the raising and
// lowering operators to live. The layout is GSHIndices', which the transform
// already writes, so this is a type over storage rather than a new storage
// scheme.
//
// Two things it records that a buffer cannot. The upper index is a
// compile-time property, so the operators that change it change the type. And
// the reduced m >= 0 storage of a real-valued field is a type distinction
// rather than a convention: a RealValued expansion holds only the
// non-negative orders, the rest being fixed by f_{l,-m} = (-1)^m conj(f_{lm}).
//
// The degrees start at |N| because d^l_{mN} vanishes identically below it, so
// there is no coefficient there to hold.
template <std::ptrdiff_t _N, AngularGrid _Grid,
          RealOrComplexValued _Value = ComplexValued,
          typename _Element = std::complex<typename _Grid::Real>>
class SpinExpansionBase {
 public:
  using Int = std::ptrdiff_t;

  static constexpr Int UpperIndex = _N;
  using Value = _Value;
  using GridType = _Grid;
  using Real = typename _Grid::Real;
  using Complex = std::complex<Real>;

  // Coefficients are complex whatever the field is; what a real field changes
  // is how many of them there are.
  using Scalar = Complex;
  using MRange =
      std::conditional_t<std::same_as<Value, RealValued>, NonNegative, All>;

  static_assert(std::same_as<Value, ComplexValued> or UpperIndex == 0,
                "A spin-weighted field can be real-valued only at upper index "
                "zero, so only there can its expansion use the reduced "
                "m >= 0 storage");

  SpinExpansionBase() = delete;

  SpinExpansionBase(GridType grid, Int lMax, std::span<_Element> data)
      : _grid{std::move(grid)},
        _indices{lMax, lMax, UpperIndex},
        _data{data} {
    if (lMax < std::abs(UpperIndex)) {
      throw std::invalid_argument(
          "An expansion's degree cannot be below its upper index, since the "
          "harmonics there do not exist");
    }
    if (_data.size() != static_cast<std::size_t>(_indices.Size())) {
      throw std::invalid_argument(
          "An expansion over " + std::to_string(_data.size()) +
          " coefficients does not hold the " + std::to_string(_indices.Size()) +
          " this degree and upper index need");
    }
  }

  const GridType& Grid() const { return _grid; }

  auto MaxDegree() const { return _indices.MaxDegree(); }
  auto MinDegree() const { return _indices.MinDegree(); }
  auto Degrees() const { return _indices.Degrees(); }
  auto Orders(Int l) const { return GSHSubIndices<MRange>(l, MaxDegree()).Orders(); }
  auto Size() const { return static_cast<Int>(_data.size()); }
  auto Data() const { return _data; }

  Complex operator[](Int l, Int m) const { return _data[Index(l, m)]; }

  _Element& operator[](Int l, Int m)
  requires(not std::is_const_v<_Element>)
  {
    return _data[Index(l, m)];
  }

 private:
  GridType _grid;
  GSHIndices<MRange> _indices;
  std::span<_Element> _data;

  std::size_t Index(Int l, Int m) const {
    assert(l >= MinDegree() && l <= MaxDegree());
    return static_cast<std::size_t>(_indices.Index(l, m));
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

// The owning expansion: the same interface over storage it holds itself.
template <std::ptrdiff_t _N, AngularGrid _Grid,
          RealOrComplexValued _Value = ComplexValued>
class SpinExpansion {
 public:
  using Int = std::ptrdiff_t;

  static constexpr Int UpperIndex = _N;
  using Value = _Value;
  using GridType = _Grid;
  using Real = typename _Grid::Real;
  using Complex = std::complex<Real>;
  using ViewType = SpinExpansionView<_N, _Grid, _Value>;
  using ConstViewType = ConstSpinExpansionView<_N, _Grid, _Value>;
  using MRange =
      std::conditional_t<std::same_as<_Value, RealValued>, NonNegative, All>;

  static_assert(std::same_as<_Value, ComplexValued> or UpperIndex == 0);

  SpinExpansion() = delete;

  SpinExpansion(GridType grid, Int lMax)
      : _grid{std::move(grid)},
        _indices{Checked(lMax), lMax, UpperIndex},
        _data(static_cast<std::size_t>(_indices.Size())) {}

  const GridType& Grid() const { return _grid; }
  auto MaxDegree() const { return _indices.MaxDegree(); }
  auto MinDegree() const { return _indices.MinDegree(); }
  auto Degrees() const { return _indices.Degrees(); }
  auto Orders(Int l) const {
    return GSHSubIndices<MRange>(l, MaxDegree()).Orders();
  }
  auto Size() const { return static_cast<Int>(_data.size()); }
  auto Data() { return std::span<Complex>(_data); }
  auto Data() const { return std::span<const Complex>(_data); }

  auto View() { return ViewType(_grid, MaxDegree(), Data()); }
  auto View() const { return ConstViewType(_grid, MaxDegree(), Data()); }

  Complex operator[](Int l, Int m) const { return _data[Index(l, m)]; }
  Complex& operator[](Int l, Int m) { return _data[Index(l, m)]; }

 private:
  GridType _grid;
  GSHIndices<MRange> _indices;
  FFTWpp::vector<Complex> _data;

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
    return static_cast<std::size_t>(_indices.Index(l, m));
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
