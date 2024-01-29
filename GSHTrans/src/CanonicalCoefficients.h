#ifndef GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H
#define GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H

#include <FFTWpp/All>
#include <algorithm>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/tuple/tuple.hpp>
#include <complex>
#include <functional>
#include <iostream>
#include <ranges>
#include <type_traits>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {

//----------------------------------------------------//
//       Base class for canonical coefficients        //
//----------------------------------------------------//

template <typename Derived>
class CanonicalCoefficientBase {
  using Int = std::ptrdiff_t;

 public:
  auto DataView() { return GetDerived()._DataView(); }
  auto& Grid() const { return GetDerived()._Grid(); }

  auto begin() { return DataView().begin(); }
  auto end() { return DataView().end(); }
  auto size() { return DataView().size(); }

  auto MinDegree() { return _GSHView().MinDegree(); }
  auto MaxDegree() { return _GSHView().MaxDegree(); }
  auto UpperIndex() { return _GSHView().UpperIndex(); }
  auto Degrees() { return _GSHView().Degrees(); }

  auto& operator[](Int i) const {
    assert(i >= 0 && i < size());
    return DataView()[i];
  }

  auto& operator[](Int i) {
    assert(i >= 0 && i < size());
    return DataView()[i];
  }

  auto operator()(Int l) { return _GSHView()(l); }

 private:
  auto& GetDerived() { return static_cast<Derived&>(*this); }
  auto& GetDerived() const { return static_cast<const Derived&>(*this); }

  auto _GSHView() {
    auto lMax = GetDerived()._MaxDegree();
    auto n = GetDerived()._UpperIndex();
    return GSHTrans::GSHView<typename Derived::value_type,
                             typename Derived::MRange>(lMax, lMax, n, begin());
  }
};

//----------------------------------------------------//
//             Canonical coefficients class           //
//----------------------------------------------------//

template <typename GSHGrid, RealOrComplexValued ValueType>
class CanonicalCoefficient : public CanonicalCoefficientBase<
                                 CanonicalCoefficient<GSHGrid, ValueType>> {
  using Int = std::ptrdiff_t;
  using Real = typename GSHGrid::real_type;
  using Complex = typename GSHGrid::complex_type;
  using Vector = FFTWpp::vector<Complex>;

 public:
  using value_type = Complex;
  using CoefficientType = ValueType;
  using MRange =
      std::conditional_t<std::same_as<ValueType, RealValued>, NonNegative, All>;

  CanonicalCoefficient(Int lMax, Int n, GSHGrid& grid)
      : _lMax{lMax},
        _n{n},
        _grid{grid},
        _data{Vector(GSHIndices<MRange>(_lMax, _lMax, n).size())} {
    assert(_lMax >= std::abs(_n) && _lMax <= _grid.MaxDegree());
    assert(std::ranges::any_of(_grid.UpperIndices(),
                               [n](auto np) { return n == np; }));
  }

  CanonicalCoefficient(Int n, GSHGrid& grid)
      : CanonicalCoefficient(grid.MaxDegree(), n, grid) {}

  CanonicalCoefficient(CanonicalCoefficient&) = default;
  CanonicalCoefficient(CanonicalCoefficient&&) = default;

  template <typename Derived>
  CanonicalCoefficient(CanonicalCoefficientBase<Derived>& other)
      : _lMax{other.MaxDegree()},
        _n{other.UpperIndex()},
        _grid{other.grid()},
        _data{Vector(other.begin(), other.end())} {}

  template <typename Derived>
  CanonicalCoefficient(CanonicalCoefficientBase<Derived>&& other)
      : CanonicalCoefficient(other) {}

  template <typename Derived>
  requires std::same_as<typename Derived::CoefficientType, CoefficientType> and
      std::same_as<typename Derived::value_type, value_type>
  auto& operator=(CanonicalCoefficientBase<Derived>& other) {
    assert(_lMax == other.MaxDegree());
    assert(std::abs(_n) == other.UpperIndex());
    _n = other.UpperIndex();
    std::ranges::copy(other.DataView(), _data.begin());
    return *this;
  }

  template <typename Derived>
  requires std::same_as<typename Derived::CoefficientType, CoefficientType> and
      std::same_as<typename Derived::value_type, value_type>
  auto& operator=(CanonicalCoefficientBase<Derived>&& other) {
    *this = other;
    return *this;
  }

 private:
  Int _lMax;
  Int _n;
  GSHGrid& _grid;
  Vector _data;

  auto _DataView() { return boost::sub_range<Vector>(_data); }
  auto& _Grid() const { return _grid; }
  auto _MaxDegree() { return _lMax; }
  auto _UpperIndex() { return _n; }

  friend class CanonicalCoefficientBase<
      CanonicalCoefficient<GSHGrid, ValueType>>;
};

//-----------------------------------------------------------------//
//    Canonical coeffcient class that stores a view to its data    //
//-----------------------------------------------------------------//

template <typename GSHGrid, RealOrComplexValued ValueType,
          std::ranges::common_range View>
requires std::same_as < RemoveComplex<std::ranges::range_value_t<View>>,
typename GSHGrid::real_type > class CanonicalCoefficientView
    : public CanonicalCoefficientBase<
          CanonicalCoefficientView<GSHGrid, ValueType, View>> {
  using Int = std::ptrdiff_t;

 public:
  using value_type = std::ranges::range_value_t<View>;
  using CoefficientType = ValueType;
  using MRange =
      std::conditional_t<std::same_as<ValueType, RealValued>, NonNegative, All>;

  CanonicalCoefficientView(Int lMax, Int n, GSHGrid& grid, View view, ValueType)
      : _lMax{lMax}, _n{n}, _grid{grid}, _view{view} {
    assert(_lMax >= std::abs(_n) && _lMax <= _grid.MaxDegree());
    assert(std::ranges::any_of(_grid.UpperIndices(),
                               [n](auto np) { return n == np; }));
  }

  template <typename Derived>
  requires std::same_as<typename Derived::CoefficientType, CoefficientType> and
      std::same_as<typename Derived::value_type, value_type>
  auto& operator=(CanonicalCoefficientBase<Derived>& other) {
    assert(_lMax == other.MaxDegree());
    assert(std::abs(_n) == other.UpperIndex());
    _n = other.UpperIndex();
    std::ranges::copy(other.DataView(), _view.begin());
    return *this;
  }

 private:
  Int _lMax;
  Int _n;
  GSHGrid& _grid;
  View _view;

  auto _DataView() { return _view; }
  auto& _Grid() const { return _grid; }
  auto _MaxDegree() { return _lMax; }
  auto _UpperIndex() { return _n; }

  friend class CanonicalCoefficientBase<
      CanonicalCoefficientView<GSHGrid, ValueType, View>>;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H
