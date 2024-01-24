#ifndef GSH_TRANS_FIELDS_GUARD_H
#define GSH_TRANS_FIELDS_GUARD_H

#include <FFTWpp/All>
#include <algorithm>
#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>
#include <numeric>
#include <ranges>

#include "Concepts.h"
#include "FieldsBase.h"
#include "Indexing.h"

namespace GSHTrans {

/*

template <typename Grid, typename Range>
class CanonicalComponent
  : public CanonicalComponentBase<CanonicalComponent<Grid, Range>> {
using Int = std::ptrdiff_t;

public:
CanonicalComponent(Grid& grid, Range& range) : _grid{grid}, _range{range} {}

template <typename OtherDerived>
CanonicalComponent& operator=(CanonicalComponentBase<OtherDerived>&& other) {
  this->_grid = other.GridReference();
  std::ranges::copy(other, this->begin());
  return *this;
}

auto& operator[](Int i) { return _range[i]; }

auto& operator()(Int iTheta, Int iPhi) {
  auto i = this->NumberOfLongitudes() * iTheta + iPhi;
  return operator[](i);
}

private:
Grid& _grid;
Range& _range;

auto& _GridReference() const { return _grid; }
//  auto _View() { return boost::sub_range<Range>(_range); }
auto _View() { return std::ranges::views::all(_range); }

friend class CanonicalComponentBase<CanonicalComponent<Grid, Range>>;
};

// View to a canonical component produced by the action of a unary function.
template <typename View, typename Function>
class CanonicalComponentUnaryFunction
  : public CanonicalComponentBase<
        CanonicalComponentUnaryFunction<View, Function>> {
public:
CanonicalComponentUnaryFunction(View& f, Function&& function)
    : _f{f}, _function{function} {}

private:
View _f;
Function& _function;

auto& _GridReference() const { return _f.GridReference(); }
// auto _View() { return _f.View() | boost::adaptors::transformed(_function);
// }

auto _View() { return _f.View() | std::ranges::views::transform(_function); }

friend class CanonicalComponentBase<
    CanonicalComponentUnaryFunction<View, Function>>;
};

// Overload the unary minus operator.
template <typename Derived>
auto operator-(CanonicalComponentBase<Derived> f) {
return CanonicalComponentUnaryFunction(f, [](auto x) { return -x; });
}

/*
// Overload scalar multiplication.
template <typename Derived, typename Scalar>
requires std::integral<Scalar> or RealOrComplexFloatingPoint<Scalar>
auto operator*(CanonicalComponentBase<Derived> f, Scalar a) {
using F = decltype(f(0, 0));
return CanonicalComponentUnaryFunction(
  f, [a](auto x) { return static_cast<F>(a) * x; });
}

template <typename Derived, typename Scalar>
requires std::integral<Scalar> or RealOrComplexFloatingPoint<Scalar>
auto operator*(Scalar a, CanonicalComponentBase<Derived> f) { return f * a; }


// View to a canonical component produced by the action of a binary function.
template <typename Field1, typename Field2, typename Function>
class CanonicalComponentBinaryFunction
  : public CanonicalComponentBase<
        CanonicalComponentBinaryFunction<Field1, Field2, Function>> {
public:
CanonicalComponentBinaryFunction(Field1& field1, Field2& field2,
                                 Function&& function)
    : _field1{field1}, _field2{field2}, _function{function} {}

private:
Field1& _field1;
Field2& _field2;
Function& _function;

auto& _GridReference() const { return _field1.GridReference(); }

auto _View() {
  return boost::combine(_field1.View(), _field2.View()) |
         boost::adaptors::transformed([this](auto pair) {
           auto x = boost::get<0>(pair);
           auto y = boost::get<1>(pair);
           return _function(x, y);
         });
}

friend class CanonicalComponentBase<
    CanonicalComponentBinaryFunction<Field1, Field2, Function>>;
};

*/

  }  // namespace GSHTrans

#endif  // GSH_TRANS_FIELDS_GUARD_H
