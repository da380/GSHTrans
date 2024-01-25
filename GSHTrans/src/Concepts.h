#ifndef GSH_TRANS_CONCEPTS_GUARD_H
#define GSH_TRANS_CONCEPTS_GUARD_H

#include <complex>
#include <concepts>
#include <iterator>
#include <type_traits>

namespace GSHTrans {

// Tag classes and concepts for order storage.
struct All {};
struct NonNegative {};
struct Single {};

template <typename Indices>
concept IndexRange = std::same_as<Indices, All> or
    std::same_as<Indices, NonNegative> or std::same_as<Indices, Single>;

template <typename Indices>
concept OrderIndexRange =
    std::same_as<Indices, All> or std::same_as<Indices, NonNegative>;

template <typename Indices>
concept AngleIndexRange =
    std::same_as<Indices, NonNegative> or std::same_as<Indices, Single>;

// Tag classes and concepts for normalisations.
struct Ortho {};
struct FourPi {};

template <typename Norm>
concept Normalisation = std::same_as<Norm, Ortho> or std::same_as<Norm, FourPi>;

// Tag classes and concepts for transformation types.
struct C2C {
  using IndexRange = All;

  template <std::floating_point Real>
  using Scalar = std::complex<Real>;
};
struct R2C {
  using IndexRange = NonNegative;

  template <std::floating_point Real>
  using Scalar = Real;
};

template <typename Type>
concept TransformType = std::same_as<Type, C2C> or std::same_as<Type, R2C>;

// Concepts for real or complex floating point types.
template <typename T>
struct IsComplexFloatingPoint : public std::false_type {};

template <typename T>
struct IsComplexFloatingPoint<std::complex<T>>
    : public std::bool_constant<std::is_floating_point_v<T>> {};

template <typename T>
concept RealFloatingPoint = std::floating_point<T>;

template <typename T>
concept ComplexFloatingPoint =
    IsComplexFloatingPoint<std::remove_const_t<T>>::value;

template <typename T>
concept RealOrComplexFloatingPoint =
    RealFloatingPoint<T> or ComplexFloatingPoint<T>;

template <typename T>
concept IntegerOrRealOrComplexFloatingPoint =
    std::integral<T> or RealOrComplexFloatingPoint<T>;

template <typename T>
struct RemoveComplexHelper {
  using value_type = T;
};

template <typename T>
struct RemoveComplexHelper<std::complex<T>> {
  using value_type = T;
};

template <typename T>
using RemoveComplex = typename RemoveComplexHelper<T>::value_type;

// Concepts for iterators with real or complex floating point values.
template <typename T>
concept RealFloatingPointIterator = requires() {
  requires std::random_access_iterator<T>;
  requires RealFloatingPoint<std::iter_value_t<T>>;
};

template <typename T>
concept ComplexFloatingPointIterator = requires() {
  requires std::random_access_iterator<T>;
  requires ComplexFloatingPoint<std::iter_value_t<T>>;
};

template <typename T>
concept RealOrComplexFloatingPointIterator = requires() {
  requires std::random_access_iterator<T>;
  requires RealOrComplexFloatingPoint<std::iter_value_t<T>>;
};

// Concepts for ranges with real or complex floating point values.
template <typename T>
concept RealFloatingPointRange = requires() {
  requires std::ranges::random_access_range<T>;
  requires RealFloatingPoint<std::ranges::range_value_t<T>>;
};

template <typename T>
concept ComplexFloatingPointRange = requires() {
  requires std::ranges::random_access_range<T>;
  requires ComplexFloatingPoint<std::ranges::range_value_t<T>>;
};

template <typename T>
concept RealOrComplexFloatingPointRange = requires() {
  requires std::ranges::random_access_range<T>;
  requires RealOrComplexFloatingPoint<std::ranges::range_value_t<T>>;
};

// Concepts for scalar-valued functions.
template <typename Func, typename Real>
concept ScalarFunction2D = requires(Real theta, Real phi, Real w, Func f) {
  RealFloatingPoint<Real>;
  std::invocable<Func, Real, Real>;
  { f(theta, phi) } -> RealOrComplexFloatingPoint<>;
  { f(theta, phi) * w } -> RealOrComplexFloatingPoint<>;
};

// Concepts for grid classes.
template <typename Grid>
concept SphereGrid = requires(Grid grid) {
  typename Grid::difference_type;
  { grid.MaxDegree() } -> std::same_as<typename Grid::difference_type>;
  { grid.MaxUpperIndex() } -> std::same_as<typename Grid::difference_type>;
};

}  // namespace GSHTrans

#endif  //  GSH_TRANS_CONCEPTS_GUARD_H
