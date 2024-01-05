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

template <typename Orders>
concept IndexRange =
    std::same_as<Orders, All> or std::same_as<Orders, NonNegative>;

// Tag classes and concepts for normalisations.
struct Ortho {};
struct FourPi {};

template <typename Norm>
concept Normalisation = std::same_as<Norm, Ortho> or std::same_as<Norm, FourPi>;

// Tag classes and concepts for transformation types.
struct C2C {
  using IndexRange = All;
};
struct R2C {
  using IndexRange = NonNegative;
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

}  // namespace GSHTrans

#endif  //  GSH_TRANS_CONCEPTS_GUARD_H
