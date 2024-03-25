#ifndef GSH_TRANS_CONCEPTS_GUARD_H
#define GSH_TRANS_CONCEPTS_GUARD_H

#include <complex>
#include <concepts>
#include <iterator>
#include <type_traits>

namespace GSHTrans {

//-------------------------------------------------------------------------//
//                                 Tag classes                              //
//--------------------------------------------------------------------------//

// Matrix storage options.
struct ColumnMajor {};
struct RowMajor {};

template <typename T>
concept MatrixStorage =
    std::same_as<T, ColumnMajor> or std::same_as<T, RowMajor>;

// Index storage options.
struct All {};
struct NonNegative {};
struct Single {};
struct Multiple {};

struct UpperIndexFirst {};
struct AngleFirst {};

template <typename Indices>
concept IndexRange =
    std::same_as<Indices, All> or std::same_as<Indices, NonNegative> or
    std::same_as<Indices, Single>;

template <typename Indices>
concept OrderIndexRange =
    std::same_as<Indices, All> or std::same_as<Indices, NonNegative>;

template <typename Indices>
concept AngleIndexRange =
    std::same_as<Indices, Multiple> or std::same_as<Indices, Single>;

// Normalisation options.
struct Ortho {};
struct FourPi {};

template <typename Norm>
concept Normalisation = std::same_as<Norm, Ortho> or std::same_as<Norm, FourPi>;

// Value type options.
struct RealValued {};
struct ComplexValued {};

template <typename T>
concept RealOrComplexValued =
    std::same_as<T, RealValued> or std::same_as<T, ComplexValued>;

// Concepts for fields.
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
struct RemoveComplexHelper {
  using value_type = T;
};

template <typename T>
struct RemoveComplexHelper<std::complex<T>> {
  using value_type = T;
};

template <typename T>
using RemoveComplex = typename RemoveComplexHelper<T>::value_type;

}  // namespace GSHTrans

#endif  //  GSH_TRANS_CONCEPTS_GUARD_H
