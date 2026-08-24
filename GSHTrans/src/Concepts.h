#ifndef GSH_TRANS_CONCEPTS_GUARD_H
#define GSH_TRANS_CONCEPTS_GUARD_H

/**
 * @file Concepts.h
 * @brief The tag types that select the library's storage options, and the
 * numeric concepts everything else is written against.
 *
 * @details The tags are empty types used as template arguments. They are tags
 * rather than enumerators because the choices they express are made at compile
 * time and change the type of what is built: an expansion storing only
 * non-negative orders is a different type from one storing all of them, and
 * keeping that distinction in the type system is the point.
 */

#include <NumericConcepts/Numeric.hpp>
#include <concepts>
#include <cstddef>

namespace GSHTrans {

//-------------------------------------------------------------------------//
//                    Which orders and upper indices exist                  //
//--------------------------------------------------------------------------//

/**
 * @brief Every value in range.
 * @details Orders @f$-l \le m \le l@f$, or upper indices
 * @f$-n_{\max} \le n \le n_{\max}@f$.
 */
struct All {};

/**
 * @brief Non-negative values only.
 * @details For orders this is the reduced storage of a real-valued field,
 * whose negative orders follow from
 * @f$f_{l,-m} = (-1)^m \overline{f_{lm}}@f$. For upper indices it is a grid
 * that serves real scalars alone.
 */
struct NonNegative {};

/** @brief One value only, the largest in range. */
struct Single {};

/** @brief More than one value. */
struct Multiple {};

/** @brief Which orders a coefficient block stores. */
template <typename Indices>
concept OrderIndexRange =
    std::same_as<Indices, All> or std::same_as<Indices, NonNegative>;

/** @brief Which upper indices a grid or a table of Wigner values covers. */
template <typename Indices>
concept IndexRange =
    std::same_as<Indices, All> or std::same_as<Indices, NonNegative> or
    std::same_as<Indices, Single>;

/** @brief Whether a table of Wigner values holds one colatitude or many. */
template <typename Indices>
concept AngleIndexRange =
    std::same_as<Indices, Multiple> or std::same_as<Indices, Single>;

//-------------------------------------------------------------------------//
//                          Real-valued or complex                          //
//--------------------------------------------------------------------------//

/**
 * @brief The field's samples are real.
 * @details Available only at upper index zero. Real-valuedness is not
 * preserved by the frame rotation @f$e_\pm \mapsto e^{\mp i\psi}e_\pm@f$, so
 * it is not a property any component of any tensor can have at
 * @f$N \ne 0@f$.
 */
struct RealValued {};

/** @brief The field's samples are complex. */
struct ComplexValued {};

/** @brief Either of the two value kinds. */
template <typename T>
concept RealOrComplexValued =
    std::same_as<T, RealValued> or std::same_as<T, ComplexValued>;

//-------------------------------------------------------------------------//
//                      Numeric concepts, from elsewhere                    //
//--------------------------------------------------------------------------//

// These are NumericConcepts' definitions rather than this library's. The point
// of that dependency is that the projects built on it agree on what "real",
// "complex" and "precision" mean, and a second, subtly different set of
// definitions here would defeat it.
//
// They are re-spelled rather than imported unqualified because NumericConcepts
// names them Real and Complex, while this library uses Real and Complex as
// member type aliases in nearly every class. Importing those names would have
// them shadowed at exactly the points where the concept is most likely to be
// wanted, and the resulting errors would be obscure. Anything new should use
// the NumericConcepts spelling directly, qualified.

/** @brief A real floating-point type. */
template <typename T>
concept RealFloatingPoint = NumericConcepts::Real<T>;

/** @brief A `std::complex` over a real floating-point type. */
template <typename T>
concept ComplexFloatingPoint = NumericConcepts::Complex<T>;

/** @brief Either a real floating-point type or a complex one. */
template <typename T>
concept RealOrComplexFloatingPoint = NumericConcepts::RealOrComplex<T>;

/** @brief The underlying real type of a possibly-complex scalar. */
using NumericConcepts::RemoveComplex;

/**
 * @brief A callable giving a scalar at a point of the sphere.
 *
 * @details What a field's constructor samples, and what an interpolant
 * models, so that one can be handed straight to the other and remeshing is
 * a single expression.
 *
 * @tparam Function The callable, invoked as `f(theta, phi)`.
 * @tparam Real The precision of the two angles.
 * @tparam Scalar The value type the result must convert to.
 */
template <typename Function, typename Real, typename Scalar>
concept ScalarFunctionS2 = requires(Function f, Real theta, Real phi) {
  requires RealFloatingPoint<Real>;
  requires RealOrComplexFloatingPoint<Scalar>;
  { f(theta, phi) } -> std::convertible_to<Scalar>;
};

}  // namespace GSHTrans

#endif  //  GSH_TRANS_CONCEPTS_GUARD_H
