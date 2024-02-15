#ifndef GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H
#define GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H

#include <FFTWpp/All>
#include <algorithm>
#include <complex>
#include <iostream>
#include <memory>
#include <ranges>

#include "Concepts.h"
#include "Grid.h"
#include "Indexing.h"

namespace GSHTrans {

template <typename Derived>
class CanonicalCoefficientBase
    : public std::ranges::view_interface<CanonicalCoefficientBase<Derived>> {
 public:
 private:
  auto& _Derived() const { return static_cast<const Derived&>(*this); }
  auto& _Derived() { return static_cast<Derived&>(*this); }
};

template <typename Grid, RealOrComplexValued Type>
requires std::derived_from<Grid, GridBase<Grid>>
class CanonicalCoefficient
    : public CanonicalCoefficientBase<CanonicalCoefficient<Grid, Type>>,
      public GSHIndices<std::conditional_t<std::same_as<Type, ComplexValued>,
                                           All, NonNegative>>

{};

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H
