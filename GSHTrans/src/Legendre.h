#ifndef GHS_TRANS_LEGNEDRE_GUARD_H
#define GHS_TRANS_LEGNEDRE_GUARD_H

#include <concepts>

#include "Wigner.h"

namespace GSHTrans {

template <std::floating_point Real, OrderRange Orders,
          Normalisation Norm = Ortho>
class AssociatedLegendre : public Wigner<Real, Orders, Norm> {
 public:
  using typename Wigner<Real, Orders, Norm>::const_iterator;
  using typename Wigner<Real, Orders, Norm>::difference_type;
  using typename Wigner<Real, Orders, Norm>::iterator;
  using typename Wigner<Real, Orders, Norm>::size_type;
  using typename Wigner<Real, Orders, Norm>::value_type;

  AssociatedLegendre() = default;

  AssociatedLegendre(difference_type lMax, difference_type mMax, Real theta)
      : Wigner<Real, Orders, Norm>(lMax, mMax, 0, theta) {}

  template <RealFloatingPointIterator Iterator>
  AssociatedLegendre(difference_type lMax, difference_type mMax,
                     Iterator thetaStart, Iterator thetaFinish)
      : Wigner<Real, Orders, Norm>(lMax, mMax, 0, thetaStart, thetaFinish) {}

  template <RealFloatingPointRange Range>
  AssociatedLegendre(difference_type lMax, difference_type mMax,
                     difference_type n, Range&& theta)
      : Wigner<Real, Orders, Norm>(lMax, mMax, 0, theta) {}
};

template <std::floating_point Real, Normalisation Norm = Ortho>
class Legendre : public Wigner<Real, All, Norm> {
 public:
  using typename Wigner<Real, All, Norm>::const_iterator;
  using typename Wigner<Real, All, Norm>::difference_type;
  using typename Wigner<Real, All, Norm>::iterator;
  using typename Wigner<Real, All, Norm>::size_type;
  using typename Wigner<Real, All, Norm>::value_type;

  Legendre() = default;

  Legendre(difference_type lMax, Real theta)
      : Wigner<Real, All, Norm>(lMax, 0, 0, theta) {}

  template <RealFloatingPointIterator Iterator>
  Legendre(difference_type lMax, Iterator thetaStart, Iterator thetaFinish)
      : Wigner<Real, All, Norm>(lMax, 0, 0, thetaStart, thetaFinish) {}

  template <RealFloatingPointRange Range>
  Legendre(difference_type lMax, difference_type n, Range&& theta)
      : Wigner<Real, All, Norm>(lMax, 0, 0, theta) {}

  auto operator()(difference_type l) {
    return Wigner<Real, All, Norm>::operator()(l, 0);
  }

  auto operator()(difference_type i, difference_type l) {
    return Wigner<Real, All, Norm>::operator()(i, l, 0);
  }
};

}  // namespace GSHTrans

#endif  // GHS_TRANS_LEGNEDRE_GUARD_H
