#ifndef GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
#define GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H

#include <GaussQuad/All>

#include "Concepts.h"
#include "Indexing.h"
#include "Legendre.h"
#include "Wigner.h"

namespace GSHTrans {

template <RealFloatingPoint Real = double, TransformType Type = C2C>
class GaussLegendreGrid {
 public:
  GaussLegendreGrid() = default;

  GaussLegendreGrid(int lMax, int nMax)
      : _lMax{lMax},
        _nMax{nMax},
        _quad{GaussQuad::GaussLegendreQuadrature1D<Real>(_lMax + 1)} {}

 private:
  // Store maximum orders and upper indices.
  int _lMax;
  int _nMax;

  // Store quadrature points and weights.
  GaussQuad::Quadrature1D<Real> _quad;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
