#ifndef GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
#define GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H

#include <GaussQuad/All>

#include "Concepts.h"
#include "Indexing.h"
#include "Legendre.h"
#include "Wigner.h"

namespace GSHTrans {

template <RealFloatingPoint Real, TransformType Type,
          Normalisation Norm = Ortho>
class GaussLegendreGrid {
  using MRange = Type::IndexRange;
  using NRange = Type::IndexRange;
  using WignerType = Wigner<Real, MRange, NRange, Norm>;

 public:
  GaussLegendreGrid() = default;

  GaussLegendreGrid(int lMax, int nMax)
      : _lMax{lMax},
        _nMax{nMax},
        _quad{GaussQuad::GaussLegendreQuadrature1D<Real>(_lMax + 1)} {
    assert(MaxDegree() >= 0);
    assert(MaxUpperIndex() >= 0 && MaxUpperIndex() <= MaxDegree());

    // Transform the quadrature points.
    _quad.Transform([](auto x) { return std::acos(-x); },
                    [](auto x) { return 1; });

    // Compute the Wigner values
    _wigner = std::make_shared<WignerType>(_lMax, _lMax, _nMax, _quad.Points());
  }

  auto MaxDegree() const { return _lMax; }
  auto MaxUpperIndex() const { return _nMax; }

 private:
  // Store maximum orders and upper indices.
  int _lMax;
  int _nMax;

  // Quadrature points and weights.
  GaussQuad::Quadrature1D<Real> _quad;

  // Wigner values
  std::shared_ptr<WignerType> _wigner;
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_GAUSS_LEGENDRE_GRID_GUARD_H
