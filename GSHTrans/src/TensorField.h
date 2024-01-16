#ifndef GSH_TRANS_TENSOR_FIELD_GUARD_H
#define GSH_TRANS_TENSOR_FIELD_GUARD_H

#include <Eigen/Core>

#include "Concepts.h"
#include "Indexing.h"

namespace GSHTrans {}  // namespace GSHTrans

template <typename Grid>
class ScalarField {
  using Scalar = Grid::scalar_type;
  using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

 public:
  ScalarField(std::shared_ptr<Grid> grid)
      : _grid{grid},
        _data{Vector(_grid->NumberOfCoLatitudes() *
                     _grid->NumberOfLongitudes())} {}

 private:
  std::shared_ptr<Grid> _grid;
  Vector _data;
};

#endif  // GSH_TRANS_TENSOR_FIELD_GUARD_H
