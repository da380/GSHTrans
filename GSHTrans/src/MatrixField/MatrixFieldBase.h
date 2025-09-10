#ifndef GSH_TRANS_MATRIX_FIELD_BASE_GUARD_H
#define GSH_TRANS_MATRIX_FIELD_BASE_GUARD_H

#include <concepts>
#include <iostream>

#include "../Concepts.h"
#include "../FieldBase.h"
#include "../GridBase.h"
#include "../ScalarField/ScalarFieldBase.h"
#include "CanonicalMatrix.h"

namespace GSHTrans {

template <typename Derived>
class MatrixFieldBase : public FieldBase<MatrixFieldBase<Derived>> {};

}  // namespace GSHTrans

#endif