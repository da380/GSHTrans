#pragma once

/**
 * @file GSHTrans.hpp
 * @brief The whole library.
 *
 * @details This is the header to include. The pieces stack -- layered fields
 * on tensors and expansions, those on spin fields, those on the core -- so
 * "everything" is what nearly every caller needs, and one name for it is
 * easier to remember than a taxonomy of modules.
 *
 * A caller who wants less can have it: `GSHTrans/Core.hpp` is the grids, the
 * transforms, the Wigner functions and the 3-j symbols with none of the field
 * algebra, and every header in the tree is self-sufficient, so any one of
 * them can be included by its own name.
 */

#include "Core.hpp"

// Spin fields: a field of definite upper index on the sphere, and its lazy
// algebra.
#include "SpinField/SpinField.hpp"
#include "SpinField/SpinFieldNodes.hpp"
#include "SpinField/SpinFieldOverloads.hpp"
#include "SpinField/SpinFieldView.hpp"
#include "SpinField/SpinWeighted.hpp"

// Tensor fields: canonical components addressed by multi-index, and the
// storage that follows from their symmetries. A single component is a spin
// field.
#include "Tensor/BundleMaps.hpp"
#include "Tensor/MultiIndex.hpp"
#include "Tensor/Orbits.hpp"
#include "Tensor/TensorExpr.hpp"
#include "Tensor/TensorField.hpp"

// The spectral side: a field's coefficients as an object rather than a raw
// buffer, and the operators that connect different upper indices.
#include "Expansion/BundleMaps.hpp"
#include "Expansion/ContravariantDerivative.hpp"
#include "Expansion/Eth.hpp"
#include "Expansion/Interpolate.hpp"
#include "Expansion/IntrinsicDerivative.hpp"
#include "Expansion/SpinExpansion.hpp"
#include "Expansion/TensorExpansion.hpp"

// Three-dimensional fields: a radial grid crossed with an angular one, held as
// a stack of angular slices. Two-dimensional is the primitive; a slice is an
// ordinary spin field.
#include "Layered/LayeredGradient.hpp"
#include "Layered/LayeredSpinField.hpp"
#include "Layered/LayeredTensorField.hpp"
#include "Layered/RadialDerivatives.hpp"
#include "Layered/RadialGrid.hpp"
#include "Layered/RadialMajor.hpp"
#include "Layered/RadialOperator.hpp"
#include "Layered/RadialResample.hpp"
#include "Layered/RadialSplineDerivative.hpp"
