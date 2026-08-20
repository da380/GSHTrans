#ifndef GSH_TRANS_RADIAL_OPERATOR_GUARD_H
#define GSH_TRANS_RADIAL_OPERATOR_GUARD_H

#include <omp.h>

#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "../Concepts.h"
#include "LayeredSpinField.h"
#include "RadialGrid.h"

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                              The radial seam                              //
//--------------------------------------------------------------------------//

// Where this library stops and the application's radial discretisation begins.
//
// The angular half is the library's: the transform, the index algebra, the
// surface gradient. The radial half is not, and deliberately so. A
// finite-difference derivative is banded, a spectral-element one is
// block-diagonal, a solve may carry a factorisation the caller wants to reuse,
// and none of that is a spherical-harmonic question. What the library owns is
// the *seam*: gathering a radial line out of a radius-major stack, handing it
// to whatever the caller supplied, and scattering the answer back.
//
// The seam is a callable rather than a matrix, which is the choice
// field-algebra-plan.md section 17.2 settled and for the reason it gives: a
// caller who has factorised something wants to apply the factorisation, and a
// matrix interface forecloses that. It takes contiguous spans and not strided
// ones, because a contiguous span is what LAPACK, a band solver and a plain
// loop all want, and the gather that buys it is a fraction of the transform it
// sits beside.
//
// A radial operator maps one radial line to another. Both spans have length
// nR, and it must not assume they alias or that they do not.
template <typename Op, typename Scalar>
concept RadialOperator =
    requires(const Op& op, std::span<const Scalar> in, std::span<Scalar> out) {
      op(in, out);
    };

// The two stack types present the radial axis under the same names, and this
// is all of it that a radial operator sees: nR lines, SliceSize() apart. A
// field's slice holds angular points and an expansion's holds coefficients;
// the radial axis does not distinguish them, so neither does this.
template <typename Stack>
concept LayeredStack = requires(const Stack& stack) {
  { stack.NumberOfRadii() } -> std::convertible_to<std::ptrdiff_t>;
  { stack.SliceSize() } -> std::convertible_to<std::ptrdiff_t>;
  { stack.Radial() };
  { stack.Data() };
  { stack.SameShape() } -> std::same_as<Stack>;
};

// Apply a radial operator at every angular index, or at every coefficient.
//
// `in` and `out` must be distinct objects: each line is gathered before the
// operator runs and scattered after, so an in-place call would be correct only
// by accident of the gather buffer, and requiring distinctness says so rather
// than relying on it. Use the returning form, which allocates the result.
//
// Threading is over lines, which is the axis with the most of them and the one
// with no dependence between iterations. The gather buffers are thread-local
// and grow to fit, for the same reason the transform's work buffers are: the
// alternative is an allocation per line.
template <LayeredStack Stack, typename Op>
requires RadialOperator<Op, typename std::remove_cvref_t<
                                decltype(std::declval<Stack&>().Data())>::
                                element_type>
void ApplyRadially(const Stack& in, Stack& out, const Op& op,
                   Execution policy = Execution::Sequential()) {
  using Int = std::ptrdiff_t;
  using Scalar = typename std::remove_cvref_t<
      decltype(std::declval<Stack&>().Data())>::element_type;

  if (in.Radial().Identity() != out.Radial().Identity()) {
    throw std::invalid_argument(
        "A radial operator maps a stack to one on the same radial grid");
  }
  if (in.SliceSize() != out.SliceSize()) {
    throw std::invalid_argument(
        "A radial operator acts along the radial axis alone, so the slices of "
        "its argument and its result must be the same size");
  }
  if (static_cast<const void*>(in.Data().data()) ==
      static_cast<const void*>(out.Data().data())) {
    throw std::invalid_argument(
        "ApplyRadially gathers each radial line before applying the operator, "
        "so its argument and result must be distinct");
  }

  const auto nR = in.NumberOfRadii();
  const auto lines = in.SliceSize();
  const auto source = in.Data();
  auto target = out.Data();

  const auto run = [&](Int j) {
    thread_local auto gathered = std::vector<Scalar>{};
    thread_local auto applied = std::vector<Scalar>{};
    const auto n = static_cast<std::size_t>(nR);
    if (gathered.size() < n) gathered.resize(n);
    if (applied.size() < n) applied.resize(n);

    for (auto i = Int{0}; i < nR; i++) {
      gathered[static_cast<std::size_t>(i)] =
          source[static_cast<std::size_t>(i * lines + j)];
    }
    op(std::span<const Scalar>(gathered.data(), n),
       std::span<Scalar>(applied.data(), n));
    for (auto i = Int{0}; i < nR; i++) {
      target[static_cast<std::size_t>(i * lines + j)] =
          applied[static_cast<std::size_t>(i)];
    }
  };

  const auto threads =
      policy.IsParallel() && !omp_in_parallel()
          ? (policy.Threads() > 0 ? policy.Threads() : omp_get_max_threads())
          : 1;

  if (threads == 1) {
    for (auto j = Int{0}; j < lines; j++) run(j);
  } else {
#pragma omp parallel for schedule(static) num_threads(threads)
    for (Int j = 0; j < lines; j++) run(j);
  }
}

template <LayeredStack Stack, typename Op>
auto ApplyRadially(const Stack& in, const Op& op,
                   Execution policy = Execution::Sequential()) {
  auto out = in.SameShape();
  ApplyRadially(in, out, op, policy);
  return out;
}

//--------------------------------------------------------------------------//
//                          Quadrature over the radii                        //
//--------------------------------------------------------------------------//

// Integrate a stack over radius, leaving one slice.
//
// The weights are whatever rule the caller built the radial grid with, applied
// as they stand. In particular the r^2 of the volume element is *not* inserted
// here: a caller integrating over the ball folds it into the weights, and one
// integrating a radial profile does not want it. Guessing which is meant would
// be wrong half the time, so the grid's weights are taken literally.
//
// The result is a plain buffer of SliceSize() values, because what a slice
// *is* differs between the two stack types and the sum does not care.
template <LayeredStack Stack>
auto IntegrateRadially(const Stack& stack) {
  using Int = std::ptrdiff_t;
  using Scalar = typename std::remove_cvref_t<
      decltype(std::declval<const Stack&>().Data())>::value_type;

  if (!stack.Radial().HasWeights()) {
    throw std::invalid_argument(
        "Integrating over radius needs the radial grid's weights, and this "
        "grid carries points alone");
  }

  const auto nR = stack.NumberOfRadii();
  const auto lines = stack.SliceSize();
  const auto data = stack.Data();
  const auto weights = stack.Radial().Weights();

  auto total = std::vector<Scalar>(static_cast<std::size_t>(lines), Scalar{});
  for (auto i = Int{0}; i < nR; i++) {
    const auto w = weights[static_cast<std::size_t>(i)];
    for (auto j = Int{0}; j < lines; j++) {
      total[static_cast<std::size_t>(j)] +=
          w * data[static_cast<std::size_t>(i * lines + j)];
    }
  }
  return total;
}

}  // namespace GSHTrans

#endif  // GSH_TRANS_RADIAL_OPERATOR_GUARD_H
