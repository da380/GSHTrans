#ifndef GSH_TRANS_RADIAL_RESAMPLE_GUARD_H
#define GSH_TRANS_RADIAL_RESAMPLE_GUARD_H

// Resampling a layered field onto a different set of radii.
//
// The whole of this header is conditional on GSHTRANS_WITH_INTERPOLATION,
// which is on by default. Including it without the dependency is not an error
// and gives nothing: the option is what decides whether the facility exists,
// and a caller who has turned it off has said they do not want it.
//
// It is the one place the library uses `Interpolation` at run time. The three
// ready-made derivatives do not, deliberately -- they are called once per
// radial line in the inner loop of a matrix-free solve, where an interpolant
// constructed per line would cost more in allocation than the transform it
// sits beside (field-algebra-plan.md section 19.5 [R4]).
//
// Resampling inverts that argument. Remeshing is something a caller does
// between solves and not inside one, so an interpolant per line is
// affordable, and paying for it buys the whole menu -- linear, natural cubic
// spline, Akima -- for the price of a policy argument rather than three
// implementations.

#ifdef GSHTRANS_HAVE_INTERPOLATION

#include <omp.h>

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <Interpolation/AkimaSpline.hpp>
#include <Interpolation/CubicSpline.hpp>
#include <Interpolation/Linear.hpp>

#include "../Concepts.h"
#include "../Policies.h"
#include "RadialGrid.h"
#include "RadialOperator.h"

namespace GSHTrans {

// Which interpolant carries a radial line onto the new radii.
//
// A policy value with named constructors, as `Execution`, `Batch`, `Chunking`
// and `WignerValues` are, and for the reason Policies.h gives: what it
// describes is a property of the problem rather than of the mathematics, so
// putting it in a type would make the mathematics carry it.
//
// Linear is the safe one -- it cannot overshoot, so a monotone profile stays
// monotone, which for a density or a modulus is often what matters more than
// smoothness. CubicSpline is smooth and can overshoot near a sharp change.
// Akima is the compromise: local, so a bad patch stays local, and much less
// prone to overshoot than a global spline.
class RadialInterpolation {
 public:
  static RadialInterpolation Linear() { return RadialInterpolation(Kind::Line); }
  static RadialInterpolation CubicSpline() {
    return RadialInterpolation(Kind::Cubic);
  }
  static RadialInterpolation Akima() { return RadialInterpolation(Kind::Akima); }

  bool IsLinear() const { return _kind == Kind::Line; }
  bool IsCubicSpline() const { return _kind == Kind::Cubic; }
  bool IsAkima() const { return _kind == Kind::Akima; }

  bool operator==(const RadialInterpolation&) const = default;

 private:
  enum class Kind { Line, Cubic, Akima };
  explicit RadialInterpolation(Kind kind) : _kind{kind} {}
  Kind _kind;
};

namespace ResampleDetails {

// One line, through one interpolant. The interpolant borrows both ranges, so
// it must not outlive them, and here it does not: it is built and used inside
// this call.
template <typename Interpolant, typename Real, typename Scalar>
void Fit(std::span<const Real> from, std::span<const Scalar> values,
         std::span<const Real> onto, std::span<Scalar> out) {
  const Interpolant interpolant{from, values};
  for (std::size_t k = 0; k < onto.size(); k++) {
    out[k] = interpolant(onto[k]);
  }
}

}  // namespace ResampleDetails

// The same field, on different radii.
//
// **Not a RadialOperator, and it could not be one.** An operator maps a stack
// to one of the same shape; this changes the length of the radial axis, so it
// needs its own seam -- `SameShapeOn` rather than `SameShape` -- and cannot go
// through `ApplyRadially`. That is the whole reason it is a free function here
// rather than another entry in RadialDerivatives.h.
//
// Extrapolation is refused rather than silently performed. Every interpolant
// here will happily return a number outside the range it was fitted on, and
// that number is worth nothing; a caller who wants to extend a model beyond
// its outermost radius is doing something the model does not say, and should
// say it themselves.
template <LayeredStack Stack, typename Real = typename std::remove_cvref_t<
                                  decltype(std::declval<const Stack&>()
                                               .Radial())>::Real>
auto Resample(const Stack& in, RadialGrid<Real> onto,
              RadialInterpolation scheme = RadialInterpolation::CubicSpline(),
              Execution policy = Execution::Sequential()) {
  using Int = std::ptrdiff_t;
  using Scalar = typename std::remove_cvref_t<
      decltype(std::declval<const Stack&>().Data())>::value_type;

  const auto from = in.Radial().Radii();
  const auto target = onto.Radii();

  if (target.front() < from.front() || target.back() > from.back()) {
    throw std::invalid_argument(
        "Resampling interpolates and does not extrapolate, and these radii "
        "reach outside the ones the field is given on");
  }

  const auto nOld = in.NumberOfRadii();
  const auto nNew = onto.NumberOfRadii();
  const auto lines = in.SliceSize();
  const auto source = in.Data();

  auto out = in.SameShapeOn(std::move(onto));
  auto destination = out.Data();

  const auto run = [&](Int j) {
    // Gathered per line, as ApplyRadially does and for the same reason: the
    // interpolants want a contiguous range, and the stack is radius-major.
    thread_local auto gathered = std::vector<Scalar>{};
    thread_local auto applied = std::vector<Scalar>{};
    if (gathered.size() < static_cast<std::size_t>(nOld)) {
      gathered.resize(static_cast<std::size_t>(nOld));
    }
    if (applied.size() < static_cast<std::size_t>(nNew)) {
      applied.resize(static_cast<std::size_t>(nNew));
    }

    for (auto i = Int{0}; i < nOld; i++) {
      gathered[static_cast<std::size_t>(i)] =
          source[static_cast<std::size_t>(i * lines + j)];
    }

    const auto values =
        std::span<const Scalar>(gathered.data(), static_cast<std::size_t>(nOld));
    auto answer =
        std::span<Scalar>(applied.data(), static_cast<std::size_t>(nNew));

    // The branch is on the policy and not on the data, so it is the same for
    // every line and costs a predicted jump.
    if (scheme.IsLinear()) {
      ResampleDetails::Fit<Interpolation::Linear<std::span<const Real>,
                                                 std::span<const Scalar>>>(
          from, values, target, answer);
    } else if (scheme.IsAkima()) {
      ResampleDetails::Fit<Interpolation::AkimaSpline<std::span<const Real>,
                                                      std::span<const Scalar>>>(
          from, values, target, answer);
    } else {
      ResampleDetails::Fit<Interpolation::CubicSpline<std::span<const Real>,
                                                      std::span<const Scalar>>>(
          from, values, target, answer);
    }

    for (auto i = Int{0}; i < nNew; i++) {
      destination[static_cast<std::size_t>(i * lines + j)] =
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

  return out;
}

}  // namespace GSHTrans

#endif  // GSHTRANS_HAVE_INTERPOLATION

#endif  // GSH_TRANS_RADIAL_RESAMPLE_GUARD_H
