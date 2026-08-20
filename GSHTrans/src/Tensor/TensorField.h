#ifndef GSH_TRANS_TENSOR_FIELD_GUARD_H
#define GSH_TRANS_TENSOR_FIELD_GUARD_H

#include <FFTWpp/Core>

#include <array>
#include <complex>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

#include "../Concepts.h"
#include "../SpinField/SpinFieldOverloads.h"
#include "../SpinField/SpinFieldView.h"
#include "../SpinField/SpinWeighted.h"
#include "MultiIndex.h"
#include "Orbits.h"

namespace GSHTrans {

// A tensor field on the sphere: one contiguous buffer, handing out its
// canonical components as phase-1 nodes.
//
// One buffer rather than a tuple of separately allocated fields, because the
// components are what get transformed and the transform wants to see them as a
// batch (core-plan.md [C9]). A component is therefore a *view* into the
// buffer, which is why phase 1 made views admissible wherever an owning field
// is, and why operator[] returns by value on every node.
//
// What is stored is one component per orbit of the symmetry group, computed by
// Orbits.h. Everything else is derived: a component related to a stored one by
// a symmetric permutation *is* that same view, one related by an
// antisymmetric permutation is the view negated, and one whose orbit vanishes
// is not representable at all. So Component<...>() returns different types for
// different components and traversal over all of them is a compile-time loop.
// The field-algebra plan records that consequence for phase 4's reality
// reduction; it arrives here already, because antisymmetry has it too.
template <std::ptrdiff_t _Rank, TensorSymmetry<_Rank> _Symmetry,
          TensorReality _Reality, AngularGrid _Grid>
class TensorField {
 public:
  using Int = std::ptrdiff_t;

  static constexpr Int Rank = _Rank;
  using Symmetry = _Symmetry;
  using Reality = _Reality;
  using GridType = _Grid;
  using Real = typename _Grid::Real;
  using Complex = std::complex<Real>;

  // Every component of a tensor is a complex field in phase 2. Reality makes
  // the all-zero component real-valued, and that is phase 4's business: see
  // Orbits.h, where the switch lives.
  using Scalar = Complex;

  static constexpr auto& Orbits = TensorOrbits<Rank, Symmetry, Reality>;
  static constexpr Int Components = MultiIndex<Rank>::Size;
  static constexpr Int StoredComponents = Orbits.storedCount;

  // The buffer's component order, which is **by upper index** and not by flat
  // multi-index.
  //
  // This is a storage decision and it is forced by the transform. A batch
  // shares grid, degree and upper index, and is described by (count, stride,
  // dist) -- a uniform spacing (core-plan.md [C9]). In flat order the stored
  // components carrying one upper index are scattered at no fixed spacing
  // once there is any symmetry, so no single descriptor covers them and the
  // batching that step F exists for would be unreachable. Ordered by upper
  // index they are contiguous, and one batch per upper index describes the
  // whole tensor.
  //
  // Orbits.h stays in flat order, which is pure combinatorics. The layout
  // belongs here.
  struct Layout {
    std::array<Int, StoredComponents> flatOfSlot{};
    std::array<Int, StoredComponents> upperIndexOfSlot{};
    std::array<Int, 2 * Rank + 1> firstSlotAt{};
    std::array<Int, 2 * Rank + 1> countAt{};
  };

  static constexpr Layout ComponentLayout = [] {
    auto layout = Layout{};
    auto slot = Int{0};
    for (auto n = -Rank; n <= Rank; n++) {
      layout.firstSlotAt[n + Rank] = slot;
      for (auto i = Int{0}; i < StoredComponents; i++) {
        const auto flat = Orbits.stored[i];
        if (MultiIndex<Rank>::FromFlat(flat).UpperIndex() != n) continue;
        layout.flatOfSlot[slot] = flat;
        layout.upperIndexOfSlot[slot] = n;
        slot++;
      }
      layout.countAt[n + Rank] = slot - layout.firstSlotAt[n + Rank];
    }
    return layout;
  }();

  // Where a stored component sits in the buffer, by flat multi-index.
  static constexpr Int SlotOfFlat(Int flat) {
    for (auto slot = Int{0}; slot < StoredComponents; slot++) {
      if (ComponentLayout.flatOfSlot[slot] == flat) return slot;
    }
    return -1;
  }

  // The flat component index of a multi-index given as template arguments.
  template <Int... Alphas>
  static constexpr Int FlatOf =
      MultiIndex<Rank>(std::array<Int, Rank>{Alphas...}).Flat();

  // The upper index a component carries, which is the spin weight of its
  // field and so a compile-time quantity (eq:N).
  template <Int... Alphas>
  static constexpr Int UpperIndexOf =
      MultiIndex<Rank>(std::array<Int, Rank>{Alphas...}).UpperIndex();

  // Whether a component vanishes identically, which happens when a
  // permutation maps it to itself with a sign of -1. Exposed so that a
  // compile-time traversal can skip those rather than failing on them.
  //
  // Whether a component can be read at all, and whether it can be written.
  // These are the conditions on the accessors below, named so that a
  // compile-time traversal can ask before it asks for the component, and so
  // that the negative cases can be tested. static_assert(!requires { ... })
  // is the idiom phase 1 established, and it works only when the constraint is
  // a requires-clause: an assertion inside the body is a hard error that no
  // requires-expression can see, which makes the negative test vacuous.
  //
  // The wrong number of indices has to be rejected *before* the multi-index is
  // formed, or the failure is a hard error inside std::array rather than an
  // unsatisfied constraint -- so the pack size is checked with `if constexpr`
  // and the discarded branch never forms one.
  template <Int... Alphas>
  static constexpr bool VanishesFn() {
    if constexpr (sizeof...(Alphas) != Rank) {
      return false;
    } else {
      return Orbits.constraint[FlatOf<Alphas...>] == ComponentConstraint::Zero;
    }
  }

  template <Int... Alphas>
  static constexpr bool RepresentsFn() {
    if constexpr (sizeof...(Alphas) != Rank) {
      return false;
    } else {
      return not VanishesFn<Alphas...>();
    }
  }

  template <Int... Alphas>
  static constexpr bool WritableFn() {
    if constexpr (sizeof...(Alphas) != Rank) {
      return false;
    } else {
      return RepresentsFn<Alphas...>() and Orbits.sign[FlatOf<Alphas...>] == 1;
    }
  }

  template <Int... Alphas>
  static constexpr bool Vanishes = VanishesFn<Alphas...>();

  template <Int... Alphas>
  static constexpr bool Represents = RepresentsFn<Alphas...>();

  template <Int... Alphas>
  static constexpr bool Writable = WritableFn<Alphas...>();

  TensorField() = delete;

  explicit TensorField(GridType grid)
      : _grid{std::move(grid)},
        _data(static_cast<std::size_t>(StoredComponents) *
              static_cast<std::size_t>(_grid.FieldSize())) {
    // A rank-p tensor has components at every upper index from -p to p, so a
    // grid that does not carry them cannot hold one. Checked here rather than
    // left to the first component that asks, which would fail late and name
    // the component instead of the grid.
    if (_grid.MaxUpperIndex() < Rank) {
      throw std::invalid_argument(
          "A rank-" + std::to_string(Rank) +
          " tensor has components at upper index " + std::to_string(Rank) +
          ", but this grid carries only " +
          std::to_string(_grid.MaxUpperIndex()));
    }
  }

  const GridType& Grid() const { return _grid; }

  auto FieldSize() const { return static_cast<Int>(_grid.FieldSize()); }
  auto Size() const { return static_cast<Int>(_data.size()); }

  // The whole buffer, in [component][iTheta][iPhi] order.
  auto Data() { return std::span<Scalar>(_data); }
  auto Data() const { return std::span<const Scalar>(_data); }

  //------------------------------------------------------------------------//
  //                            Component access                             //
  //------------------------------------------------------------------------//

  // The component with this multi-index, as a phase-1 node.
  //
  // Read-only, and available for every component the tensor can represent:
  // stored ones as a view, sign-reversed ones as that view negated.
  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const {
    constexpr auto flat = FlatOf<Alphas...>;
    constexpr auto N = UpperIndexOf<Alphas...>;
    constexpr auto sign = Orbits.sign[flat];

    auto view = ConstSpinFieldView<N, GridType, ComplexValued>(
        _grid, StoredSpan<flat>());

    // Negated with std::move so that the expression node owns the view rather
    // than referring to this local one. A phase-1 node holds an lvalue
    // terminal by reference, which is exactly right at a call site and
    // exactly wrong here.
    if constexpr (sign == 1) {
      return view;
    } else {
      return -std::move(view);
    }
  }

  // The component as writable storage.
  //
  // Only where the component *is* the stored one up to a symmetric
  // permutation. Writing through a sign-reversed component would mean
  // negating on the way in, and a view cannot express that -- so rather than
  // returning a proxy that behaves unlike every other node in the library,
  // this is a compile error naming the component to write instead. For a
  // symmetric tensor every representable component is writable, which is the
  // common case.
  template <Int... Alphas>
  requires Writable<Alphas...>
  auto Component() {
    constexpr auto flat = FlatOf<Alphas...>;
    constexpr auto N = UpperIndexOf<Alphas...>;
    return SpinFieldView<N, GridType, ComplexValued>(_grid,
                                                     StoredSpan<flat>());
  }

  //------------------------------------------------------------------------//
  //                    What the transform layer will need                   //
  //------------------------------------------------------------------------//

  // Which stored components carry a given upper index: a contiguous run of
  // slots, which is what makes them a batch. Returns (first slot, count).
  static constexpr auto StoredAtUpperIndex(Int n) {
    if (n < -Rank || n > Rank) return std::pair(Int{0}, Int{0});
    return std::pair(ComponentLayout.firstSlotAt[n + Rank],
                     ComponentLayout.countAt[n + Rank]);
  }

  //------------------------------------------------------------------------//
  //                             The transform                               //
  //------------------------------------------------------------------------//

  // How many coefficients a transform at this degree produces: one block per
  // stored component, each sized by that component's upper index, in the same
  // order as the components themselves.
  //
  // The blocks are not all the same length, since the coefficient count
  // depends on the upper index. That is why this is computed rather than
  // being StoredComponents times something.
  Int CoefficientSize(Int lMax) const {
    auto total = Int{0};
    for (auto slot = Int{0}; slot < StoredComponents; slot++) {
      total += static_cast<Int>(
          _grid.CoefficientSize(lMax, ComponentLayout.upperIndexOfSlot[slot]));
    }
    return total;
  }

  // Transform every stored component, batching those that share an upper
  // index.
  //
  // This is the first consumer of the batched primitive step F was built for,
  // and the reason the buffer is ordered by upper index: each group is a
  // contiguous run on both sides, so one Batch::Contiguous describes it and
  // the Wigner block for that upper index is streamed once for the whole
  // group rather than once per component. A rank-2 tensor has three
  // components at N = 0, so that is three fields for the price of one pass.
  //
  // Derived components are not transformed. They are determined by the stored
  // ones, and transforming them would be doing the same work twice and
  // storing the answer twice.
  void ForwardTransformation(Int lMax, std::span<Complex> out,
                             Execution policy = Execution::Sequential()) const {
    CheckCoefficients(out.size(), lMax);
    const auto fieldSize = FieldSize();
    auto offset = std::size_t{0};
    for (auto n = -Rank; n <= Rank; n++) {
      const auto [first, count] = StoredAtUpperIndex(n);
      if (count == 0) continue;
      const auto coefficientSize =
          static_cast<Int>(_grid.CoefficientSize(lMax, n));
      auto fields = Data().subspan(static_cast<std::size_t>(first * fieldSize),
                                   static_cast<std::size_t>(count * fieldSize));
      auto block =
          out.subspan(offset, static_cast<std::size_t>(count * coefficientSize));
      _grid.ForwardTransformation(lMax, n, fields,
                                  Batch::Contiguous(count, fieldSize), block,
                                  Batch::Contiguous(count, coefficientSize),
                                  policy);
      offset += static_cast<std::size_t>(count * coefficientSize);
    }
  }

  void InverseTransformation(Int lMax, std::span<const Complex> in,
                             Execution policy = Execution::Sequential()) {
    CheckCoefficients(in.size(), lMax);
    const auto fieldSize = FieldSize();
    auto offset = std::size_t{0};
    for (auto n = -Rank; n <= Rank; n++) {
      const auto [first, count] = StoredAtUpperIndex(n);
      if (count == 0) continue;
      const auto coefficientSize =
          static_cast<Int>(_grid.CoefficientSize(lMax, n));
      auto block =
          in.subspan(offset, static_cast<std::size_t>(count * coefficientSize));
      auto fields = Data().subspan(static_cast<std::size_t>(first * fieldSize),
                                   static_cast<std::size_t>(count * fieldSize));
      _grid.InverseTransformation(lMax, n, block,
                                  Batch::Contiguous(count, coefficientSize),
                                  fields, Batch::Contiguous(count, fieldSize),
                                  policy);
      offset += static_cast<std::size_t>(count * coefficientSize);
    }
  }

 private:
  GridType _grid;
  FFTWpp::vector<Scalar> _data;

  void CheckCoefficients(std::size_t given, Int lMax) const {
    const auto needed = static_cast<std::size_t>(CoefficientSize(lMax));
    if (given != needed) {
      throw std::invalid_argument(
          "Tensor coefficient range has size " + std::to_string(given) +
          ", but this tensor needs " + std::to_string(needed) + " at degree " +
          std::to_string(lMax));
    }
  }

  // The stored component's samples. The representative is looked up at
  // compile time; only the multiplication by the field size is left to run
  // time, and that because the grid is a runtime object.
  template <Int Flat>
  std::span<Scalar> StoredSpan() {
    constexpr auto slot = SlotOfFlat(Orbits.representative[Flat]);
    static_assert(slot >= 0);
    const auto size = static_cast<std::size_t>(_grid.FieldSize());
    return std::span<Scalar>(_data).subspan(
        static_cast<std::size_t>(slot) * size, size);
  }

  template <Int Flat>
  std::span<const Scalar> StoredSpan() const {
    constexpr auto slot = SlotOfFlat(Orbits.representative[Flat]);
    static_assert(slot >= 0);
    const auto size = static_cast<std::size_t>(_grid.FieldSize());
    return std::span<const Scalar>(_data).subspan(
        static_cast<std::size_t>(slot) * size, size);
  }
};

// The ranks the library exposes, named for readability at the call site.
template <TensorReality Reality, AngularGrid Grid>
using VectorField = TensorField<1, NoSymmetry<1>, Reality, Grid>;

template <typename Symmetry, TensorReality Reality, AngularGrid Grid>
using Rank2TensorField = TensorField<2, Symmetry, Reality, Grid>;

template <TensorReality Reality, AngularGrid Grid>
using ElasticTensorField = TensorField<4, ElasticSymmetry, Reality, Grid>;

}  // namespace GSHTrans

#endif  // GSH_TRANS_TENSOR_FIELD_GUARD_H
