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
#include "../Policies.h"
#include "../SpinField/SpinFieldOverloads.h"
#include "../SpinField/SpinFieldView.h"
#include "../SpinField/SpinWeighted.h"
#include "MultiIndex.h"
#include "Orbits.h"

namespace GSHTrans {

// How a tensor's components are arranged in its buffer.
//
// ComponentMajor keeps each component contiguous, which is what a transform
// wants and what a component-at-a-time traversal wants. PointMajor keeps all
// of a point's components together, which is what applying a rank-4 tensor
// pointwise to a rank-2 one wants -- the operation reads 81 numbers at one
// point and none at any other.
//
// Both are offered rather than one being required, because the transform's
// batch descriptor covers either: (stride 1, dist FieldSize) for the first
// and (stride nStored, dist 1) for the second (core-plan.md [C9]). Which is
// faster is a measurement, not a precondition, and the field-algebra plan
// dropped the repack requirement on the strength of exactly that.
struct ComponentMajor {};
struct PointMajor {};

template <typename L>
concept TensorLayout =
    std::same_as<L, ComponentMajor> or std::same_as<L, PointMajor>;

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
          TensorReality _Reality, AngularGrid _Grid,
          TensorLayout _Layout = ComponentMajor>
class TensorField {
 public:
  using Int = std::ptrdiff_t;

  static constexpr Int Rank = _Rank;
  using Symmetry = _Symmetry;
  using Reality = _Reality;
  using GridType = _Grid;
  using LayoutPolicy = _Layout;

  static constexpr bool IsComponentMajor =
      std::same_as<_Layout, ComponentMajor>;
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
  // Reality splits the storage in two, and that is the other thing the layout
  // has to carry.
  //
  // Once the reality condition is one of the generators, an orbit that
  // contains its own negation pins its component to a single real number --
  // real, or purely imaginary if a permutation sign gets in the way. Storing
  // those as complex fields would waste half of each and, worse, would let a
  // caller write an imaginary part into a component that cannot have one. So
  // they live in a buffer of their own.
  //
  // They are always at upper index zero: permutation preserves the slot sum
  // and negation reverses it, so a component fixed by the combination
  // satisfies N = -N. That is what makes a real-valued field admissible for
  // them at all -- phase 1 forbids RealValued anywhere else -- and it means
  // the real buffer is one transform group rather than several.
  //
  // The arithmetic comes out exactly right: two reals per complex component
  // and one per constrained one is 3^p, the real degrees of freedom of a real
  // rank-p tensor. Nine for rank 2, six symmetric, three antisymmetric, ten
  // for symmetric rank 3.
  struct Layout {
    std::array<Int, StoredComponents> flatOfSlot{};
    std::array<Int, StoredComponents> upperIndexOfSlot{};
    std::array<bool, StoredComponents> realOfSlot{};
    std::array<Int, 2 * Rank + 1> firstSlotAt{};
    std::array<Int, 2 * Rank + 1> countAt{};
    Int complexCount{};
    Int realCount{};
  };

  static constexpr Layout ComponentLayout = [] {
    auto layout = Layout{};
    const auto constrained = [](Int flat) {
      return Orbits.constraint[flat] != ComponentConstraint::None;
    };

    // The complex components first, grouped by upper index so that each group
    // is a batch, then the constrained ones, which are all at zero.
    auto slot = Int{0};
    for (auto n = -Rank; n <= Rank; n++) {
      layout.firstSlotAt[n + Rank] = slot;
      for (auto i = Int{0}; i < StoredComponents; i++) {
        const auto flat = Orbits.stored[i];
        if (constrained(flat)) continue;
        if (MultiIndex<Rank>::FromFlat(flat).UpperIndex() != n) continue;
        layout.flatOfSlot[slot] = flat;
        layout.upperIndexOfSlot[slot] = n;
        layout.realOfSlot[slot] = false;
        slot++;
      }
      layout.countAt[n + Rank] = slot - layout.firstSlotAt[n + Rank];
    }
    layout.complexCount = slot;

    for (auto i = Int{0}; i < StoredComponents; i++) {
      const auto flat = Orbits.stored[i];
      if (!constrained(flat)) continue;
      layout.flatOfSlot[slot] = flat;
      layout.upperIndexOfSlot[slot] = 0;
      layout.realOfSlot[slot] = true;
      slot++;
    }
    layout.realCount = slot - layout.complexCount;
    return layout;
  }();

  static constexpr Int ComplexComponents = ComponentLayout.complexCount;
  static constexpr Int RealComponents = ComponentLayout.realCount;

  // The real numbers one grid point of this tensor costs.
  static constexpr Int RealsPerPoint = 2 * ComplexComponents + RealComponents;

  static constexpr Int UpperIndexOfFlat(Int flat) {
    return MultiIndex<Rank>::FromFlat(flat).UpperIndex();
  }

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
      return RepresentsFn<Alphas...>() and
             Orbits.sign[FlatOf<Alphas...>] == 1 and
             not Orbits.conjugate[FlatOf<Alphas...>];
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
        _data(static_cast<std::size_t>(ComplexComponents) *
              static_cast<std::size_t>(_grid.FieldSize())),
        _real(static_cast<std::size_t>(RealComponents) *
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

  // The complex buffer, in [component][iTheta][iPhi] order, and the real one
  // holding the components the reality condition pins to a single real number.
  // The second is empty unless reality is being reduced on.
  auto Size() const { return static_cast<Int>(_data.size()); }
  auto RealSize() const { return static_cast<Int>(_real.size()); }

  auto Data() { return std::span<Scalar>(_data); }
  auto Data() const { return std::span<const Scalar>(_data); }
  auto RealData() { return std::span<Real>(_real); }
  auto RealData() const { return std::span<const Real>(_real); }

  //------------------------------------------------------------------------//
  //                            Component access                             //
  //------------------------------------------------------------------------//

  // The component with this multi-index, as a phase-1 node.
  //
  // Read-only, and available for every component the tensor can represent.
  // What comes back depends on how the component is related to the one stored
  // for it, and there are now four cases rather than two:
  //
  //   stored, sign +1        the view itself
  //   sign -1               that view negated, an expression
  //   related by reality    conj of the view, at the reversed upper index
  //   pinned by its orbit   a real-valued view, or i times one
  //
  // The reality case is the one phase 1 was made to accommodate. The relation
  // is T^{-alpha} = (-1)^N conj(T^{alpha}) (eq:reality), and conj reverses the
  // upper index -- which is why getting that wrong was one of the three
  // structural defects in the layer phase 1 replaced. The (-1)^N is already
  // folded into the orbit table's sign.
  //
  // Note what this means for a grid: the derived partner of a stored
  // component at N has upper index -N, so on a grid carrying only
  // non-negative upper indices half of these could not be *terminals*. They
  // are expressions, and the field-algebra plan's section 3.7 put the grid's
  // N-support check on terminals and views alone for exactly this case.
  template <Int... Alphas>
  requires Represents<Alphas...>
  auto Component() const {
    constexpr auto flat = FlatOf<Alphas...>;
    constexpr auto sign = Orbits.sign[flat];
    constexpr auto conjugated = Orbits.conjugate[flat];
    constexpr auto constraint = Orbits.constraint[flat];
    constexpr auto scale = static_cast<Real>(sign);

    if constexpr (constraint == ComponentConstraint::None) {
      // The stored component's own upper index, which is this one's when the
      // relation does not conjugate and its negative when it does.
      constexpr auto stored = UpperIndexOfFlat(Orbits.representative[flat]);
      auto view = ConstSpinFieldView<stored, GridType, ComplexValued>(
          _grid, StoredSpan<flat>(), ComponentStride);

      // Moved in so that the expression node owns the view rather than
      // referring to this local one. A phase-1 node holds an lvalue terminal
      // by reference, which is right at a call site and wrong here.
      if constexpr (conjugated) {
        return scale * conj(std::move(view));
      } else if constexpr (sign == 1) {
        return view;
      } else {
        return -std::move(view);
      }
    } else {
      // A pinned component: one real number per point. Real means the value
      // is that number; Imaginary means it is i times it, which is what an
      // antisymmetric real tensor's self-paired component is.
      auto view = ConstSpinFieldView<0, GridType, RealValued>(
          _grid, RealStoredSpan<flat>(), RealComponentStride);
      constexpr auto turn = conjugated ? -scale : scale;
      if constexpr (constraint == ComponentConstraint::Real) {
        return turn * std::move(view);
      } else {
        return Complex{0, turn} * std::move(view);
      }
    }
  }

  // The component as writable storage.
  //
  // Only where the component *is* the stored one: writing through a
  // sign-reversed or conjugated component would mean negating or conjugating
  // on the way in, and a view cannot express that -- so rather than returning
  // a proxy that behaves unlike every other node in the library, this is a
  // compile error naming the component to write instead.
  //
  // A pinned component is writable as the single real number it is. For an
  // Imaginary one that number is the coefficient of i, which is the only
  // sensible reading and is why the accessor's value kind says RealValued.
  template <Int... Alphas>
  requires Writable<Alphas...>
  auto Component() {
    constexpr auto flat = FlatOf<Alphas...>;
    constexpr auto N = UpperIndexOf<Alphas...>;
    constexpr auto constraint = Orbits.constraint[flat];
    if constexpr (constraint == ComponentConstraint::None) {
      return SpinFieldView<N, GridType, ComplexValued>(_grid, StoredSpan<flat>(),
                                                       ComponentStride);
    } else {
      return SpinFieldView<0, GridType, RealValued>(
          _grid, RealStoredSpan<flat>(), RealComponentStride);
    }
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
    for (auto slot = Int{0}; slot < ComplexComponents; slot++) {
      total += static_cast<Int>(
          _grid.CoefficientSize(lMax, ComponentLayout.upperIndexOfSlot[slot]));
    }
    // A pinned component is a real field, so it uses the reduced m >= 0
    // storage -- which halves its coefficients too, and is the same saving in
    // the spectral domain that the real buffer is in the spatial one.
    total += RealComponents * static_cast<Int>(_grid.RealCoefficientSize(lMax));
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
      auto fields = FieldGroup(first, count);
      auto block =
          out.subspan(offset, static_cast<std::size_t>(count * coefficientSize));
      _grid.ForwardTransformation(lMax, n, fields, FieldBatch(count), block,
                                  Batch::Contiguous(count, coefficientSize),
                                  policy);
      offset += static_cast<std::size_t>(count * coefficientSize);
    }

    if constexpr (RealComponents > 0) {
      // One more batch: the pinned components, all at upper index zero and
      // all real-valued, so they go through the transform's real path and
      // land in its reduced m >= 0 storage.
      const auto coefficientSize =
          static_cast<Int>(_grid.RealCoefficientSize(lMax));
      auto fields = RealFieldGroup();
      auto block = out.subspan(
          offset, static_cast<std::size_t>(RealComponents * coefficientSize));
      _grid.ForwardTransformation(lMax, 0, fields, RealFieldBatch(), block,
                                  Batch::Contiguous(RealComponents,
                                                    coefficientSize),
                                  policy);
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
      auto fields = FieldGroup(first, count);
      _grid.InverseTransformation(lMax, n, block,
                                  Batch::Contiguous(count, coefficientSize),
                                  fields, FieldBatch(count), policy);
      offset += static_cast<std::size_t>(count * coefficientSize);
    }

    if constexpr (RealComponents > 0) {
      const auto coefficientSize =
          static_cast<Int>(_grid.RealCoefficientSize(lMax));
      auto block = in.subspan(
          offset, static_cast<std::size_t>(RealComponents * coefficientSize));
      auto fields = RealFieldGroup();
      _grid.InverseTransformation(lMax, 0, block,
                                  Batch::Contiguous(RealComponents,
                                                    coefficientSize),
                                  fields, RealFieldBatch(), policy);
    }
  }

 private:
  GridType _grid;
  FFTWpp::vector<Scalar> _data;
  FFTWpp::vector<Real> _real;

  // The window of the buffer holding a run of `count` components starting at
  // `first`, and the descriptor that reads them.
  //
  // This is where the layout stops mattering to anything downstream. Both
  // arrangements are one Batch: components laid end to end are (stride 1,
  // dist FieldSize), and components interleaved point by point are (stride
  // nStored, dist 1) with the run's own components picked out of the wider
  // interleaving -- which is exactly the case Batch::Interleaved documents,
  // where the count is smaller than the stride and the other components are
  // no business of the call.
  auto FieldGroup(Int first, Int count) const {
    const auto fieldSize = FieldSize();
    if constexpr (IsComponentMajor) {
      return Data().subspan(static_cast<std::size_t>(first * fieldSize),
                            static_cast<std::size_t>(count * fieldSize));
    } else {
      return Data().subspan(
          static_cast<std::size_t>(first),
          static_cast<std::size_t>((fieldSize - 1) * ComplexComponents +
                                   count));
    }
  }

  auto FieldGroup(Int first, Int count) {
    const auto fieldSize = FieldSize();
    if constexpr (IsComponentMajor) {
      return Data().subspan(static_cast<std::size_t>(first * fieldSize),
                            static_cast<std::size_t>(count * fieldSize));
    } else {
      return Data().subspan(
          static_cast<std::size_t>(first),
          static_cast<std::size_t>((fieldSize - 1) * ComplexComponents +
                                   count));
    }
  }

  // The pinned components' window and descriptor. They are one group, since
  // they all sit at upper index zero.
  auto RealFieldGroup() {
    const auto fieldSize = FieldSize();
    if constexpr (IsComponentMajor) {
      return RealData().subspan(
          0, static_cast<std::size_t>(RealComponents * fieldSize));
    } else {
      return RealData().subspan(
          0, static_cast<std::size_t>((fieldSize - 1) * RealComponents +
                                      RealComponents));
    }
  }

  auto RealFieldGroup() const {
    const auto fieldSize = FieldSize();
    if constexpr (IsComponentMajor) {
      return RealData().subspan(
          0, static_cast<std::size_t>(RealComponents * fieldSize));
    } else {
      return RealData().subspan(
          0, static_cast<std::size_t>((fieldSize - 1) * RealComponents +
                                      RealComponents));
    }
  }

  Batch RealFieldBatch() const {
    if constexpr (IsComponentMajor) {
      return Batch::Contiguous(RealComponents, FieldSize());
    } else {
      return Batch::Interleaved(RealComponents, RealComponents);
    }
  }

  Batch FieldBatch(Int count) const {
    if constexpr (IsComponentMajor) {
      return Batch::Contiguous(count, FieldSize());
    } else {
      return Batch::Interleaved(count, ComplexComponents);
    }
  }

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
  // Where a stored component's samples live, and how far apart. Contiguous
  // and one apart in ComponentMajor; starting at the component's slot and one
  // buffer-width apart in PointMajor.
  //
  // The two buffers interleave independently, since they hold different
  // numbers of components.
  static constexpr Int ComponentStride =
      IsComponentMajor ? Int{1} : ComplexComponents;
  static constexpr Int RealComponentStride =
      IsComponentMajor ? Int{1} : RealComponents;

  template <Int Flat>
  std::span<Scalar> StoredSpan() {
    return std::span<Scalar>(_data).subspan(SpanOffset<Flat, false>(),
                                            SpanExtent(ComponentStride));
  }

  template <Int Flat>
  std::span<const Scalar> StoredSpan() const {
    return std::span<const Scalar>(_data).subspan(SpanOffset<Flat, false>(),
                                                  SpanExtent(ComponentStride));
  }

  template <Int Flat>
  std::span<Real> RealStoredSpan() {
    return std::span<Real>(_real).subspan(SpanOffset<Flat, true>(),
                                          SpanExtent(RealComponentStride));
  }

  template <Int Flat>
  std::span<const Real> RealStoredSpan() const {
    return std::span<const Real>(_real).subspan(SpanOffset<Flat, true>(),
                                                SpanExtent(RealComponentStride));
  }

  // Slots are numbered across both buffers, complex first, so a real
  // component's index within its own buffer is its slot less the complex
  // count.
  template <Int Flat, bool InRealBuffer>
  std::size_t SpanOffset() const {
    constexpr auto slot = SlotOfFlat(Orbits.representative[Flat]);
    static_assert(slot >= 0);
    constexpr auto index = InRealBuffer ? slot - ComplexComponents : slot;
    static_assert(index >= 0);
    if constexpr (IsComponentMajor) {
      return static_cast<std::size_t>(index) *
             static_cast<std::size_t>(_grid.FieldSize());
    } else {
      return static_cast<std::size_t>(index);
    }
  }

  // The elements a strided component is spread over: the last sample's offset
  // plus one, which is less than the whole buffer by the slots that follow.
  std::size_t SpanExtent(Int stride) const {
    return static_cast<std::size_t>((_grid.FieldSize() - 1) * stride + 1);
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
