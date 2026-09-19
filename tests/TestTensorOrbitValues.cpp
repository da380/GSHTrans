#include <gtest/gtest.h>

#include <GSHTrans/All>
#include <array>
#include <complex>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "TestRandom.h"

// The values a tensor hands back for its derived components.
//
// A tensor stores one representative per orbit and derives the rest, and the
// derivation is written once, in Orbits.h. What is checked here is not that
// helper against a second copy of itself but the *defining relations*, which
// are independent of how any component is obtained:
//
//   permutation  T^{pi(alpha)} = s T^{alpha}             for every generator
//   reality      T^{-alpha}    = (-1)^N conj(T^{alpha})  for a real tensor
//
// and their spectral forms,
//
//   C^{pi(alpha)}_{lm} = s C^{alpha}_{lm}
//   C^{-alpha}_{lm}    = (-1)^m conj(C^{alpha}_{l,-m}),
//
// the second because conj of a field of upper index N has coefficients
// (-1)^{m+N} conj(f_{l,-m}). Every component of every type is read, so a sign
// that is wrong on one member of one orbit -- which is what this file was
// written after -- has nowhere to hide. Rank 2 cannot show such an error: no
// pinned member there is reached through a conjugating step. Rank 3 and the
// elastic tensor can, which is why they are here and with values in them.

using namespace GSHTrans;

namespace {

using Int = std::ptrdiff_t;
using Real = double;
using Complex = std::complex<Real>;
using Grid = GaussLegendreGrid<Real, All, All>;

constexpr Int lMax = 4;
auto TestGrid() { return Grid(lMax, 4, FFTWpp::Estimate); }

// A Riemann-like symmetry: antisymmetric in each pair, symmetric under their
// exchange. Here because it is a GeneratedBy rather than a named symmetry, and
// because its orbits mix signs with the reality relation.
using Riemann = GeneratedBy<4, SlotPermutation<4>{{1, 0, 2, 3}, -1},
                            SlotPermutation<4>{{0, 1, 3, 2}, -1},
                            SlotPermutation<4>{{2, 3, 0, 1}, 1}>;

//--------------------------------------------------------------------------//
//                 Reading a component chosen by its flat index              //
//--------------------------------------------------------------------------//

template <typename T, auto Slots, std::size_t... I>
Complex ValueOf(const T& t, Int iTheta, Int iPhi, std::index_sequence<I...>) {
  if constexpr (T::template Represents<Slots[I]...>) {
    return Complex{t.template Component<Slots[I]...>()[iTheta, iPhi]};
  } else {
    return Complex{};
  }
}

// Every component at one point, by flat index. Read through a const reference,
// which is the accessor that derives.
template <typename T>
auto ValuesAt(const T& t, Int iTheta, Int iPhi) {
  auto values = std::array<Complex, T::Components>{};
  [&]<std::size_t... F>(std::index_sequence<F...>) {
    ((values[F] = ValueOf<T, T::Index::FromFlat(Int{F}).Slots()>(
          t, iTheta, iPhi, std::make_index_sequence<T::Rank>{})),
     ...);
  }(std::make_index_sequence<T::Components>{});
  return values;
}

template <typename E, auto Slots, std::size_t... I>
Complex CoefficientOf(const E& e, Int l, Int m, std::index_sequence<I...>) {
  return e.template Coefficient<Slots[I]...>(l, m);
}

template <typename T, typename E>
auto CoefficientsAt(const E& e, Int l, Int m) {
  auto values = std::array<Complex, T::Components>{};
  [&]<std::size_t... F>(std::index_sequence<F...>) {
    ((values[F] = CoefficientOf<E, T::Index::FromFlat(Int{F}).Slots()>(
          e, l, m, std::make_index_sequence<T::Rank>{})),
     ...);
  }(std::make_index_sequence<T::Components>{});
  return values;
}

template <typename T>
void FillAtRandom(T& t, GSHTransTest::Generator& gen) {
  auto dist = std::uniform_real_distribution<Real>(-1, 1);
  for (auto& x : t.Data()) x = Complex{dist(gen), dist(gen)};
  for (auto& x : t.RealData()) x = dist(gen);
}

template <typename T>
constexpr bool IsReal = std::same_as<typename T::Reality, RealTensor>;

//--------------------------------------------------------------------------//
//                               The type list                               //
//--------------------------------------------------------------------------//

template <typename T>
class TensorOrbitValues : public ::testing::Test {};

using TensorTypes = ::testing::Types<
    TensorField<1, NoSymmetry<1>, RealTensor, Grid>,
    TensorField<2, NoSymmetry<2>, RealTensor, Grid>,
    TensorField<2, Symmetric<2>, RealTensor, Grid>,
    TensorField<2, Antisymmetric<2>, RealTensor, Grid>,
    TensorField<3, Symmetric<3>, RealTensor, Grid>,
    TensorField<3, Antisymmetric<3>, RealTensor, Grid>,
    TensorField<4, ElasticSymmetry, RealTensor, Grid>,
    TensorField<4, ElasticSymmetry, RealTensor, Grid, PointMajor>,
    TensorField<4, ElasticSymmetry, ComplexTensor, Grid>,
    TensorField<4, Symmetric<4>, RealTensor, Grid>,
    TensorField<4, Riemann, RealTensor, Grid>,
    TensorField<2, Symmetric<2>, RealTensor, Grid, ComponentMajor,
                TangentialSlots>>;
TYPED_TEST_SUITE(TensorOrbitValues, TensorTypes);

}  // namespace

//--------------------------------------------------------------------------//
//                                  Spatial                                  //
//--------------------------------------------------------------------------//

TYPED_TEST(TensorOrbitValues, EveryComponentObeysEveryRelation) {
  using T = TypeParam;
  const auto seed = GSHTransTest::TestSeed();
  auto gen = GSHTransTest::MakeGenerator(seed);

  auto grid = TestGrid();
  auto t = T(grid);
  FillAtRandom(t, gen);

  const auto values = ValuesAt(std::as_const(t), 1, 2);

  for (auto flat = Int{0}; flat < T::Components; flat++) {
    const auto index = T::Index::FromFlat(flat);
    for (const auto& generator : T::Symmetry::Generators()) {
      const auto to = index.Permuted(generator.image).Flat();
      EXPECT_EQ(values[to], static_cast<Real>(generator.sign) * values[flat])
          << "permutation, from flat " << flat << " to " << to << ", seed "
          << seed;
    }
    if constexpr (IsReal<T>) {
      const auto to = index.Negated().Flat();
      const auto sign = static_cast<Real>(MinusOneToPower(index.UpperIndex()));
      EXPECT_EQ(values[to], sign * std::conj(values[flat]))
          << "reality, from flat " << flat << " to " << to << ", seed " << seed;
    }
  }
}

//--------------------------------------------------------------------------//
//                                  Spectral                                 //
//--------------------------------------------------------------------------//

TYPED_TEST(TensorOrbitValues, EveryCoefficientObeysEveryRelation) {
  using T = TypeParam;
  const auto seed = GSHTransTest::TestSeed();
  auto gen = GSHTransTest::MakeGenerator(seed);

  auto grid = TestGrid();
  auto t = T(grid);
  FillAtRandom(t, gen);
  const auto e = Expand(std::as_const(t), lMax);
  using E = std::remove_cvref_t<decltype(e)>;

  for (auto l = Int{0}; l <= lMax; l++) {
    for (auto m = -l; m <= l; m++) {
      const auto here = CoefficientsAt<T, E>(e, l, m);
      const auto mirror = CoefficientsAt<T, E>(e, l, -m);
      for (auto flat = Int{0}; flat < T::Components; flat++) {
        const auto index = T::Index::FromFlat(flat);
        for (const auto& generator : T::Symmetry::Generators()) {
          const auto to = index.Permuted(generator.image).Flat();
          EXPECT_EQ(here[to], static_cast<Real>(generator.sign) * here[flat])
              << "permutation, from flat " << flat << " to " << to << " at ("
              << l << ", " << m << "), seed " << seed;
        }
        if constexpr (IsReal<T>) {
          const auto to = index.Negated().Flat();
          const auto sign = static_cast<Real>(MinusOneToPower(m));
          EXPECT_EQ(here[to], sign * std::conj(mirror[flat]))
              << "reality, from flat " << flat << " to " << to << " at (" << l
              << ", " << m << "), seed " << seed;
        }
      }
    }
  }
}

// The two sides against each other. Each relation above is absolute, but
// nothing in them ties a spatial component to the spectral one of the same
// name. This does: the coefficient the expansion reports for a component must
// be the transform of the field the tensor reports for it. The transform is
// linear and the relations are pointwise, so this is exact to rounding however
// rough the data.
namespace {

template <typename T, typename E, auto Slots, std::size_t... I>
void ExpectCoefficientsAreTheTransform(const T& t, const E& e,
                                       std::index_sequence<I...>) {
  if constexpr (T::template Represents<Slots[I]...>) {
    constexpr auto N = T::template UpperIndexOf<Slots[I]...>;
    auto component = t.template Component<Slots[I]...>();
    auto field = SpinField<N, Grid>(t.Grid());
    for (auto [iTheta, iPhi] : t.Grid().PointIndices()) {
      field[iTheta, iPhi] = Complex{component[iTheta, iPhi]};
    }
    const auto direct = Expand(field, lMax);
    for (auto l = (N < 0 ? -N : N); l <= lMax; l++) {
      for (auto m = -l; m <= l; m++) {
        const auto got = e.template Coefficient<Slots[I]...>(l, m);
        EXPECT_NEAR(std::abs(got - direct[l, m]), 0, 1e-13)
            << "component of upper index " << N << " at (" << l << ", " << m
            << ")";
      }
    }
  }
}

}  // namespace

TYPED_TEST(TensorOrbitValues, ACoefficientIsTheTransformOfItsComponent) {
  using T = TypeParam;
  auto gen = GSHTransTest::MakeGenerator(GSHTransTest::TestSeed());

  auto grid = TestGrid();
  auto t = T(grid);
  FillAtRandom(t, gen);
  const auto e = Expand(std::as_const(t), lMax);
  using E = std::remove_cvref_t<decltype(e)>;

  [&]<std::size_t... F>(std::index_sequence<F...>) {
    (ExpectCoefficientsAreTheTransform<T, E,
                                       T::Index::FromFlat(Int{F}).Slots()>(
         std::as_const(t), e, std::make_index_sequence<T::Rank>{}),
     ...);
  }(std::make_index_sequence<T::Components>{});
}

//--------------------------------------------------------------------------//
//                         The layered type against the flat                 //
//--------------------------------------------------------------------------//

// The layered tensor defers its combinatorics to the flat one, so the two must
// agree component for component, in value and in kind. The reads go through a
// const reference because that is the accessor that derives, and the
// comparison is made on a component held in a named variable because that is
// how it is used -- an expression that refers to a view which has since died
// reads as garbage here and as a stack-use-after-return under the address
// sanitiser.
namespace {

template <typename L, auto Slots, std::size_t... I>
void FillStack(L& layered, GSHTransTest::Generator& gen,
               std::index_sequence<I...>) {
  if constexpr (L::template Writable<Slots[I]...>) {
    auto dist = std::uniform_real_distribution<Real>(-1, 1);
    using Scalar = std::ranges::range_value_t<
        decltype(layered.template ComponentStack<Slots[I]...>().Data())>;
    for (auto& x : layered.template ComponentStack<Slots[I]...>().Data()) {
      if constexpr (std::same_as<Scalar, Real>) {
        x = dist(gen);
      } else {
        x = Complex{dist(gen), dist(gen)};
      }
    }
  }
}

template <typename L, auto Slots, std::size_t... I>
void CopySlice(L& layered, typename L::Flat& flat, Int i,
               std::index_sequence<I...>) {
  if constexpr (L::template Writable<Slots[I]...>) {
    auto from = layered.template Component<Slots[I]...>(i);
    auto to = flat.template Component<Slots[I]...>();
    for (auto [iTheta, iPhi] : flat.Grid().PointIndices()) {
      to[iTheta, iPhi] = from[iTheta, iPhi];
    }
  }
}

template <typename L, auto Slots, std::size_t... I>
void ExpectSliceMatches(const L& layered, const typename L::Flat& flat, Int i,
                        std::index_sequence<I...>) {
  if constexpr (L::template Represents<Slots[I]...>) {
    auto got = layered.template Component<Slots[I]...>(i);
    auto wanted = flat.template Component<Slots[I]...>();
    EXPECT_TRUE((std::same_as<typename decltype(got)::Value,
                              typename decltype(wanted)::Value>))
        << "the two types disagree on whether this component is real-valued";
    for (auto [iTheta, iPhi] : flat.Grid().PointIndices()) {
      EXPECT_EQ((Complex{got[iTheta, iPhi]}), (Complex{wanted[iTheta, iPhi]}));
    }
  }
}

template <typename L, auto Slots, std::size_t... I>
Complex LayeredCoefficientOf(const L& e, Int i, Int l, Int m,
                             std::index_sequence<I...>) {
  return e.template Coefficient<Slots[I]...>(i, l, m);
}

template <typename T>
class LayeredOrbitValues : public ::testing::Test {};

using LayeredTypes =
    ::testing::Types<LayeredTensorField<1, NoSymmetry<1>, RealTensor, Grid>,
                     LayeredTensorField<2, Symmetric<2>, RealTensor, Grid>,
                     LayeredTensorField<2, Antisymmetric<2>, RealTensor, Grid>,
                     LayeredTensorField<3, Antisymmetric<3>, RealTensor, Grid>,
                     LayeredTensorField<4, ElasticSymmetry, RealTensor, Grid>>;
TYPED_TEST_SUITE(LayeredOrbitValues, LayeredTypes);

}  // namespace

TYPED_TEST(LayeredOrbitValues, EverySliceIsTheFlatTensorsComponent) {
  using L = TypeParam;
  using F = typename L::Flat;
  auto gen = GSHTransTest::MakeGenerator(GSHTransTest::TestSeed());

  auto grid = TestGrid();
  auto radial = RadialGrid<Real>(std::vector<Real>{0.5, 0.75, 1.0});
  auto layered = L(radial, grid);

  constexpr auto all = std::make_index_sequence<F::Components>{};
  constexpr auto rank = std::make_index_sequence<F::Rank>{};

  [&]<std::size_t... K>(std::index_sequence<K...>) {
    (FillStack<L, F::Index::FromFlat(Int{K}).Slots()>(layered, gen, rank), ...);
  }(all);

  for (auto i : layered.RadiusIndices()) {
    auto flat = F(grid);
    [&]<std::size_t... K>(std::index_sequence<K...>) {
      (CopySlice<L, F::Index::FromFlat(Int{K}).Slots()>(layered, flat, i, rank),
       ...);
    }(all);
    [&]<std::size_t... K>(std::index_sequence<K...>) {
      (ExpectSliceMatches<L, F::Index::FromFlat(Int{K}).Slots()>(
           std::as_const(layered), std::as_const(flat), i, rank),
       ...);
    }(all);
  }
}

TYPED_TEST(LayeredOrbitValues, EveryCoefficientObeysEveryRelation) {
  using L = TypeParam;
  using F = typename L::Flat;
  const auto seed = GSHTransTest::TestSeed();
  auto gen = GSHTransTest::MakeGenerator(seed);

  auto grid = TestGrid();
  auto radial = RadialGrid<Real>(std::vector<Real>{0.5, 1.0});
  auto layered = L(radial, grid);

  constexpr auto all = std::make_index_sequence<F::Components>{};
  constexpr auto rank = std::make_index_sequence<F::Rank>{};
  [&]<std::size_t... K>(std::index_sequence<K...>) {
    (FillStack<L, F::Index::FromFlat(Int{K}).Slots()>(layered, gen, rank), ...);
  }(all);

  const auto e = Expand(std::as_const(layered), lMax);
  using E = std::remove_cvref_t<decltype(e)>;

  const auto at = [&](Int i, Int l, Int m) {
    auto values = std::array<Complex, F::Components>{};
    [&]<std::size_t... K>(std::index_sequence<K...>) {
      ((values[K] = LayeredCoefficientOf<E, F::Index::FromFlat(Int{K}).Slots()>(
            e, i, l, m, rank)),
       ...);
    }(all);
    return values;
  };

  for (auto i : layered.RadiusIndices()) {
    for (auto l = Int{0}; l <= lMax; l++) {
      for (auto m = -l; m <= l; m++) {
        const auto here = at(i, l, m);
        const auto mirror = at(i, l, -m);
        for (auto flat = Int{0}; flat < F::Components; flat++) {
          const auto index = F::Index::FromFlat(flat);
          for (const auto& generator : F::Symmetry::Generators()) {
            const auto to = index.Permuted(generator.image).Flat();
            EXPECT_EQ(here[to], static_cast<Real>(generator.sign) * here[flat])
                << "permutation, flat " << flat << " to " << to << ", seed "
                << seed;
          }
          const auto to = index.Negated().Flat();
          const auto sign = static_cast<Real>(MinusOneToPower(m));
          EXPECT_EQ(here[to], sign * std::conj(mirror[flat]))
              << "reality, flat " << flat << " to " << to << ", seed " << seed;
        }
      }
    }
  }
}

//--------------------------------------------------------------------------//
//                      The composite the layer exists for                   //
//--------------------------------------------------------------------------//

// A real elastic tensor applied to a real strain, against the double sum
// written out:
//
//   s^{ij} = sum_{a,b} (-1)^{a+b} c^{ijab} e^{-a,-b},
//
// which is what contracting slot k with m and then l with n means under the
// canonical metric. The oracle reads both tensors one component at a time and
// does its own arithmetic, so it shares nothing with the expression layer but
// the component accessor -- and that is what the tests above pin.
//
// This is here rather than beside the algebra tests because it is the case
// that went wrong: most of the components this sum reads from a real elastic
// tensor are derived ones, and the ones on real orbits came back negated.
TEST(RealElasticTensor, AppliedToAStrainIsTheDoubleSum) {
  using C = TensorField<4, ElasticSymmetry, RealTensor, Grid>;
  using E = TensorField<2, Symmetric<2>, RealTensor, Grid>;
  auto gen = GSHTransTest::MakeGenerator(GSHTransTest::TestSeed());

  auto grid = TestGrid();
  auto c = C(grid);
  auto e = E(grid);
  FillAtRandom(c, gen);
  FillAtRandom(e, gen);

  const auto& elastic = c;
  const auto& strain = e;
  auto stress = Contract<2, 3>(Contract<2, 4>(TensorProduct(elastic, strain)));
  static_assert(decltype(stress)::Rank == 2);

  constexpr auto letters = std::array<Int, 3>{-1, 0, 1};
  const auto iTheta = Int{2};
  const auto iPhi = Int{3};
  const auto cAt = ValuesAt(elastic, iTheta, iPhi);
  const auto eAt = ValuesAt(strain, iTheta, iPhi);

  const auto expected = [&](Int i, Int j) {
    auto sum = Complex{};
    for (auto a : letters) {
      for (auto b : letters) {
        const auto cFlat = C::Index(std::array<Int, 4>{i, j, a, b}).Flat();
        const auto eFlat = E::Index(std::array<Int, 2>{-a, -b}).Flat();
        sum +=
            static_cast<Real>(MinusOneToPower(a + b)) * cAt[cFlat] * eAt[eFlat];
      }
    }
    return sum;
  };

  const auto expectNear = [&](Complex got, Complex wanted) {
    EXPECT_NEAR(std::abs(got - wanted), 0, 1e-14);
  };
  expectNear((stress.Component<-1, -1>()[iTheta, iPhi]), expected(-1, -1));
  expectNear((stress.Component<-1, 0>()[iTheta, iPhi]), expected(-1, 0));
  expectNear((stress.Component<-1, 1>()[iTheta, iPhi]), expected(-1, 1));
  expectNear((stress.Component<0, -1>()[iTheta, iPhi]), expected(0, -1));
  expectNear((stress.Component<0, 0>()[iTheta, iPhi]), expected(0, 0));
  expectNear((stress.Component<0, 1>()[iTheta, iPhi]), expected(0, 1));
  expectNear((stress.Component<1, -1>()[iTheta, iPhi]), expected(1, -1));
  expectNear((stress.Component<1, 0>()[iTheta, iPhi]), expected(1, 0));
  expectNear((stress.Component<1, 1>()[iTheta, iPhi]), expected(1, 1));

  // A real tensor applied to a real tensor is real, and the stress of a
  // symmetric strain through an elastic tensor is symmetric. Neither is built
  // into the expression, so both are checks on the values.
  const auto at = [&](Int i, Int j) { return expected(i, j); };
  for (auto i : letters) {
    for (auto j : letters) {
      expectNear(at(i, j), at(j, i));
      expectNear(at(-i, -j), static_cast<Real>(MinusOneToPower(i + j)) *
                                 std::conj(at(i, j)));
    }
  }
}
