#include <gtest/gtest.h>

#include <array>

#include "CheckAdditionTheorem.h"
#include "CheckLegendre.h"

namespace {

void CheckSingleUpperIndexAccess(std::ptrdiff_t n) {
  using namespace GSHTrans;

  constexpr std::ptrdiff_t lMax = 5;
  constexpr std::ptrdiff_t mMax = 3;

  auto singleAngle =
      Wigner<double, FourPi, All, Single, Single>(lMax, mMax, n, 0.7);
  auto explicitSingleAngle = singleAngle[n, 0];
  for (auto l : singleAngle.Degrees()) {
    for (auto m : explicitSingleAngle[l].Orders()) {
      EXPECT_DOUBLE_EQ(singleAngle[l][m], explicitSingleAngle[l][m]);
    }
  }

  const auto angles = std::array{0.2, 0.7, 1.3};
  auto multipleAngles =
      Wigner<double, FourPi, All, Single, Multiple>(lMax, mMax, n, angles);
  for (auto iTheta : multipleAngles.AngleIndices()) {
    auto implicitView = multipleAngles[iTheta];
    auto explicitView = multipleAngles[n, iTheta];
    for (auto l : explicitView.Degrees()) {
      for (auto m : explicitView[l].Orders()) {
        EXPECT_DOUBLE_EQ(implicitView[l][m], explicitView[l][m]);
      }
    }
  }
}

}  // namespace

// Compare values for n = 0 to the std library function.
TEST(Wigner, CheckLegendreDouble) {
  int i = CheckLegendre<double>();
  EXPECT_EQ(i, 0);
}

TEST(Wigner, CheckLegendreLongDouble) {
  int i = CheckLegendre<long double>();
  EXPECT_EQ(i, 0);
}

// Check the addition theorem is satisfied.
TEST(Wigner, CheckAdditionTheoremDouble) {
  int i = CheckAdditionTheorem<double>();
  EXPECT_EQ(i, 0);
}

TEST(Wigner, CheckAdditionTheoremLongDouble) {
  int i = CheckAdditionTheorem<long double>();
  EXPECT_EQ(i, 0);
}

TEST(Wigner, SinglePositiveUpperIndexAccess) {
  CheckSingleUpperIndexAccess(2);
}

TEST(Wigner, SingleNegativeUpperIndexAccess) {
  CheckSingleUpperIndexAccess(-2);
}

TEST(Wigner, SingleMaximumUpperIndexAccess) {
  CheckSingleUpperIndexAccess(5);
  CheckSingleUpperIndexAccess(-5);
}
