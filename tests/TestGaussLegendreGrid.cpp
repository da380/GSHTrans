#include <gtest/gtest.h>

#include "CheckCoeff2Coeff.h"

TEST(GaussLegendreGrid, Coeff2CoeffDoubleR2C) {
  using Scalar = double;
  bool result = Coeff2Coeff<Scalar, All, All>();
  EXPECT_FALSE(result);
}

TEST(GaussLegendreGrid, Coeff2CoeffLongDoubleR2C) {
  using Scalar = long double;
  bool result = Coeff2Coeff<Scalar, All, All>();
  EXPECT_FALSE(result);
}

TEST(GaussLegendreGrid, Coeff2CoeffDoubleC2C) {
  using Scalar = std::complex<double>;
  bool result = Coeff2Coeff<Scalar, All, All>();
  EXPECT_FALSE(result);
}

TEST(GaussLegendreGrid, Coeff2CoeffLongDoubleC2C) {
  using Scalar = std::complex<long double>;
  bool result = Coeff2Coeff<Scalar, All, All>();
  EXPECT_FALSE(result);
}
