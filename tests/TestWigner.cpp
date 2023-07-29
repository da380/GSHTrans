#include <gtest/gtest.h>

#include "CheckAdditionTheorem.h"
#include "CheckUpperIndexZero.h"

// Compare values for n = 0 to the std library function.
TEST(Wigner, CheckUpperIndexZeroDouble) {
  int i = CheckUpperIndexZero<double>();
  EXPECT_EQ(i, 0);
}

TEST(Wigner, CheckUpperIndexZeroLongDouble) {
  int i = CheckUpperIndexZero<long double>();
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
