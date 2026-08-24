// 22 -- Wigner 3-j symbols, and the check that says when to trust them
//
// The 3-j symbols are the coupling coefficients of three angular momenta, and
// what a spherical-harmonic library wants them for is Gaunt integrals and
// mode coupling. The library offers three shapes:
//
//   Wigner3jSymbol(l1,l2,l3, m1,m2,m3)  one symbol, for convenience
//   Wigner3jMatrix(l1,l2,l3)            the whole (m1, m3) plane at fixed
//                                       degrees, with m2 fixed by the rule
//   Wigner3jStack(l1,l3)                one such matrix per middle degree l2
//
// The table is the primitive: a symbol on its own costs a whole table, so ask
// for the table if you want more than one.
//
// **The important thing in this example is the last section.** The recursion
// behind these values is accurate over a wide range and breaks down near
// *stretched* triangles -- where one degree approaches the sum of the other
// two -- at high degree. It breaks down silently: the numbers stay finite and
// the table keeps its shape while the values become meaningless. So every
// table checks itself against the completeness relation and refuses rather
// than returning nonsense.
//
// See docs/3j-plan.md.

#include <GSHTrans/Core>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main() {
  using namespace GSHTrans;

  std::cout << std::scientific << std::setprecision(3);

  //------------------------------------------------------------------------//
  //                        A table, and what is in it                       //
  //------------------------------------------------------------------------//

  const auto table = Wigner3jMatrix<double>(3, 4, 5);
  std::cout << "A table for degrees (3, 4, 5) holds " << table.size()
            << " values,\n"
            << "over m1 in [-3, 3] and m3 in [-5, 5], with m2 = -(m1 + m3).\n";

  std::cout << "\n  (3 4 5; 1 -3 2) = " << table(1, -3, 2) << "\n"
            << "  the same, by (m1, m3) = " << table(1, 2) << "\n"
            << "  and the free function = "
            << Wigner3jSymbol<double>(3, 4, 5, 1, -3, 2) << "\n";

  // Queries that break a selection rule are zero rather than an error: that
  // is the mathematically consistent value, and it keeps a coupling sum from
  // needing a guard at every term.
  std::cout << "\nOutside the selection rules the symbol is zero, not an "
               "error:\n"
            << "  orders not summing to zero = " << table(1, 1, 1) << "\n"
            << "  an order beyond its degree = " << table(7, 0) << "\n";

  //------------------------------------------------------------------------//
  //                       Two closed forms, as a check                      //
  //------------------------------------------------------------------------//

  // (j j 0; m -m 0) = (-1)^{j-m} / sqrt(2j+1), which the library reproduces
  // exactly.
  {
    constexpr auto j = 7;
    const auto zeroDegree = Wigner3jMatrix<double>(j, j, 0);
    auto worst = 0.0;
    for (auto m : zeroDegree.M1Axis()) {
      const auto want =
          ((j - m) % 2 == 0 ? 1.0 : -1.0) / std::sqrt(2.0 * j + 1);
      worst = std::max(worst, std::abs(zeroDegree(m, 0) - want));
    }
    std::cout << "\n(j j 0; m -m 0) against its closed form, j = " << j
              << ": " << worst << "\n";
    if (!(worst < 1e-14)) return 1;
  }

  // The fully stretched symbol, where every factorial cancels:
  //   (l1 l2 l1+l2; l1 l2 -(l1+l2)) = 1 / sqrt(2 l3 + 1)
  {
    constexpr auto l1 = 5;
    constexpr auto l2 = 8;
    constexpr auto l3 = l1 + l2;
    const auto stretched = Wigner3jMatrix<double>(l1, l2, l3);
    const auto want = 1 / std::sqrt(2.0 * l3 + 1);
    const auto error = std::abs(stretched(l1, l2, -l3) - want);
    std::cout << "the fully stretched symbol against 1/sqrt(2 l3 + 1): "
              << error << "\n";
    if (!(error < 1e-14)) return 1;
  }

  //------------------------------------------------------------------------//
  //                    A stack, over the middle degree                      //
  //------------------------------------------------------------------------//

  {
    const auto stack = Wigner3jStack<double>(4, 6);
    std::cout << "\nA stack for (l1, l3) = (4, 6) covers l2 in ["
              << stack.L2Axis().Min() << ", " << stack.L2Axis().Max()
              << "],\nwhich is the triangle range, as " << stack.size()
              << " tables.\n";
  }

  //------------------------------------------------------------------------//
  //             The completeness relation, and the refusal                  //
  //------------------------------------------------------------------------//

  // Summing the squares over the whole plane gives one, exactly, for any
  // triple satisfying the triangle rule. It needs no reference to compare
  // against, which is what lets the library check itself.
  {
    auto sum = 0.0;
    for (auto value : table) sum += value * value;
    std::cout << "\nThe completeness relation on (3, 4, 5): sum of squares = "
              << sum << ", departure from one = " << std::abs(sum - 1)
              << "\n";
  }

  // Near a stretched triangle at high degree the recursion loses the values
  // altogether. Without the self-check this would hand back a table of
  // numbers of order 1e112, finite and correctly shaped.
  try {
    const auto doomed = Wigner3jMatrix<double>(40, 40, 80);
    std::cout << "\n(40, 40, 80) should have been refused\n";
    (void)doomed;
    return 1;
  } catch (const std::runtime_error& error) {
    std::cout << "\n(40, 40, 80) is refused, and says why:\n  "
              << error.what() << "\n";
  }

  std::cout << "\nWhere it is accurate, it is accurate to rounding; where it "
               "is not,\nit says so. The working range is wide -- (l,l,l) "
               "holds to l = 64 and\nbeyond -- and it narrows as a triangle "
               "approaches stretched.\n";

  return 0;
}
