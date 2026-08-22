// 14 -- The generalized Legendre functions
//
// Underneath every transform is a table of d^l_{Nm}(theta). This is how to
// get at them directly, what the stored value actually is, and how to check
// it against the two relations that pin it.
//
// It also writes the values out, so they can be plotted. If gnuplot is on the
// path the example draws the figure itself; otherwise the data file is there
// to be plotted however you like.

#include <GSHTrans/All>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Int = std::ptrdiff_t;

  constexpr auto lMax = Int{8};
  constexpr auto nMax = Int{2};

  // A table over a set of colatitudes. The template arguments say: all orders
  // m, all upper indices n, and many angles rather than one. Building it runs
  // the recursion once for each (n, theta).
  const auto nTheta = Int{257};
  auto theta = std::vector<Real>(nTheta);
  for (Int i = 0; i < nTheta; i++) {
    theta[i] = std::numbers::pi * i / (nTheta - 1);
  }
  auto d = Wigner<Real, All, All, Multiple>(lMax, lMax, nMax, theta);

  // What is stored is *not* d^l_{Nm} but the orthonormalised scalar harmonic
  //
  //     X^N_{lm}(theta) = sqrt((2l+1)/4pi) d^l_{Nm}(theta),
  //
  // which is Dahlen & Tromp (C.117). The upper index is the *first* subscript
  // of d: d^l_{mN} is a different function, differing by (-1)^{N-m}.
  const auto view = d[1, 128];        // upper index 1, the colatitude at pi/2
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "at theta = pi/2, upper index N = 1:\n";
  for (auto m : view[3].Orders()) {
    std::cout << "  X^1_{3," << std::setw(2) << m << "} = " << std::setw(10)
              << (view[3][m]) << "\n";
  }
  std::cout << "\n";

  // Two relations worth checking by hand, because between them they pin the
  // whole table.
  //
  // The addition theorem, D&T (C.127): sum_m P^N_{lm} P^{N'}_{lm} = delta.
  // In terms of what is stored the right-hand side carries the square of the
  // normalisation, so it is delta * (2l+1)/4pi.
  const auto l = Int{5};
  const auto iTheta = Int{97};
  const auto expected = (2 * l + 1) / (4 * std::numbers::pi);
  for (auto n : {Int{0}, Int{1}, Int{2}}) {
    auto sum = Real{0};
    auto cross = Real{0};
    for (auto m : d[n, iTheta][l].Orders()) {
      sum += d[n, iTheta][l][m] * d[n, iTheta][l][m];
      cross += d[n, iTheta][l][m] * d[0, iTheta][l][m];
    }
    std::cout << "addition theorem, N = N' = " << n << ": " << sum
              << "   expected " << expected << "\n";
    if (n != 0) {
      std::cout << "                  N = " << n << ", N' = 0: " << cross
                << "   expected 0\n";
    }
  }
  std::cout << "\n";

  // The involution P^{-N}_{l,-m} = (-1)^{m+N} P^N_{lm}, D&T (C.118). This is
  // one of the two symmetries that could halve the stored table; neither is
  // used, for reasons the reference document gives.
  auto worst = Real{0};
  for (auto n : d.UpperIndices()) {
    for (auto ll : d.Degrees(n)) {
      for (auto m : d[n, iTheta][ll].Orders()) {
        const auto lhs = d[-n, iTheta][ll][-m];
        const auto rhs = MinusOneToPower(m + n) * d[n, iTheta][ll][m];
        worst = std::max(worst, std::abs(lhs - rhs));
      }
    }
  }
  std::cout << "involution P^{-N}_{l,-m} = (-1)^{m+N} P^N_{lm}, worst error "
            << std::scientific << worst << "\n\n";

  // Write the values out. One column per upper index, at fixed degree and
  // order, over the whole range of colatitude.
  const auto degree = Int{6};
  const auto order = Int{2};
  const auto path = std::string("wigner-l6-m2.dat");
  {
    auto out = std::ofstream(path);
    out << "# theta";
    for (auto n : d.UpperIndices()) out << "  X^" << n << "_{6,2}";
    out << "\n" << std::fixed << std::setprecision(10);
    for (auto i : d.AngleIndices()) {
      out << theta[i];
      for (auto n : d.UpperIndices()) out << " " << d[n, i][degree][order];
      out << "\n";
    }
  }
  std::cout << "wrote " << path << ": " << nTheta << " colatitudes, "
            << d.NumberOfUpperIndices() << " upper indices\n";

  // Plot it if gnuplot is there. The functions are real, oscillate more with
  // degree, and are pushed away from the poles as |m| grows -- which is the
  // fall-off like sin^|m|(theta) that a polar truncation would exploit.
  if (std::system("command -v gnuplot > /dev/null 2>&1") == 0) {
    auto script = std::ofstream("wigner.gp");
    script << "set terminal pngcairo size 900,540\n"
              "set output 'wigner-l6-m2.png'\n"
              "set title 'X^N_{6,2}(theta), the stored Wigner values'\n"
              "set xlabel 'theta'; set ylabel 'X^N_{6,2}'\n"
              "set xrange [0:pi]; set grid; set key outside\n"
              "set xtics ('0' 0, 'pi/2' pi/2, 'pi' pi)\n"
              "plot ";
    Int column = 2;
    for (auto n : d.UpperIndices()) {
      script << (column > 2 ? ", " : "") << "'" << path << "' using 1:"
             << column << " with lines lw 2 title 'N = " << n << "'";
      column++;
    }
    script << "\n";
    script.close();
    if (std::system("gnuplot wigner.gp") == 0) {
      std::cout << "wrote wigner-l6-m2.png\n";
    }
  } else {
    std::cout << "gnuplot not found; the data file is there to plot\n";
  }
}
