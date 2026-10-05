// 14 -- The generalized Legendre functions
//
// What this shows. Underneath every transform is a table of the Wigner
// d-functions d^l_{Nm}(theta), the colatitude part of the generalised
// spherical harmonics:
//
//     Y^N_{lm}(theta, phi) = sqrt((2l+1)/4pi) d^l_{Nm}(theta) exp(i m phi).
//
// The forward transform is an FFT in longitude followed by a quadrature sum
// over colatitude against these values (docs/gshtrans-reference.tex, "The
// transform"). This example builds a table directly with the Wigner class,
// says exactly what is stored, checks it against two of the relations that
// pin the convention, and writes values out for plotting.
//
// Read first. 01 (the grid) and 12 (expansions). The grid of the other
// examples holds a table of exactly this kind at its own Gauss-Legendre
// colatitudes; here the colatitudes are chosen freely.
//
// Introduced. Wigner<Real, MRange, NRange, AngleRange>, indexed
// d[n, iTheta][l][m]; UpperIndices(), Degrees(n), AngleIndices(); the stored
// normalisation X^N_{lm}; the addition theorem and the (N, m) -> (-N, -m)
// symmetry; MinusOneToPower.
//
// The output. The values X^1_{3,m}(pi/2) for every m; the addition theorem
// holding for N = N' and vanishing for N != N'; the symmetry holding to
// rounding over the whole table; and a data file wigner-l6-m2.dat (plus a
// PNG if gnuplot is on the path, otherwise the data file is there to be
// plotted however you like).

#include <GSHTrans/GSHTrans.hpp>
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

  // A table over 257 equally spaced colatitudes from 0 to pi inclusive, so
  // index 128 is the equator. The template arguments say: all orders m, all
  // upper indices n (from -nMax to nMax), and many angles rather than one
  // (Single would hold one). The constructor takes (lMax, mMax, nMax, angles).
  // Building it runs a recursion in degree once for each (n, theta).
  const auto nTheta = Int{257};
  auto theta = std::vector<Real>(nTheta);
  for (Int i = 0; i < nTheta; i++) {
    theta[i] = std::numbers::pi * i / (nTheta - 1);
  }
  auto d = Wigner<Real, All, All, Multiple>(lMax, lMax, nMax, theta);

  // What is stored is *not* d^l_{Nm} but the orthonormalised colatitude part
  // of the harmonic,
  //
  //     X^N_{lm}(theta) = sqrt((2l+1)/4pi) d^l_{Nm}(theta),
  //
  // which is Dahlen & Tromp (C.117), their P^N_{lm} being d^l_{Nm}. The upper
  // index is the *first* subscript of d: d^l_{mN} is a different function,
  // differing by (-1)^{N-m}, so swapping them is a sign error whenever N - m
  // is odd. Woodhouse and Phinney & Burridge use the same index order but no
  // sqrt((2l+1)/4pi); the reference's "Translation table" compares the three.
  //
  // d[n, iTheta] is the block of every (l, m) at one upper index and one
  // colatitude; [l] then selects a degree and [m] an order.
  const auto view = d[1, 128];  // upper index 1, the colatitude at pi/2
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "at theta = pi/2, upper index N = 1:\n";
  for (auto m : view[3].Orders()) {
    std::cout << "  X^1_{3," << std::setw(2) << m << "} = " << std::setw(10)
              << (view[3][m]) << "\n";
  }
  std::cout << "\n";

  // Two relations worth checking by hand. The library's own tests pin the
  // convention more tightly (tests/CheckWignerConvention.hpp compares the
  // l = 1 values with closed forms); these are the ones a reader can check at
  // any degree.
  //
  // The addition theorem, D&T (C.127): sum_m P^N_{lm} P^{N'}_{lm} = delta_{NN'}
  // at every colatitude. In terms of what is stored the right-hand side
  // carries the square of the normalisation, so it is
  // delta_{NN'} (2l+1)/4pi.
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

  // The involution P^{-N}_{l,-m} = (-1)^{m+N} P^N_{lm}, D&T (C.118), checked
  // over every stored (N, l, m) at one colatitude. MinusOneToPower(k) is the
  // library's (-1)^k for an integer k.
  //
  // This symmetry is not used to shrink the stored table. Its companion, the
  // reflection d^l_{Nm}(pi - theta) = (-1)^{l+N} d^l_{N,-m}(theta), is used,
  // by the matrix transform kernel only, which stores non-negative orders and
  // recovers the rest from it (example 20).
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

  // Write the values out: X^N_{6,2}(theta) for N = -2 .. 2, one column per
  // upper index, over the whole range of colatitude.
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

  // Plot it if gnuplot is there. The functions are real and have l - max(|m|,
  // |N|) zeros inside (0, pi). Near the poles d^l_{Nm} behaves like
  // (sin theta/2)^{|m-N|} at theta = 0 and (cos theta/2)^{|m+N|} at theta =
  // pi, so with m = 2 only the N = 2 curve is non-zero at the north pole and
  // only N = -2 at the south pole. Read the other way, at the north pole
  // only the order m = N survives, which is why a spin-N field there is
  // c exp(i N phi), a function of the direction of approach (example 21).
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
      script << (column > 2 ? ", " : "") << "'" << path
             << "' using 1:" << column << " with lines lw 2 title 'N = " << n
             << "'";
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
