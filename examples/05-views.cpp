// 05 -- Views
//
// A field need not own its storage. A view is a pointer into someone else's
// buffer plus a grid handle, and it is admissible everywhere an owning field
// is -- which is what makes a tensor able to hold one buffer and hand out its
// components, and what makes a 3D field able to hand out radial slices.
//
// What this shows
//   Wrapping a caller-owned buffer as a field, contiguous or strided,
//   writable or read-only, and using it in the field algebra.
//
// Assumes
//   Examples 02-04 (fields, upper indices, expressions).
//
// Introduced
//   SpinFieldView<N, Grid>(grid, span, stride), ConstSpinFieldView.
//
// Output
//   An inner product of two views, and the caller's buffer showing the
//   writes made through them.
//
// The tensor fields of example 08 hand out their components as views like
// these, and the layered fields of example 16 hand out radial slices.

#include <GSHTrans/GSHTrans.hpp>
#include <cmath>
#include <complex>
#include <iostream>
#include <span>
#include <vector>

int main() {
  using namespace GSHTrans;

  using Real = double;
  using Complex = std::complex<Real>;
  using Grid = GaussLegendreGrid<Real, All, All>;

  auto grid = Grid(8, 2);
  const auto size = static_cast<std::size_t>(grid.FieldSize());

  // Storage the library did not allocate: two fields laid end to end, as a
  // caller with their own memory would have them. Each field's samples are in
  // the library's order, iTheta * nPhi + iPhi.
  auto storage = std::vector<Complex>(2 * size);
  for (std::size_t i = 0; i < storage.size(); i++) {
    storage[i] = Complex{std::cos(0.1 * i), std::sin(0.2 * i)};
  }

  // Each half is a field of upper index 1. The view checks, and throws
  // std::invalid_argument otherwise, that the grid carries its upper index
  // and that the span covers every grid point. It does not own the buffer,
  // which must outlive it.
  auto first = SpinFieldView<1, Grid>(grid, std::span(storage).first(size));
  auto second = SpinFieldView<1, Grid>(grid, std::span(storage).last(size));

  static_assert(SpinWeighted<decltype(first)>);
  static_assert(decltype(first)::UpperIndex == 1);

  // They compose with owning fields and with each other. A view is a handle
  // (a grid handle, a span and a stride), so an expression holds it by value
  // and stays valid if the named view goes out of scope first; it is the
  // buffer underneath that must stay alive.
  std::cout << "<first, second> = " << Integrate(conj(first) * second) << "\n";

  // Writing through a view writes the caller's buffer.
  first[0, 0] = Complex{42.0, 0.0};
  std::cout << "storage[0] is now " << storage[0] << "\n";

  // A view may also be strided, which is how a tensor stored point by point
  // hands out a component whose samples are not contiguous. The stride is a
  // constructor argument and defaults to one. Here three fields are
  // interleaved and the view takes the middle one: sample k of the view is
  // element 1 + 3 k of the buffer.
  auto interleaved = std::vector<Complex>(3 * size);
  auto middle =
      SpinFieldView<1, Grid>(grid, std::span(interleaved).subspan(1), 3);
  middle[0, 0] = Complex{1.0, 0.0};
  middle[0, 1] = Complex{2.0, 0.0};
  std::cout << "strided view wrote elements 1 and 4: " << interleaved[1] << " "
            << interleaved[4] << "\n";

  // A read-only view over storage nobody may write through: it takes a span
  // of const elements, and has no writable operator[].
  auto locked = ConstSpinFieldView<1, Grid>(
      grid, std::span<const Complex>(storage).first(size));
  std::cout << "read-only view at the same point " << (locked[0, 0]) << "\n";
}
