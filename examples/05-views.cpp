// 05 -- Views
//
// A field need not own its storage. A view is a pointer into someone else's
// buffer plus a grid handle, and it is admissible everywhere an owning field
// is -- which is what makes a tensor able to hold one buffer and hand out its
// components, and what makes a 3D field able to hand out radial slices.

#include <GSHTrans/All>
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
  // caller with their own memory would have them.
  auto storage = std::vector<Complex>(2 * size);
  for (std::size_t i = 0; i < storage.size(); i++) {
    storage[i] = Complex{std::cos(0.1 * i), std::sin(0.2 * i)};
  }

  // Each half is a field of upper index 1.
  auto first = SpinFieldView<1, Grid>(grid, std::span(storage).first(size));
  auto second = SpinFieldView<1, Grid>(grid, std::span(storage).last(size));

  static_assert(SpinWeighted<decltype(first)>);
  static_assert(decltype(first)::UpperIndex == 1);

  // They compose with owning fields and with each other, and they are
  // terminals, so an expression holds an lvalue one by reference.
  std::cout << "<first, second> = "
            << Integrate(conj(first) * second) << "\n";

  // Writing through a view writes the caller's buffer.
  first[0, 0] = Complex{42.0, 0.0};
  std::cout << "storage[0] is now " << storage[0] << "\n";

  // A view may also be strided, which is how a tensor stored point by point
  // hands out a component whose samples are not contiguous. The stride is a
  // constructor argument and defaults to one.
  auto interleaved = std::vector<Complex>(3 * size);
  auto middle = SpinFieldView<1, Grid>(grid, std::span(interleaved).subspan(1),
                                       3);
  middle[0, 0] = Complex{1.0, 0.0};
  middle[0, 1] = Complex{2.0, 0.0};
  std::cout << "strided view wrote elements 1 and 4: " << interleaved[1] << " "
            << interleaved[4] << "\n";

  // A read-only view over storage nobody may write through.
  auto locked = ConstSpinFieldView<1, Grid>(
      grid, std::span<const Complex>(storage).first(size));
  std::cout << "read-only view at the same point " << (locked[0, 0]) << "\n";
}
