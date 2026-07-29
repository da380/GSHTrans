
#include <iostream>

#include "GSHTrans/All"

int main() {
  auto s = GSHTrans::Wigner3jSymbol<double>(10, 0, 10, 0, 0, 0);

  std::cout << s << std::endl;
}
