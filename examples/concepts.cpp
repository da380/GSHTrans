
#include <iostream>

#include "NumericConcepts/Numeric.hpp"

template <NumericConcepts::Numeric T>
void print_numeric_value(T value) {
  std::cout << "Numeric value: " << value << std::endl;
}

int main() {
  print_numeric_value(42);     // OK: int is Integral
  print_numeric_value(3.14f);  // OK: float is Real
  print_numeric_value(
      std::complex(1.0, 2.0));  // OK: std::complex<double> is Complex
  // print_numeric_value("hello");          // Fails: const char* is not Numeric
}
