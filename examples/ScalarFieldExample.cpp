#include <GSHTrans/All>
#include <array>
#include <iostream>
#include <ranges>
#include <tuple>

using namespace GSHTrans;
int main() {
  using namespace GSHTrans;

  for (auto index : TensorIndices<10>()) {
    for (auto alpha : index) std::cout << alpha << " ";
    std::cout << std::endl;
  }
}