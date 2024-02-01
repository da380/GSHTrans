#include <GSHTrans/All>
#include <algorithm>
#include <cassert>
#include <concepts>
#include <functional>
#include <iostream>
#include <ranges>
#include <type_traits>
#include <vector>

template <typename Iter>
class Container : public std::ranges::view_interface<Container<Iter>> {
  using Int = std::ptrdiff_t;

 public:
  using value_type = std::iter_value_t<Iter>;

  Container(Iter start, Iter finish) : _start{start}, _finish{finish} {}

  Container(Iter start, Int size) : Container(start, std::next(start, size)) {}

  Container(const Container&) = default;
  Container(Container&&) = default;

  Container& operator=(const Container&) = default;
  Container& operator=(Container&&) = default;

  template <std::ranges::view View>
  requires std::convertible_to<std::ranges::range_value_t<View>, value_type>
  auto& operator=(View&& view) {
    assert(view.size() == this->size());
    std::ranges::copy(view, _start);
    return *this;
  }

  auto begin() { return _start; }
  auto end() { return _finish; }

 private:
  Iter _start;
  Iter _finish;
};

int main() {
  using namespace GSHTrans::Views;

  auto x = std::vector<double>(10);

  auto i = 0;
  for (auto& val : x) val = i++;

  auto v = Container(x.begin(), x.end());

  for (auto [x1, x2] : std::ranges::views::zip(x, 4 * v / 7))
    std::cout << x1 << " " << x2 << std::endl;
}