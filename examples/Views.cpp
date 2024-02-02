#include <GSHTrans/All>
#include <algorithm>
#include <cassert>
#include <concepts>
#include <functional>
#include <iostream>
#include <ranges>
#include <type_traits>
#include <vector>

class Dummy {
 public:
  void Print() { std::cout << "Hi\n"; }
};

template <std::ranges::view View>
class Container : public std::ranges::view_interface<Container<View>>,
                  public Dummy {
  using Int = std::ptrdiff_t;

 public:
  using value_type = std::ranges::range_value_t<View>;

  Container(View view, Int i) : _view{std::move(view)}, _i{i} {}

  Container(const Container&) = default;
  Container(Container&&) = default;

  Container& operator=(const Container&) = default;
  Container& operator=(Container&&) = default;

  template <std::ranges::view OtherView>
  requires std::convertible_to<std::ranges::range_value_t<OtherView>,
                               value_type>
  auto& operator=(OtherView&& view) {
    assert(view.size() == this->size());
    std::ranges::copy(view, _view.begin());
    return *this;
  }

  auto begin() { return _view.begin(); }
  auto end() { return _view.end(); }

  auto I() const { return _i; }

 private:
  View _view;
  Int _i;
};

class ToContainer : public std::ranges::range_adaptor_closure<ToContainer> {
  using Int = std::ptrdiff_t;

 public:
  ToContainer(Int i) : _i{i} {}

  template <std::ranges::view View>
  auto operator()(View v) {
    return Container(v, _i);
  }

 private:
  Int _i;
};

template <std::ranges::view View>
auto operator-(Container<View> v) {
  return v | std::ranges::views::transform([](auto x) { return -x; }) |
         ToContainer(v.I());
}

template <std::ranges::view View, GSHTrans::Field S>
auto operator*(Container<View> v, S s) {
  using T = std::ranges::range_value_t<View>;
  return v | std::ranges::views::transform([s](auto x) { return T(s) * x; }) |
         ToContainer(v.I());
}

template <std::ranges::view View, GSHTrans::Field S>
auto operator*(S s, Container<View> v) {
  return v * s;
}

template <std::ranges::random_access_range R>
Container(R&& r, std::ptrdiff_t i) -> Container<std::ranges::views::all_t<R>>;

int main() {
  auto x = std::vector<double>(10);

  auto i = 0;
  for (auto& val : x) val = i++;

  auto v = Container(x, 1);

  auto w = -v * 2;

  w.Print();

  std::cout << w.I() << std::endl;

  for (auto [val1, val2] : std::ranges::views::zip(2 * v, w))
    std::cout << val1 << " " << val2 << std::endl;
}