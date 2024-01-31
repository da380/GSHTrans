#include <concepts>
#include <functional>
#include <iostream>
#include <ranges>
#include <type_traits>
#include <vector>

template <typename T>
class AffineTransformation {
 public:
  AffineTransformation() : _scale{1}, _shift{0} {}
  AffineTransformation(T scale, T shift) : _scale{scale}, _shift{shift} {}
  // AffineTransformation(AffineTransformation&) = default;
  //  AffineTransformation(AffineTransformation&&) = default;
  T operator()(T x) const { return _scale * x + _shift; }

 private:
  T _scale;
  T _shift;
};

template <typename T>
class Vector {
  using Int = std::ptrdiff_t;

 public:
  explicit Vector(Int n) : _data{Vector(n)} {}

 private:
  std::vector<T> _data;
};

template <std::ranges::view View, typename Function>
auto UnaryTransform(View view, Function f) {
  return view | std::ranges::views::transform(f);
}

template <std::ranges::view View, typename T>
auto operator+(View view, T t) {
  return UnaryTransform(view, AffineTransformation<T>(1, t));
  // return UnaryTransform(view, [t](auto x) { return x + t; });
}

int main() {
  auto x = std::vector<double>(10);

  auto view = std::ranges::views::all(x);

  for (auto val : (view + 1) + 3) std::cout << val << std::endl;
}