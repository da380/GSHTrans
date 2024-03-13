#ifndef GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H
#define GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H

#include <algorithm>
#include <complex>
#include <iostream>
#include <memory>
#include <ranges>

#include "Concepts.h"
#include "GridBase.h"
#include "Indexing.h"

namespace GSHTrans {

template <typename Derived>
class CanonicalCoefficientBase {
  using Int = std::ptrdiff_t;

 public:
  // Data access methods.
  auto View() const { return _Derived().View(); }
  auto View() { return _Derived().View(); }
  auto begin() { return _Derived().View().begin(); }
  auto end() { return _Derived().View().end(); }
  auto size() { return _Derived().View().size(); }
  auto UpperIndex() const { return _Derived().UpperIndex(); }

  // Grid access methods.
  auto Grid() const { return _Derived().Grid(); }

  auto UpperIndex() const { return _Derived().UpperIndex(); }
  auto MaxDegree() const { return Grid().MaxDegree(); }

 private:
  auto& _Derived() const { return static_cast<const Derived&>(*this); }
  auto& _Derived() { return static_cast<Derived&>(*this); }
};

}  // namespace GSHTrans

#endif  // GSH_TRANS_CANONICAL_COEFFICIENTS_GUARD_H
