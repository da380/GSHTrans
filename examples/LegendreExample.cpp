
#include <concepts>
#include <fstream>
#include <iostream>
#include <numbers>


#include "GSHT.h"


int main() {
  using Float = double;

  using std::cout;
  using std::endl;
  cout.setf(std::ios_base::scientific);
  cout.setf(std::ios_base::showpos);
  cout.precision(16);

  int n = 100;
  int l = 200;
  std::vector<Float> thetas(n);  

  Float th{0.0};
  Float dth = std::numbers::pi_v<Float>/static_cast<Float>(n-1);
  for (auto& val : thetas)
    {
      val = th;
      th += dth;
    }
  

  for(auto theta : thetas)
    {

      
    // Construct the iterator.
    GSHT::LegendreIterator p(theta);

    // Take l steps.
    p + l;

    // Get a reference the current values
    GSHT::LegendreIterator<Float>::reference plm = *p;

    // Build vector of values using STL function.
    std::vector<Float> plm2(l + 1);
    std::generate(plm2.begin(), plm2.end(), [l, m = 0, theta]() mutable {
      return std::sph_legendre(l, m++, theta);
    });

    // Compute the difference.
    std::transform(plm.begin(), plm.end(), plm2.begin(), plm.begin(),
                   std::minus<>());

    // Write out the maximum value
    auto max = std::abs(*std::max_element(plm.begin(), plm.end()));

    cout << theta  << " " << max << endl;
  }
}
