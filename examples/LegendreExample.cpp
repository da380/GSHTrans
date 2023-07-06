
#include <iostream>
#include <numbers>

#include "Legendre.h"

int main() {

  using Float = double;
  
  Float theta = 0.3*std::numbers::pi_v<Float>;

  GSHT::LegendreValue p(theta);

  using std::cout;
  using std::endl;
  cout.setf(std::ios_base::scientific);
  cout.setf(std::ios_base::showpos);
  cout.precision(16);

  int L = 4;
  for(int l = 0; l <= L; l++)
    {
      ++p;
      for (int m = 0; m <= l; m++)
	cout << p(m) << " ";
      cout << endl;
    }
  
  
  
  
}
