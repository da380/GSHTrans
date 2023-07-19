
#include <cmath>
#include <concepts>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>


#include <GSHT/Core>


int main() {
  using Float = double;

  using std::cout;
  using std::endl;
  cout.setf(std::ios_base::scientific);
  cout.setf(std::ios_base::showpos);
  cout.precision(std::numeric_limits<Float>::digits10);

  
  using namespace GSHT;

  Float theta = 0.45 * std::numbers::pi_v<Float>;
  int L = 1000;
  int M = L;
  
  LegendreValues p(theta,L,M);

  Float max = 0.0;
  for(int  l = 0; l <= L; l++)
    {

      for(int m = 0; m <= std::min(l,M); m++)
	{
	  if(auto error = std::abs(p(l,m) - std::sph_legendre(l,m,theta)); error > max)
	    max = error;
	}
    }
  cout << max << endl;
 


}
