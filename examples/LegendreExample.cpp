
#include <concepts>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>


#include "GSHT.h"


int main() {
  using Float = double;

  using std::cout;
  using std::endl;
  cout.setf(std::ios_base::scientific);
  cout.setf(std::ios_base::showpos);
  cout.precision(std::numeric_limits<Float>::digits10);

  
  using namespace GSHT;

  Float theta = 0.2;
  int L = 2;
  int M = L;
  
  LegendreValues p(theta,L,M);

  for(int  l = 0; l <= L; l++)
    {

      for(int m = 0; m <= std::min(l,M); m++)
	{
	  cout << p(l,m) << " ";
	}

      cout << endl;
      
    }
  


}
