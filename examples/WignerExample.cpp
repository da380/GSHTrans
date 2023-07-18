
#include <cmath>
#include <concepts>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <tuple>


#include "GSHT.h"


int main() {
  using Float = double;

  using std::cout;
  using std::endl;
  cout.setf(std::ios_base::scientific);
  cout.precision(std::numeric_limits<Float>::digits10);  
  using namespace GSHT;


  int L = 5;
  int M = 5;
  int N = 2;
  Float theta = 0.2;
  WignerValues d(L,M,N,theta);

  
  for ( int n = 0; n <= N; n++)
    {
      
      for (int l  = n; l <= L; l++)
	{

	  auto start = d.beginForDegreeAtUpperIndex(l,n);
	  auto finish = d.endForDegreeAtUpperIndex(l,n);
	  while(start != finish){
	    auto [nn,ll,mm] = *start++;
	    cout << nn << " " << ll << " " << mm  << endl;
	  }
	  cout << "---------------------------------------\n";
	  
	  
	  
	}
      cout << "==============================================\n";

      
    }

  
  

}
