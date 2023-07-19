
#include <cmath>
#include <concepts>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <tuple>


#include <GSHT/Core>


int main() {
  using Float = double;

  using std::cout;
  using std::endl;
  cout.setf(std::ios_base::scientific);
  cout.precision(std::numeric_limits<Float>::digits10);  
  using namespace GSHT;


  int L = 5;
  int M = 3;
  int N = 2;
  Float theta = 0.2;
  WignerValues d(L,M,N,theta);

  
  for (int n : UpperIndices(0,N))
    {
      
      for (int l : Degrees(n,L))
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
