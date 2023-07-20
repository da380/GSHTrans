
#include <cmath>
#include <concepts>
#include <execution>
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


  int L = 10;
  int M = 4;
  int N = 2;
  Float theta = 0.2;
  WignerValues d(L,M,N,theta);




}
