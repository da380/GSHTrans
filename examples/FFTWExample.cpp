
#include <complex>
#include <iostream>
#include <iterator>
#include <list>
#include <vector>

#include "FFTW.h"





int main()
{

  int n = 10;
  
  std::vector<std::complex<double>> in(n), out(n);
  
  FFTW::Plan plan(n,in.begin(),out.begin(),FFTW::DirectionFlag::Forward,
	     FFTW::PlanFlag::Estimate);

  std::cout << plan.GetNorm() << std::endl;


  
}

