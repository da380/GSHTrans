
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
  using namespace GSHT;

  Float theta = 10;
  int L = 4;
  int M = 1;
  
  LegendreValues plm(theta,L,M);


  for(int l = 0; l <= L; l++)
    {
      LegendreValues<Float>::iterator start = plm.begin(l);
      LegendreValues<Float>::iterator finish = plm.end(l);
      if(l < L)
	{	  
	  cout << *start <<  " " << *finish << endl;
	}
      else
	{
	  cout << *start << endl;
	}
    }



}
