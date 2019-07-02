#include <array>
#include <iostream>

#include "lagrange.h"

int main()
{
  std::array<double,3> pts{0., 0.5, 1.};
  PolynomialsLagrange<2,double> pols(pts);

  std::cout << pols.derivatives<0>(.2)[0][0] << " "
	    << pols.derivatives<0>(.2)[0][1] << " "
	    << pols.derivatives<0>(.3)[0][1] << std::endl;
}
