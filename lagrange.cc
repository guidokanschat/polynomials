#include <array>
#include <iostream>

#include "lagrange.h"

int main()
{
  std::array<double,5> pts{0., 0.25, 0.5, 0.75, 1.};
  PolynomialsLagrange<4,double> pols(pts);

  for (unsigned int i=0; i<=20; ++i)
    {
      const double x = (double)i/20.;
      for (unsigned int d=0; d<pts.size(); ++d)
        std::cout << pols.derivatives<0>(x)[0][d] << " ";
      std::cout << std::endl;
    }
}
