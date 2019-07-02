#include <array>
#include <iostream>

#include "lagrange.h"
#include "roots.h"
#include "jacobi.h"

int main()
{
  auto r = roots(PolynomialsJacobiInt<5,double>(1,1));
  std::array<double,7> pts;
  pts[0] = 0;
  for (unsigned int i=0; i<r.size(); ++i)
    pts[i+1] = 0.5+0.5*r[i];
  pts[6] = 1;
  std::cout << "roots: ";
  for (auto r : pts)
    std::cout << r << " ";
  std::cout << std::endl;
  PolynomialsLagrange<6,double> pols(pts);

  for (unsigned int i=0; i<=20; ++i)
    {
      const double x = (double)i/20.;
      for (unsigned int d=0; d<pts.size(); ++d)
        std::cout << pols.derivatives<0>(x)[0][d] << " ";
      std::cout << std::endl;
    }
}
