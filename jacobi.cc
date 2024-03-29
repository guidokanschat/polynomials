#include <array>
#include <iostream>

#include "jacobi.h"
#include "roots.h"

int main()
{
  PolynomialsJacobiInt<3,double> pols(0,0);

  auto val = pols.derivatives<0>(-1.);

  for (unsigned int i=0;i<val[0].size();++i)
    std::cout << i << "\t" << (long double)val[0][i] << std::endl;

  auto r = roots(pols);

  for (unsigned int i=0;i<r.size();++i)
    std::cout << (long double)r[i] << ' ';
  std::cout << std::endl;
}
