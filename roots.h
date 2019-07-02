#include <cmath>

template <class Poly>
std::array<typename Poly::value_type,Poly::degree>
roots(const Poly& p)
{
  const unsigned int degree = Poly::degree;
  std::array<typename Poly::value_type,Poly::degree> result;
  
  typename Poly::value_type x = -1.;
  for (unsigned int d=0;d<degree;++d)
    {
      while (true)
	{
	  auto val = p.template derivatives<1>(x);
	  if (std::fabs(val[0][degree]) < 1.e-15)
	    {
	      std::cerr << "root " << d << " = " << x << std::endl;
	      result[d] = x;
	      break;
	    }
	  x -= val[0][degree]/val[1][degree];
	}
      if (d>0)
	x = 2.5*x-1.5*result[d-1];
      else
	x = 2.5*x+1.5;
      std::cerr << "seed = " << x << std::endl;      
    }
  return result;
}
