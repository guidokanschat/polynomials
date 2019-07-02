#include <cmath>
#include <limits>

namespace internal
{
  template <typename Number>
  Number get_epsilon()
  {
    return std::numeric_limits<Number>::epsilon();
  }

  template <>
  __float128 get_epsilon<__float128>()
  {
    return std::pow(2,-115);
  }
}

template <class Poly>
std::array<typename Poly::value_type,Poly::degree>
roots(const Poly& p)
{
  using vtype = typename Poly::value_type;
  const unsigned int degree = Poly::degree;
  std::array<vtype,degree> result;

  auto do_newton = [&](const vtype initial_guess,
                       const vtype *start_array,
                       const vtype *end_array)
    {
      vtype x = initial_guess;
      return x;
    };


  vtype x = -1.;
  for (unsigned int d=0;d<degree;++d)
    {
      while (true)
	{
          vtype s = 0.;
          for (unsigned int i=0; i<d; ++i)
            s += 1. / (x - result[i]);

	  auto val = p.template derivatives<1>(x);
          const vtype update = val[0][degree]/(val[0][degree] * s - val[1][degree]);
	  x += update;
	  if (std::abs(update) < 8. * internal::get_epsilon<vtype>())
            break;
	}
      result[d] = x;
      if (d<degree-1 && d>(degree-1)/2)
        x = -result[degree-d-2];
      else if (d>0)
	x = 2.*result[d]-1.*result[d-1];
      else
	x = 2.*result[d]+1.;
    }
  return result;
}
