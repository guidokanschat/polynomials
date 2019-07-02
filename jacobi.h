template <int deg, typename vtype=double>
class PolynomialsJacobiInt
{
private:
  const int a;
  const int b;
public:
  typedef vtype value_type;
  static const int degree = deg;
  
  constexpr PolynomialsJacobiInt(const int a, const int b)
    : a(a), b(b)
  {}
  
  template <int order>
  std::array<std::array<vtype,degree+1>,order+1>
  derivatives(vtype x) const
  {
    std::array<std::array<vtype,degree+1>,order+1> result;
    if (order >= 0)
      {
	std::array<vtype,degree+1>& r = result[0];
	r[0] = 1.;
	if (degree > 0)
	  r[1] = .5*((2+a+b) * x + a-b);
	for (unsigned int d=2;d<=degree;++d)
	  {
	    r[d] =
	      ((2*d+a+b-1)*((2*d+a+b)*(2*d+a+b-2)*x+a*a-b*b)*r[d-1]
	       - 2*(d+a-1)*(d+b-1)*(2*d+a+b)*r[d-2])
	      /((2*d+a+b-2)*2*d*(d+a+b));
	  }
      }
    if (order >= 1)
      {
	std::array<vtype,degree+1>& r = result[1];
	PolynomialsJacobiInt<degree,value_type> opoly(a+1,b+1);
	auto other_values = opoly.template derivatives<0>(x);
	r[0] = 0.;
	for (unsigned int d=1;d<=degree;++d)
	  {
	    r[d] = .5*(1+a+b+d)*other_values[0][d-1];
	  }
      }
    return result;
  }
};
