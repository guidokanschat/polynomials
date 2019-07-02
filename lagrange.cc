#include <array>
#include <iostream>

template <int degree, typename vtype>
class PolynomialsLagrange
{
private:
  std::array<vtype,degree+1> points;
  std::array<vtype,degree+1> weights;
public:
  template <typename array_type>
  PolynomialsLagrange(const array_type& pts)
  {
    for (unsigned int i=0;i<=degree;++i)
      {
	points[i] = pts[i];
	weights[i] = 1;
	for (unsigned int k=0; k<=degree;++k)
	  if (i!=k)
	    weights[i] *= pts[i] - pts[k];
	weights[i] = vtype(1.)/weights[i];
      }
  }
  template <int order>
  std::array<std::array<vtype,degree+1>,order+1>
  derivatives(vtype x) const
  {
    std::array<std::array<vtype,degree+1>,order+1> result;
    if (order >= 0)
      {
	vtype product = 1.;
	for (unsigned int i=0;i<=degree;++i)
	  product *= (x-points[i]);
	for (unsigned int i=0;i<=degree;++i)
	  if (x == points[i])
	    {
	      result[0][i] = 1;
	    }
	  else
	    result[0][i] = weights[i] * product / (x-points[i]);
      }
    return result;
  }
};

int main()
{
  std::array<double,3> pts{0., 0.5, 1.};
  PolynomialsLagrange<2,double> pols(pts);

  std::cout << pols.derivatives<0>(.2)[0][0] << " "
	    << pols.derivatives<0>(.2)[0][1] << " "
	    << pols.derivatives<0>(.3)[0][1] << std::endl;
}
