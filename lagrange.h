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
        if (degree > 0)
          {
            std::array<vtype,degree> accumulated_product;
            accumulated_product[0] = x - points[0];
            result[0][degree-1] = x-points[degree];
            for (unsigned int d=1; d<degree; ++d)
              {
                accumulated_product[d] = accumulated_product[d-1]*(x-points[d]);
                result[0][degree-d-1] = result[0][degree-d]*(x-points[degree-d]);
              }
            result[0][0] *= weights[0];
            for (unsigned int d=1; d<degree; ++d)
              result[0][d] *= weights[d] * accumulated_product[d-1];
            result[0][degree] = weights[degree] * accumulated_product[degree-1];
          }
        else
          result[0][0] = weights[0];
      }
    return result;
  }
};
