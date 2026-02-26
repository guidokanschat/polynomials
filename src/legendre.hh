#ifndef NEW_LEGENDRE_HH_
#define NEW_LEGENDRE_HH_

#include <algorithm>
#include <array>
#include <vector>

#include "fraction.hh"

/* Class Implementation for linear combination of legendre polynoms
 * It allows computation of the product of two such objects, and their
 * integration over the interval [-1, 1]
 */

template <typename Vtype>
class Legendre {
 private:
  typedef Fraction<int> Frac;
  std::vector<Vtype> coefficients_;
  int max_degree_;

 public:
  Legendre(const std::vector<Vtype>& kCoeff)
      : coefficients_(kCoeff), max_degree_(kCoeff.size() - 1) {}

  const std::vector<Vtype> compute_grid(Vtype start, Vtype end,
                                        int number_of_points) {
    std::vector<Vtype> res(number_of_points);
    double stepsize = (end - start) / (number_of_points - 1);
    Vtype x = start;
    for (int i = 0; i < number_of_points; ++i, x += stepsize) {
      res.at(i) = operator()(x);
    }
    return res;
  }

  /*
   * Compute the function value for a given point, based on the
   * 3-term recursion:
   * (n + 1) P_{n+1}(x) = (2n + 1)x * P_{n}(x) - n P_{n-1}(x)
   */
  const Vtype operator()(const Vtype kX) {
    if (max_degree_ == 0) {
      return coefficients_.at(0);
    }
    // TODO: replace with a better data structure
    // Store the last 3 points for the 3 term recursion
    std::array<Vtype, 3> legendre_points;
    Vtype res = 0;

    legendre_points.at(0) = 1;
    legendre_points.at(1) = kX;
    for (int i = 1; i < max_degree_; ++i) {
      legendre_points.at(2) = (((2 * i + 1) * kX * legendre_points.at(1)) -
                               (i * legendre_points.at(0))) /
                              (i + 1);
      res += coefficients_.at(i - 1) * legendre_points.at(0);
      legendre_points.at(0) = legendre_points.at(1);
      legendre_points.at(1) = legendre_points.at(2);
    }
    res += legendre_points.at(0) * coefficients_.at(max_degree_ - 1) +
           legendre_points.at(1) * coefficients_.at(max_degree_);
    return res;
  }

  /*
   * Compute the Integral [-1, 1], based on the formula
   * Integral[-1, 1] = 2 / (2 * n + 1) * kroniker_{nm}
   */
  const Vtype IntegrateProduct(const Legendre<Vtype>& kF) {
    Vtype res = 0;
    int dim = std::min(max_degree_, kF.max_degree_);
    for (int i = 0; i <= dim; ++i) {
      res += (2.0 / (2.0 * i + 1.0)) * coefficients_.at(i) *
             kF.coefficients_.at(i);
    }
    return res;
  }

  static Frac smaller_k(int n, int k) {
    return Frac((-(n - k + 2) * (n + k - 1)), (k * (k - 1)));
  }

  static Frac smaller_n(int n, int k) { return Frac(-(n + k - 1), (n - k)); }

  // na_{n, i} = (2n - 1) a_{n - 1, i - 1(x) - (n - 1)_{a_{n-2, 2}}
  std::vector<Vtype> compute_cartesian() {
    std::vector<Vtype> res(max_degree_ + 1, 0);
    // The factor for the smallest  k with non zero 0 entry
    Fraction current_base_factor(1, 1);
    res.at(0) += coefficients_.at(0);
    if (max_degree_ > 0) {
      res.at(1) += coefficients_.at(1);

      for (int n = 2; n <= max_degree_; n += 2) {
        current_base_factor *= Fraction(-n + 1, n);
        Fraction a = current_base_factor;
        Vtype c = coefficients_.at(n);
        res.at(0) += c * static_cast<Vtype>(a);
        for (int k = 2; k <= n; k += 2) {
          a *= smaller_k(n, k);
          res.at(k) += static_cast<Vtype>(a) * c;
        }
      }

      current_base_factor = Fraction(1, 1);
      for (int n = 3; n <= max_degree_; n += 2) {
        current_base_factor *= Fraction(-n, n - 1);
        Fraction a = current_base_factor;
        Vtype c = coefficients_.at(n);
        res.at(1) += c * static_cast<Vtype>(a);
        for (int k = 3; k <= n; k += 2) {
          a *= smaller_k(n, k);
          res.at(k) += static_cast<Vtype>(a) * c;
        }
      }
    }
    return res;
  }
};  // Legendre

#endif
