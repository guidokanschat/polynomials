#include "../src/new_legendre.hh"
#include "../src/fraction.hh"
#include "unoptimized/not_optimzied_legerndre.hh"
#include "helpfunctions.hh"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>


typedef double DTYPE;


int main() {
  Fraction<int> x{3, 5};
  Fraction<int> y{5, 2};

  std::cout << "x = " <<  x << std::endl;
  std::cout << "y = " <<  y << std::endl;
  std::cout << "x + y = " << x + y << std::endl;
  std::cout << "x * y = " << x * y << std::endl;
  std::cout << "x / y = " << x / y << std::endl;
  std::cout << "y / x = " << y / x << std::endl;
  std::cout << "x = " <<  x << std::endl;
  std::cout << "y = " <<  y << std::endl;

  //std::vector<double> x{1, 2, 0.5};
  //std::vector<double> x{1.7, 5.5, 0.3};
  //std::vector<DTYPE> x{1.7, 5.5, 0.3, -2.74};
//  std::vector<DTYPE> x{4.77, 5.75, 1, 1};
//  std::vector<DTYPE> x2{3.3, 1.5};
//
//  std::cout << "Input coefficients = " << x << std::endl;
//  TestLegendre<DTYPE> test{x};
//  NewLegendre<DTYPE> new_test{x};
//  NewLegendre<DTYPE> new_test2{x2};
//  std::cout << "result f(1.0):" << test(1.0) << std::endl;
//  std::cout << "result f(2.0):" << test(2.0) << std::endl;
//  std::cout << "result f(-2.0):" << test(-2.0) << std::endl;
//  std::cout << "new result f(1.0):" << new_test(1.0) << std::endl;
//  std::cout << "new result f(2.0):" << new_test(2.0) << std::endl;
//  std::cout << "new result f(-2.0):" << new_test(-2.0) << std::endl;
//
//  std::cout << new_test.ComputeGrid(-2, 2, 100) << std::endl;
//  std::cout << new_test.IntegrateProduct(new_test2) << std::endl;



  return 0;
}
