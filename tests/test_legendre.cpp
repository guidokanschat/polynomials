#include "../src/legendre.hh" // File to test

#include <catch2/catch_all.hpp>
#include <vector>

#include "explizit_legendre_functions.hh"
#include "helpfunctions_for_tests.hh"

TEMPLATE_TEST_CASE("First 5 Legendre Polynom compute correct points",
    "[Legendre]", long double, double) {
  TestType x = GENERATE(take(10, random<TestType>(-100, 100)));
  int degree = GENERATE(0, 1, 2, 3, 4, 5);

  std::vector<TestType> coefficients;
  for (int i = 0; i < degree; ++i) {
    coefficients.push_back(0);
  }
  coefficients.push_back(1);
  Legendre polynome(coefficients);

  TestType truth = legendre(degree, x);
  TestType res = polynome(x);

  CAPTURE(degree);
  CAPTURE(res);
  CAPTURE(truth);
  REQUIRE_THAT(std::abs(truth - res) / std::abs(truth), Catch::Matchers::WithinAbs(0.0, tol<TestType>()));
}

TEMPLATE_TEST_CASE("Test compute_cartesian for first 5 Legendre functions", "[Legendre]", long double, double) {
  int degree = GENERATE(0, 1, 2, 3, 4, 5);
  std::vector<TestType> coefficients;
  for (int i = 0; i < degree; ++i) {
    coefficients.push_back(0);
  }
  coefficients.push_back(1);
  Legendre polynome(coefficients);

  std::vector<TestType> truth = legendre_polynome<TestType>(degree);
  std::vector<TestType> res = polynome.compute_cartesian();
  long double diff = 0;
  for (int i = 0; i < degree; ++i) {
    diff += std::abs(truth.at(i) - res.at(i));
  }

  CAPTURE(coefficients);
  CAPTURE(truth);
  CAPTURE(res);
  REQUIRE_THAT(diff, Catch::Matchers::WithinAbs(0.0, tol<TestType>()));
}

TEMPLATE_TEST_CASE("Test compute_cartesian for combinations of the first 5 Legendre functions", "[Legendre]",long double, double) {
  int degree = 5;
  std::vector<TestType> coefficients = GENERATE(take(10, chunk(5, random<TestType>(-100, 100))));
  Legendre polynome(coefficients);

  std::vector<TestType> truth(5, 0);
  for (int i = 0; i < degree; ++i) {
    std::vector<TestType> basis_vector = legendre_polynome<TestType>(i);
    for (int ii = 0; ii <= i; ++ii) {
      truth.at(ii) += coefficients.at(i) * basis_vector.at(ii);
    }
  }
  std::vector<TestType> res = polynome.compute_cartesian();
  long double diff = 0;
  for (int i = 0; i < degree; ++i) {
    diff += std::abs(truth.at(i) - res.at(i));
  }

  CAPTURE(coefficients);
  CAPTURE(truth);
  CAPTURE(res);
  REQUIRE_THAT(diff, Catch::Matchers::WithinAbs(0.0, tol<TestType>()));
}


TEMPLATE_TEST_CASE("Test general compute_cartesian", "[Legendre]", long double, double) {
  //int degree = GENERATE(take(5, random<int>(3, 12)));
  int degree = GENERATE(take(5, random<int>(3, 9)));
  std::vector<TestType> coefficients = GENERATE(take(10, chunk(12, random<TestType>(-100, 100))));
  coefficients.resize(degree);
  TestType x = GENERATE(take(10, random<TestType>(-100, 100)));

  Legendre polynome(coefficients);
  std::vector<TestType> cartesian = polynome.compute_cartesian();

  TestType res = polynome(x);
  TestType truth = horner_scheme(cartesian, x);

  CAPTURE(truth);
  CAPTURE(res);
  CAPTURE(coefficients);
  REQUIRE_THAT(std::abs(truth - res) / std::abs(truth), Catch::Matchers::WithinAbs(0.0, tol<TestType>()));
}

TEMPLATE_TEST_CASE("Test computation of grid points", "[Legendre]", long double, double) {
  int degree = GENERATE(take(5, random<int>(3, 12)));
  std::vector<TestType> coefficients = GENERATE(take(10, chunk(12, random<TestType>(-100, 100))));
  coefficients.resize(degree);
  TestType start = GENERATE(take(10, random<TestType>(-100, 100)));
  TestType end = GENERATE(take(10, random<TestType>(-100, 100)));
  //int number_of_points = GENERATE(take(10, random<TestType>(1, 100)));
  int number_of_points = 3;
  if (start > end) {
    std::swap(start, end);
  }
  TestType stepsize = (end - start) / static_cast<TestType>(number_of_points - 1);
  Legendre polynome(coefficients);
  std::vector<TestType> res = polynome.compute_grid(start, end, number_of_points);

  for (int i = 0; i < number_of_points; ++i) {
    TestType x = start + (i * stepsize);
    TestType truth = polynome(x);
    REQUIRE_THAT(std::abs(truth - res.at(i)) / std::abs(truth), Catch::Matchers::WithinAbs(0.0, tol<TestType>()));
  }
}
TEMPLATE_TEST_CASE("Test Integration of product of 2 from the first 5 polynomials", "[Legendre]", float, double) {
  int degree1 = GENERATE(0, 1, 2, 3, 4, 5);
  int degree2 = GENERATE(0, 1, 2, 3, 4, 5);
  std::vector<TestType> c1(degree1 + 1, 0);
  c1.at(degree1) = 1;
  Legendre p1(c1);
  std::vector<TestType> c2(degree2 + 1, 0);
  c2.at(degree2) = 1;
  Legendre p2(c2);

  TestType start = -1;
  TestType end = 1;
  if (degree1 != degree2) {
    SECTION("Integration of different polynomials") {
      std::vector<TestType> coefficients1 = legendre_polynome<TestType>(degree1);
      std::vector<TestType> coefficients2 = legendre_polynome<TestType>(degree2);
      std::vector<TestType> poly = compute_integral_polynoms(coefficients1, coefficients2);
      TestType truth = horner_scheme(poly, end) - horner_scheme(poly, start);
      TestType res = p1.IntegrateProduct(p2);
      REQUIRE_THAT(std::abs(truth - res), Catch::Matchers::WithinAbs(0.0, tol<TestType>()));
    }
  } else {
    SECTION("Integration of identical polynomials") {
      std::vector<TestType> coefficients1 = legendre_polynome<TestType>(degree1);
      std::vector<TestType> coefficients2 = legendre_polynome<TestType>(degree2);
      std::vector<TestType> poly = compute_integral_polynoms(coefficients1, coefficients2);
      TestType truth = horner_scheme(poly, end) - horner_scheme(poly, start);
      TestType res = p1.IntegrateProduct(p2);
      REQUIRE_THAT(std::abs(truth - res), Catch::Matchers::WithinAbs(0.0, tol<TestType>()));
    }
  }
}

TEMPLATE_TEST_CASE("Test Integration of linear combination from the first 6 polynomials", "[Legendre]", long double, double) {
  std::vector<TestType> base_integrals(6, 0);
  TestType start = -1;
  TestType end = 1;
  for (int i = 0; i < 6; ++i) {
      std::vector<TestType> coefficients1 = legendre_polynome<TestType>(i);
      std::vector<TestType> poly = compute_integral_polynoms(coefficients1, coefficients1);
      base_integrals.at(i) = horner_scheme<TestType>(poly, end) - horner_scheme<TestType>(poly, start);
  }
  std::vector<TestType> c1 = GENERATE(take(10, chunk(6, random<TestType>(-100, 100))));
  std::vector<TestType> c2 = GENERATE(take(10, chunk(6, random<TestType>(-100, 100))));

  Legendre p1(c1);
  Legendre p2(c2);
  TestType truth = 0;
  for (int i = 0; i < 6; ++i) {
    truth += base_integrals.at(i) * c1.at(i) * c2.at(i);
  }
  TestType res = p1.IntegrateProduct(p2);
  REQUIRE_THAT(std::abs(truth - res), Catch::Matchers::WithinAbs(0.0, tol<TestType>()));
}
