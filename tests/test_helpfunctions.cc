// clang-format off
#include "helpfunctions_for_tests.hh"
// clang-format on

#include <catch2/catch_all.hpp>

#include "explizit_legendre_functions.hh"


TEMPLATE_TEST_CASE("Correct computation using horner scheme Legendre 5", "[FRACTION]",
                   double) {
  TestType x = GENERATE(take(10, random<TestType>(-100, 100)));
  std::vector<TestType> c = {0, 15.0/8.0, 0, -70.0/8.0, 0, 63.0/8.0};
  TestType res = horner_scheme(c, x);

  REQUIRE_THAT(std::abs(legendre<5>(x) - res), Catch::Matchers::WithinAbs(0.0,
                        1e-4));
}


TEMPLATE_TEST_CASE("Correct computation using horner scheme Legendre 2", "[FRACTION]",
                   double) {
  TestType x = GENERATE(take(10, random<TestType>(-100, 100)));
  std::vector<TestType> c = {-0.5, 0, 3.0/2.0};
  TestType res = horner_scheme(c, x);


  REQUIRE_THAT(std::abs(legendre<2>(x) - res), Catch::Matchers::WithinAbs(0.0,
                        1e-10));
}
