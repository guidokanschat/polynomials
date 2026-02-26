// clang-format off
#include "../src/fraction.hh"
// clang-format on

#include <catch2/catch_all.hpp>

#include "helpfunctions_for_tests.hh"

TEMPLATE_TEST_CASE("Correct initialization of the fraction", "[FRACTION]", int,
                   long, short) {
  TestType a = GENERATE(take(10, random<TestType>(-100, 100)));
  TestType b = GENERATE(take(
      10, filter([](int i) { return i != 0; }, random<TestType>(-100, 100))));
  Fraction x(a, b);

  REQUIRE(a == x.Numerator);
  REQUIRE(b == x.Denominator);
}

TEMPLATE_TEST_CASE("Correct conversion to double", "[FRACTION]", int, long,
                   short) {
  TestType a = GENERATE(take(10, random<TestType>(-100, 100)));
  TestType b = GENERATE(take(
      10, filter([](int i) { return i != 0; }, random<TestType>(-100, 100))));
  Fraction x(a, b);

  REQUIRE_THAT(std::abs((double(a) / double(b)) - static_cast<double>(x)),
               Catch::Matchers::WithinAbs(0.0, 1e-16));
}

TEMPLATE_TEST_CASE("Correct reduction of the fraction", "[FRACTION]", int, long,
                   short) {
  TestType a = GENERATE(take(10, random<TestType>(-100, 100)));
  TestType b = GENERATE(take(
      10, filter([](int i) { return i != 0; }, random<TestType>(-100, 100))));
  Fraction x(a, b);
  double base = static_cast<double>(x);
  x.reduce();

  REQUIRE_THAT(std::abs(base - static_cast<double>(x)),
               Catch::Matchers::WithinAbs(0.0, 1e-14));
}

TEMPLATE_TEST_CASE("Correct Addition of fractions, with identical Denominator",
                   "[FRACTION]", int, long, short) {
  TestType a = GENERATE(take(10, random<TestType>(-100, 100)));
  TestType b = GENERATE(take(
      10, filter([](int i) { return i != 0; }, random<TestType>(-100, 100))));
  TestType c = GENERATE(take(10, random<TestType>(-100, 100)));
  Fraction x(a, b);
  Fraction y(c, b);

  x += y;

  REQUIRE((a + c) == x.Numerator);
  REQUIRE(b == x.Denominator);
}

TEMPLATE_TEST_CASE("Correct Addition of fractions, with different Denominator",
                   "[FRACTION]", int, long, short) {
  TestType a = GENERATE(take(10, random<TestType>(-100, 100)));
  TestType b = GENERATE(take(
      10, filter([](int i) { return i != 0; }, random<TestType>(-100, 100))));
  TestType c = GENERATE(take(10, random<TestType>(-100, 100)));
  TestType d = GENERATE(take(
      10, filter([](int i) { return i != 0; }, random<TestType>(-100, 100))));

  if (b == d) {
    d += 1;
  }

  Fraction x(a, b);
  Fraction y(c, d);

  x += y;

  REQUIRE(((a * d) + (b * c)) == x.Numerator);
  REQUIRE((b * d) == x.Denominator);
}

TEMPLATE_TEST_CASE("Correct Multiplication of fractions", "[FRACTION]", int,
                   long, short) {
  TestType a = GENERATE(take(10, random<TestType>(-100, 100)));
  TestType b = GENERATE(take(
      10, filter([](int i) { return i != 0; }, random<TestType>(-100, 100))));
  TestType c = GENERATE(take(10, random<TestType>(-100, 100)));
  TestType d = GENERATE(take(
      10, filter([](int i) { return i != 0; }, random<TestType>(-100, 100))));
  Fraction x(a, b);
  Fraction y(c, d);

  x *= y;

  REQUIRE((a * c) == x.Numerator);
  REQUIRE((b * d) == x.Denominator);
}

TEMPLATE_TEST_CASE("Correct Division of fractions", "[FRACTION]", int, long,
                   short) {
  TestType a = GENERATE(take(10, random<TestType>(-100, 100)));
  TestType b = GENERATE(take(
      10, filter([](int i) { return i != 0; }, random<TestType>(-100, 100))));
  TestType c = GENERATE(take(10, random<TestType>(-100, 100)));
  TestType d = GENERATE(take(
      10, filter([](int i) { return i != 0; }, random<TestType>(-100, 100))));
  Fraction x(a, b);
  Fraction y(c, d);

  x /= y;

  REQUIRE((a * d) == x.Numerator);
  REQUIRE((b * c) == x.Denominator);
}
