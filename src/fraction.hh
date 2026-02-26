#ifndef FRACTION_HH_
#define FRACTION_HH_

#include <iostream>
#include <numeric>

template <typename T>
class Fraction {
 public:
  T Numerator;
  T Denominator;

  Fraction(T a, T b) : Numerator{a}, Denominator{b} {}
  Fraction() : Numerator{0}, Denominator{1} {}

  explicit operator double() const {
    return (static_cast<double>(Numerator) / static_cast<double>(Denominator));
  }

  explicit operator long double() const {
    return (static_cast<long double>(Numerator) / static_cast<long double>(Denominator));
  }

  explicit operator float() const {
    return (static_cast<float>(Numerator) / static_cast<float>(Denominator));
  }

  void reduce() {
    T gcd = std::gcd(Numerator, Denominator);
    while (gcd > 1) {
      Numerator /= gcd;
      Denominator /= gcd;
      gcd = std::gcd(Numerator, Denominator);
    }
  }

  Fraction& operator+=(const Fraction& kX) {
    if (kX.Denominator == Denominator) {
      Numerator += kX.Numerator;
    } else {
      Numerator *= kX.Denominator;
      Numerator += kX.Numerator * Denominator;
      Denominator *= kX.Denominator;
    }
    return *this;
  }

  friend Fraction operator+(Fraction lhs, const Fraction& kRhs) {
    lhs += kRhs;
    return lhs;
  }

  Fraction& operator*=(const Fraction& kX) {
    Numerator *= kX.Numerator;
    Denominator *= kX.Denominator;
    return *this;
  }

  friend Fraction operator*(Fraction lhs, const Fraction& kRhs) {
    lhs *= kRhs;
    return lhs;
  }

  Fraction& operator/=(const Fraction& kX) {
    Numerator *= kX.Denominator;
    Denominator *= kX.Numerator;
    return *this;
  }

  friend Fraction operator/(Fraction lhs, const Fraction& kRhs) {
    lhs /= kRhs;
    return lhs;
  }
};  // Fraction

template <typename T>
std::ostream& operator<<(std::ostream& out, const Fraction<T>& kX) {
  out << kX.Numerator << " / " << kX.Denominator;
  return out;
}
#endif
