#ifndef EXPLICIT_LEGENDRE_FUNCTIONS_HH_
#define EXPLICIT_LEGENDRE_FUNCTIONS_HH_

#include <cassert>
#include <cmath>
#include <type_traits>
#include <vector>

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 0, Vtype> legendre(const Vtype kX) {
  return (Vtype{1});
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 1, Vtype> legendre(const Vtype kX) {
  return (kX);
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 2, Vtype> legendre(const Vtype kX) {
  return (1.0 / 2.0 * (3 * std::pow(kX, 2) - 1.0));
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 3, Vtype> legendre(const Vtype kX) {
  return (1.0 / 2.0 * (5 * std::pow(kX, 3) - 3.0 * kX));
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 4, Vtype> legendre(const Vtype kX) {
  return (1.0 / 8.0 * (35.0 * std::pow(kX, 4) - 30.0 * std::pow(kX, 2) + 3.0));
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 5, Vtype> legendre(const Vtype kX) {
  return (1.0 / 8.0 *
          (63.0 * std::pow(kX, 5) - 70.0 * std::pow(kX, 3) + 15.0 * kX));
}

template <typename Vtype = double>
Vtype legendre(const int Degree, const Vtype kX) {
  switch (Degree) {
    case 0:
      return legendre<0>(kX);
    case 1:
      return legendre<1>(kX);
    case 2:
      return legendre<2>(kX);
    case 3:
      return legendre<3>(kX);
    case 4:
      return legendre<4>(kX);
    case 5:
      return legendre<5>(kX);
    default:
      assert(false);
  }
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 0, std::vector<Vtype>> legendre_polynome() {
  std::vector<Vtype> res = {1};
  return res;
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 1, std::vector<Vtype>> legendre_polynome() {
  std::vector<Vtype> res = {0.0, 1.0};
  return res;
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 2, std::vector<Vtype>> legendre_polynome() {
  std::vector<Vtype> res = {-0.5, 0, 1.5};
  return res;
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 3, std::vector<Vtype>> legendre_polynome() {
  std::vector<Vtype> res = {0, -1.5, 0, 2.5};
  return res;
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 4, std::vector<Vtype>> legendre_polynome() {
  std::vector<Vtype> res = {0.375, 0, -3.75, 0, 4.375};
  return res;
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 5, std::vector<Vtype>> legendre_polynome() {
  std::vector<Vtype> res = {0, 1.875, 0, -8.75, 0, 7.875};
  return res;
}

template <typename Vtype = double>
std::vector<Vtype> legendre_polynome(const int Degree) {
  switch (Degree) {
    case 0:
      return legendre_polynome<0, Vtype>();
    case 1:
      return legendre_polynome<1, Vtype>();
    case 2:
      return legendre_polynome<2, Vtype>();
    case 3:
      return legendre_polynome<3, Vtype>();
    case 4:
      return legendre_polynome<4, Vtype>();
    case 5:
      return legendre_polynome<5, Vtype>();
    default:
      assert(false);
  }
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 0, Vtype> legendre_integral(Vtype start, Vtype end) {
  return end - start;
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 1, Vtype> legendre_integral(Vtype start, Vtype end) {
  auto f = [](Vtype x) { return 0.5 * x * x; };
  return f(end) - f(start);
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 2, Vtype> legendre_integral(Vtype start, Vtype end) {
  auto f = [](Vtype x) { return 0.5 * (std::pow(x, 3) - x); };
  return f(end) - f(start);
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 3, Vtype> legendre_integral(Vtype start, Vtype end) {
  auto f = [](Vtype x) {
    return 0.625 * std::pow(x, 4) - 0.75 * std::pow(x, 2);
  };
  return f(end) - f(start);
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 4, Vtype> legendre_integral(Vtype start, Vtype end) {
  auto f = [](Vtype x) {
    return (1.0 / 8.0) * (7 * std::pow(x, 5) - 10.0 * std::pow(x, 3) + 3 * x);
  };
  return f(end) - f(start);
}

template <int Degree, typename Vtype = double>
std::enable_if_t<Degree == 5, Vtype> legendre_integral(Vtype start, Vtype end) {
  auto f = [](Vtype x) {
    return 1.3125 * std::pow(x, 6) - 2.1875 * std::pow(x, 4) +
           0.9375 * std::pow(x, 2);
  };
  return f(end) - f(start);
}

template <typename Vtype = double>
Vtype legendre_integral(const int Degree, Vtype start, Vtype end) {
  switch (Degree) {
    case 0:
      return legendre_integral<0, Vtype>(start, end);
    case 1:
      return legendre_integral<1, Vtype>(start, end);
    case 2:
      return legendre_integral<2, Vtype>(start, end);
    case 3:
      return legendre_integral<3, Vtype>(start, end);
    case 4:
      return legendre_integral<4, Vtype>(start, end);
    case 5:
      return legendre_integral<5, Vtype>(start, end);
    default:
      assert(false);
  }
}

#endif
