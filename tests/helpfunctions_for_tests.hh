#ifndef HELPFUNCTIONS_FOR_TESTS_HH_
#define HELPFUNCTIONS_FOR_TESTS_HH_

#include <type_traits>
#include <vector>

template <class DataType>
double inline tol() {
  if constexpr (std::is_same<DataType, double>::value) {
    return 1e-8;
  }if constexpr (std::is_same<DataType, long double>::value) {
    return 1e-11;
  } else if constexpr (std::is_same<DataType, float>::value) {
    return 1e-1;
  }
  return 0;
}

template <typename T>
T horner_scheme(const std::vector<T>& kCoefficients, const T kX) {
  T res = 0;
  for (int i = kCoefficients.size() - 1; i >= 0; --i) {
    res = kCoefficients.at(i) + (kX * res);
  }
  return res;
}

template <typename T>
std::vector<T> compute_integral_polynoms(const std::vector<T>& kCoefficients1,
                            const std::vector<T>& kCoefficients2) {
    std::vector<T> res(kCoefficients1.size() + kCoefficients2.size(), 0);
    for (int i = 0; i < kCoefficients1.size(); ++i) {
      for (int ii = 0; ii < kCoefficients2.size(); ++ii) {
        res.at(i + ii + 1) += kCoefficients1.at(i) * kCoefficients2.at(ii);
      }
    }
    for (int i = 1; i < res.size(); ++i) {
      res.at(i) /= i;
    }
    return res;
}

#endif
