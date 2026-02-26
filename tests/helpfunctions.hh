#ifndef HELPFUNCTIONS_HH_
#define HELPFUNCTIONS_HH_

// TODO (Georg): maybe declaring inline functions static improves performance?

#include <cassert>
#include <iostream>
// #include <type_traits>
#include <vector>

template <class C>
std::ostream& operator<<(std::ostream& a_out, const std::vector<C>& ak_v) {
  a_out << "{";
  for (auto iter = ak_v.begin(); iter != ak_v.end(); ++iter) {
    a_out << *iter;
    if (iter + 1 != ak_v.end()) {
      a_out << ", ";
    }
  }
  a_out << "}";
  return a_out;
}

#endif
