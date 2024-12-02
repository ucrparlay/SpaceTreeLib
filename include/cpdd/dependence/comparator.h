#pragma once
#include <algorithm>
#include <cfloat>
#include <cstdio>
#include <iostream>
#include <limits>

namespace cpdd {
//*----------- double precision comparision ----------------
//* God made the integers, all else is the work of man.
//* -- Leopold Kronecker
template <typename T>
class Num_Comparator {
 public:
  static inline bool Gt(T const a, T const b) {
    if constexpr (std::is_integral_v<T>) {
      return a > b;
    } else if (std::is_floating_point_v<T>) {
      return a - b > eps;
    } else {
      return a > b;
    }
  }

  static inline bool Lt(T const a, T const b) {
    if constexpr (std::is_integral_v<T>)
      return a < b;
    else if (std::is_floating_point_v<T>)
      return a - b < -eps;
    else
      return a < b;
  }

  static inline bool Eq(T const a, T const b) {
    if constexpr (std::is_integral_v<T>)
      return a == b;
    else if (std::is_floating_point_v<T>)
      return std::abs(a - b) < eps;
    else
      return a == b;
  }

  static inline bool Geq(T const a, T const b) {
    if constexpr (std::is_integral_v<T>)
      return a >= b;
    else if (std::is_floating_point_v<T>)
      return !Lt(a, b);
    else
      return a >= b;
  }

  static inline bool Leq(T const a, T const b) {
    if constexpr (std::is_integral_v<T>)
      return a <= b;
    else if (std::is_floating_point_v<T>)
      return !Gt(a, b);
    else
      return a <= b;
  }

  static inline T min(T const a, T const b) {
    if constexpr (std::is_integral_v<T>)
      return std::min(a, b);
    else if (std::is_floating_point_v<T>)
      return Lt(a, b) ? a : b;
    else
      return std::min(a, b);
  }

  static inline T max(T const a, T const b) {
    if constexpr (std::is_integral_v<T>)
      return std::max(a, b);
    else if (std::is_floating_point_v<T>)
      return Gt(a, b) ? a : b;
    else
      return std::max(a, b);
  }

 private:
  static constexpr T eps = std::numeric_limits<T>::epsilon();
};

}  // namespace cpdd
