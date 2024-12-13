#ifndef PSTP_DEPENDENCE_COMPARATOR_H_
#define PSTP_DEPENDENCE_COMPARATOR_H_

#include <algorithm>
#include <cfloat>
#include <cstdio>
#include <iostream>
#include <limits>

namespace pstp {
//*----------- double precision comparision ----------------
//* God made the integers, all else is the work of man.
//* -- Leopold Kronecker
template <typename T>
class Num_Comparator {
 public:
  static inline bool Gt(T const a, T const b) {
    if constexpr (std::is_floating_point_v<T>) {
      return a - b > eps;
    } else {
      return a > b;
    }
  }

  static inline bool Lt(T const a, T const b) {
    if constexpr (std::is_floating_point_v<T>) {
      return a - b < -eps;
    } else {
      return a < b;
    }
  }

  static inline bool Eq(T const a, T const b) {
    if constexpr (std::is_floating_point_v<T>) {
      return std::abs(a - b) < eps;
    } else {
      return a == b;
    }
  }

  static inline bool Geq(T const a, T const b) { return !Lt(a, b); }

  static inline bool Leq(T const a, T const b) { return !Gt(a, b); }

  static inline T min(T const a, T const b) { return Lt(a, b) ? a : b; }

  static inline T max(T const a, T const b) { return Gt(a, b) ? a : b; }

 private:
  static constexpr T eps = std::numeric_limits<T>::epsilon();
};

}  // namespace pstp

#endif  // PSTP_DEPENDENCE_COMPARATOR_H_
