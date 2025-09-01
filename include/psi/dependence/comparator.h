#ifndef PSI_DEPENDENCE_COMPARATOR_H_
#define PSI_DEPENDENCE_COMPARATOR_H_

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>

namespace psi {
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

  static inline T Min(T const a, T const b) { return Lt(a, b) ? a : b; }

  static inline T Max(T const a, T const b) { return Gt(a, b) ? a : b; }

  static inline T Abs(T const a) { return Lt(a, 0) ? -a : a; }

  static inline bool IsZero(T const a) { return Eq(a, static_cast<T>(0)); }

  static inline T DivideTwoCeil(T const a) {
    if constexpr (std::is_floating_point_v<T>) {
      return std::ceil(a / 2.0);
    } else {
      return (a + 1) / 2;
    }
  }

  static inline T integer_log2_upper(T const a) {
    if constexpr (std::is_floating_point_v<T>) {
      return static_cast<T>(std::ceil(std::log2(a)));
    } else {
      if (a == 0) return 0;  // BUG: this is a bug, should return -1
      T r = 0, b = a;
      while (b >>= 1) r++;
      return r + static_cast<T>((1 << r) != a);
    }
  }

 private:
  static constexpr T eps = std::numeric_limits<T>::epsilon();
};

}  // namespace psi

#endif  // PSI_DEPENDENCE_COMPARATOR_H_
