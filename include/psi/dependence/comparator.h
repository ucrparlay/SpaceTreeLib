#ifndef PSI_DEPENDENCE_COMPARATOR_H_
#define PSI_DEPENDENCE_COMPARATOR_H_

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>

namespace psi
{
/*
 * *----------- double precision comparision ----------------
 * * God made the integers, all else is the work of man.
 * * -- Leopold Kronecker
 */
template <typename T>
class num_comparator
{
public:
	static inline bool gt(T const a, T const b)
	{
		if constexpr (std::is_floating_point_v<T>) {
			return a - b > eps;
		} else {
			return a > b;
		}
	}

	static inline bool lt(T const a, T const b)
	{
		if constexpr (std::is_floating_point_v<T>) {
			return a - b < -eps;
		} else {
			return a < b;
		}
	}

	static inline bool eq(T const a, T const b)
	{
		if constexpr (std::is_floating_point_v<T>) {
			return std::abs(a - b) < eps;
		} else {
			return a == b;
		}
	}

	static inline bool geq(T const a, T const b)
	{
		return !lt(a, b);
	}

	static inline bool leq(T const a, T const b)
	{
		return !gt(a, b);
	}

	static inline T min(T const a, T const b)
	{
		return lt(a, b) ? a : b;
	}

	static inline T max(T const a, T const b)
	{
		return gt(a, b) ? a : b;
	}

	static inline T abs(T const a)
	{
		return lt(a, 0) ? -a : a;
	}

	static inline bool is_zero(T const a)
	{
		return eq(a, static_cast<T>(0));
	}

	static inline T divide_two_ceil(T const a)
	{
		if constexpr (std::is_floating_point_v<T>) {
			return std::ceil(a / 2.0);
		} else {
			return (a + 1) / 2;
		}
	}

private:
	static constexpr T eps = std::numeric_limits<T>::epsilon();
};

} /* namespace psi */

#endif /* PSI_DEPENDENCE_COMPARATOR_H_ */
