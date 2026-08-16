#ifndef PSI_TESTS_UNIT_REFERENCE_H_
#define PSI_TESTS_UNIT_REFERENCE_H_

/*
 * Deterministic inputs and brute-force answers. The tests compare the tree
 * against these, so nothing here may use the tree.
 */

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

#include "psi/dependence/basic_point.h"

namespace psi_test
{

/* Fixed seed by default: a failing run has to be replayable. */
inline uint64_t seed_from_env()
{
	char const *s = std::getenv("PSI_TEST_SEED");
	return s ? std::strtoull(s, nullptr, 10) : 20260816u;
}

template <typename Point>
std::vector<Point> random_points(size_t n, typename Point::coord_type lo,
				 typename Point::coord_type hi, uint64_t seed)
{
	using coord_type = typename Point::coord_type;
	std::mt19937_64 gen(seed);
	std::vector<Point> v(n);
	for (size_t i = 0; i < n; i++) {
		coord_type c[Point::get_dim()];
		for (size_t d = 0; d < Point::get_dim(); d++) {
			if constexpr (std::is_integral_v<coord_type>) {
				c[d] = std::uniform_int_distribution<
					coord_type>(lo, hi)(gen);
			} else {
				c[d] = std::uniform_real_distribution<
					coord_type>(lo, hi)(gen);
			}
		}
		v[i] = Point(c);
	}
	return v;
}

/* Every point identical, or every point sharing one coordinate: the shapes
 * that used to send the split rules in circles. */
template <typename Point>
std::vector<Point> repeated_points(size_t n, typename Point::coord_type value)
{
	typename Point::coord_type c[Point::get_dim()];
	for (size_t d = 0; d < Point::get_dim(); d++)
		c[d] = value;
	return std::vector<Point>(n, Point(c));
}

template <typename Point, typename Box>
bool within(Point const &p, Box const &box)
{
	for (size_t d = 0; d < Point::get_dim(); d++) {
		if (p.pnt[d] < box.first.pnt[d] || p.pnt[d] > box.second.pnt[d])
			return false;
	}
	return true;
}

template <typename Point, typename Box>
size_t range_count(std::vector<Point> const &pts, Box const &box)
{
	size_t n = 0;
	for (auto const &p : pts)
		n += within(p, box) ? 1 : 0;
	return n;
}

template <typename Point, typename Box>
std::vector<Point> range_query(std::vector<Point> const &pts, Box const &box)
{
	std::vector<Point> out;
	for (auto const &p : pts)
		if (within(p, box))
			out.push_back(p);
	std::sort(out.begin(), out.end());
	return out;
}

template <typename Point>
auto distance_square(Point const &a, Point const &b)
{
	typename Point::dis_type s = 0;
	for (size_t d = 0; d < Point::get_dim(); d++) {
		auto diff = static_cast<typename Point::dis_type>(a.pnt[d]) -
			    static_cast<typename Point::dis_type>(b.pnt[d]);
		s += diff * diff;
	}
	return s;
}

/* The k nearest squared distances, ascending. Ties make the identity of the
 * k-th neighbour ambiguous, but its distance is not. */
template <typename Point>
std::vector<typename Point::dis_type>
knn_distances(std::vector<Point> const &pts, Point const &q, size_t k)
{
	std::vector<typename Point::dis_type> d;
	d.reserve(pts.size());
	for (auto const &p : pts)
		d.push_back(distance_square(p, q));
	std::sort(d.begin(), d.end());
	d.resize(std::min(k, d.size()));
	return d;
}

/*
 * Distances are compared, not points, and a floating-point distance summed in
 * a different order than the tree sums it differs in the last ulp.
 */
template <typename T>
bool near(T a, T b)
{
	if constexpr (std::is_integral_v<T>) {
		return a == b;
	} else {
		T const d = a > b ? a - b : b - a;
		T const m = std::max(std::abs(a), std::abs(b));
		return d <= (m > T(0) ? m : T(1)) * T(1e-12);
	}
}

template <typename T>
bool near(std::vector<T> const &a, std::vector<T> const &b)
{
	if (a.size() != b.size())
		return false;
	for (size_t i = 0; i < a.size(); i++)
		if (!near(a[i], b[i]))
			return false;
	return true;
}

template <typename Point>
std::vector<Point> sorted(std::vector<Point> v)
{
	std::sort(v.begin(), v.end());
	return v;
}

} // namespace psi_test

#endif /* PSI_TESTS_UNIT_REFERENCE_H_ */
