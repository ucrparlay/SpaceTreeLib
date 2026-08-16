#ifndef PSI_BASE_TREE_IMPL_CIRCLE_OP_HPP_
#define PSI_BASE_TREE_IMPL_CIRCLE_OP_HPP_

#include <cmath>

#include "../base_tree.h"
#include "psi/dependence/concepts.h"

namespace psi
{
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename CircleType>
inline bool base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::legal_circle(
	CircleType const &cl)
{
	return num_type::geq(cl.get_radius(), 0);
}

// NOTE: point within the circle
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename CircleType>
inline bool base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::within_circle(
	Point const &p, CircleType const &cl)
{
	coord_type r = 0;
	for (dims_type i = 0; i < num_dims; ++i) {
		r += (p.pnt[i] - cl.get_center().pnt[i]) *
		     (p.pnt[i] - cl.get_center().pnt[i]);
		if (num_type::gt(r, cl.get_radius() * cl.get_radius()))
			return false;
	}
	assert(num_type::leq(r, cl.get_radius() * cl.get_radius()));
	return true;
}

// NOTE: box within the circle
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename CircleType>
inline bool base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::within_circle(
	box_type const &box, CircleType const &cl)
{
	assert(legal_box(box));
	assert(legal_circle(cl));

	// NOTE: the logical is same as p2b_max_distance <= radius
	coord_type r = 0;
	for (dims_type i = 0; i < num_dims; ++i) {
		if (num_type::lt(cl.get_center().pnt[i],
				 (box.first.pnt[i] + box.second.pnt[i]) / 2)) {
			r += (box.second.pnt[i] - cl.get_center().pnt[i]) *
			     (box.second.pnt[i] - cl.get_center().pnt[i]);
		} else {
			r += (cl.get_center().pnt[i] - box.first.pnt[i]) *
			     (cl.get_center().pnt[i] - box.first.pnt[i]);
		}
		if (num_type::gt(r, cl.get_radius_square()))
			return false;
	}
	assert(num_type::leq(r, cl.get_radius_square()));
	return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename CircleType>
inline bool
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::circle_intersect_box(
	CircleType const &cl, box_type const &box)
{
	assert(legal_box(box));
	assert(legal_circle(cl));
	// NOTE: the logical is same as p2b_min_distance > radius
	coord_type r = 0;
	for (dims_type i = 0; i < num_dims; ++i) {
		if (num_type::lt(cl.get_center().pnt[i], box.first.pnt[i])) {
			r += (box.first.pnt[i] - cl.get_center().pnt[i]) *
			     (box.first.pnt[i] - cl.get_center().pnt[i]);
		} else if (num_type::gt(cl.get_center().pnt[i],
					box.second.pnt[i])) {
			r += (cl.get_center().pnt[i] - box.second.pnt[i]) *
			     (cl.get_center().pnt[i] - box.second.pnt[i]);
		}
		if (num_type::gt(r, cl.get_radius() * cl.get_radius()))
			return false; //* not intersect
	}
	assert(num_type::leq(r, cl.get_radius() * cl.get_radius()));
	return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename CircleType>
inline CircleType
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_circle(
	box_type const &box)
{
	assert(legal_box(box) && box != get_empty_box());

	coord_type r = 0;
	if constexpr (num_dims == 2) {
		r = num_type::divide_two_ceil(box.second.pnt[0] -
					      box.first.pnt[0]) *
			    num_type::divide_two_ceil(box.second.pnt[0] -
						      box.first.pnt[0]) +
		    num_type::divide_two_ceil(box.second.pnt[1] -
					      box.first.pnt[1]) *
			    num_type::divide_two_ceil(box.second.pnt[1] -
						      box.first.pnt[1]);
	} else if constexpr (num_dims == 3) {
		r = num_type::divide_two_ceil(box.second.pnt[0] -
					      box.first.pnt[0]) *
			    num_type::divide_two_ceil(box.second.pnt[0] -
						      box.first.pnt[0]) +
		    num_type::divide_two_ceil(box.second.pnt[1] -
					      box.first.pnt[1]) *
			    num_type::divide_two_ceil(box.second.pnt[1] -
						      box.first.pnt[1]) +
		    num_type::divide_two_ceil(box.second.pnt[2] -
					      box.first.pnt[2]) *
			    num_type::divide_two_ceil(box.second.pnt[2] -
						      box.first.pnt[2]);
	} else {
		for (dims_type i = 0; i < num_dims; ++i) {
			r += num_type::divide_two_ceil(box.second.pnt[i] -
						       box.first.pnt[i]) *
			     num_type::divide_two_ceil(box.second.pnt[i] -
						       box.first.pnt[i]);
		}
	}

	// std::cout << "r: " << std::sqrt(r) << '\n';
	// assert(within_circle(
	//     box,
	//     CircleType{get_box_center(box),
	//     CircleType::compute_radius(std::sqrt(r))}));

	return CircleType{get_box_center(box),
			  CircleType::compute_radius(std::sqrt(r))};
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename CircleType>
inline CircleType
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_circle(slice_type V)
{
	return get_circle<CircleType>(get_box(V));
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename CircleType>
inline CircleType
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_circle(
	parlay::sequence<CircleType> const &circle_seq)
{
	CircleType circle = circle_seq[0];
	for (size_t i = 1; i < circle_seq.size(); i++) {
		circle = get_circle(circle, circle_seq[i]);
	}
	return circle;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename CircleType1, typename CircleType2>
inline bool
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::circle_within_circle(
	CircleType1 const &a, CircleType2 const &b)
{
	assert(legal_circle(a) && legal_circle(b));
	return num_type::leq(a.get_radius(), b.get_radius()) &&
	       num_type::leq(
		       p2p_distance_square(a.get_center(), b.get_center()),
		       (b.get_radius() - a.get_radius()) *
			       (b.get_radius() - a.get_radius()));
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename CircleType>
inline CircleType
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_circle(
	CircleType const &a, CircleType const &b)
{
	if (circle_within_circle(a, b)) {
		return b;
	} else if (circle_within_circle(b, a)) {
		return a;
	}

	Point center;
	for (dims_type i = 0; i < num_dims; ++i) {
		center.pnt[i] =
			(a.get_center().pnt[i] + b.get_center().pnt[i]) / 2;
	}
	double radius =
		std::sqrt(p2p_distance_square(a.get_center(), b.get_center())) /
			2.0 +
		std::max(a.get_radius(), b.get_radius());

	return CircleType{center, CircleType::compute_radius(radius)};
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename CircleType>
inline CircleType
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_circle(
	Point const &p, CircleType const &cl)
{
	if (point_within_circle(p, cl)) {
		return cl;
	}
	coord_type radius = p2p_distance_square(p, cl.get_center());
	return CircleType{cl.get_center(),
			  CircleType::compute_radius(std::sqrt(radius))};
}

} // namespace psi

#endif // PSI_BASE_TREE_IMPL_CIRCLE_OP_HPP_
