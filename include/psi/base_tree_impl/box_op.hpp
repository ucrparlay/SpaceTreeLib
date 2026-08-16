#ifndef PSI_BASE_TREE_IMPL_BOX_OP_HPP_
#define PSI_BASE_TREE_IMPL_BOX_OP_HPP_

#include "psi/base_tree.h"
#include "psi/dependence/concepts.h"

namespace psi
{

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::coord_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_box_mid(
	dims_type const d, box_type const &box)
{
	if constexpr (std::is_integral_v<coord_type>) {
		/*
		 * Points on the split plane go right, so the midpoint rounds
		 * up: (1, 1) and (1, 2) only separate if the split is at 2.
		 *
		 * first + ceil((second - first) / 2), not ceil((first +
		 * second) / 2): the sum overflows once the two bounds reach
		 * 2^63, which put the plane outside the box and made
		 * spatial_median recurse forever.
		 */
		return box.first.pnt[d] +
		       num_type::divide_two_ceil(static_cast<coord_type>(
			       box.second.pnt[d] - box.first.pnt[d]));
	} else {
		return (box.first.pnt[d] + box.second.pnt[d]) / 2;
	}
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::legal_box(
	box_type const &bx)
{
	// TODO: remove it
	if (bx == get_empty_box())
		return true;

	if constexpr (num_dims == 2) {
		return !num_type::gt(bx.first.pnt[0], bx.second.pnt[0]) &&
		       !num_type::gt(bx.first.pnt[1], bx.second.pnt[1]);
	} else if constexpr (num_dims == 3) {
		return !num_type::gt(bx.first.pnt[0], bx.second.pnt[0]) &&
		       !num_type::gt(bx.first.pnt[1], bx.second.pnt[1]) &&
		       !num_type::gt(bx.first.pnt[2], bx.second.pnt[2]);
	} else {
		for (dims_type i = 0; i < num_dims; ++i) {
			if (num_type::gt(bx.first.pnt[i], bx.second.pnt[i])) {
				return false;
			}
		}
		return true;
	}
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::within_box(
	box_type const &a, box_type const &b)
{
	constexpr uint_fast8_t num_dims = Point::get_dim();
	assert(legal_box(a));
	assert(legal_box(b));

	if constexpr (num_dims == 2) {
		return !num_type::lt(a.first.pnt[0], b.first.pnt[0]) &&
		       !num_type::lt(a.first.pnt[1], b.first.pnt[1]) &&
		       !num_type::gt(a.second.pnt[0], b.second.pnt[0]) &&
		       !num_type::gt(a.second.pnt[1], b.second.pnt[1]);
	} else if constexpr (num_dims == 3) {
		return !num_type::lt(a.first.pnt[0], b.first.pnt[0]) &&
		       !num_type::lt(a.first.pnt[1], b.first.pnt[1]) &&
		       !num_type::lt(a.first.pnt[2], b.first.pnt[2]) &&
		       !num_type::gt(a.second.pnt[0], b.second.pnt[0]) &&
		       !num_type::gt(a.second.pnt[1], b.second.pnt[1]) &&
		       !num_type::gt(a.second.pnt[2], b.second.pnt[2]);
	} else {
		for (dims_type i = 0; i < num_dims; ++i) {
			if (num_type::lt(a.first.pnt[i], b.first.pnt[i]) ||
			    num_type::gt(a.second.pnt[i], b.second.pnt[i])) {
				return false;
			}
		}
		return true;
	}
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::same_box(box_type const &a,
							     box_type const &b)
{
	assert(legal_box(a));
	assert(legal_box(b));

	for (dims_type i = 0; i < num_dims; ++i) {
		if (!num_type::eq(a.first.pnt[i], b.first.pnt[i]) ||
		    !num_type::eq(a.second.pnt[i], b.second.pnt[i])) {
			return false;
		}
	}
	return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::within_box(
	Point const &p, box_type const &bx)
{
	assert(legal_box(bx));

	if constexpr (num_dims == 2) {
		return !num_type::lt(p.pnt[0], bx.first.pnt[0]) &&
		       !num_type::lt(p.pnt[1], bx.first.pnt[1]) &&
		       !num_type::gt(p.pnt[0], bx.second.pnt[0]) &&
		       !num_type::gt(p.pnt[1], bx.second.pnt[1]);
	} else if constexpr (num_dims == 3) {
		return !num_type::lt(p.pnt[0], bx.first.pnt[0]) &&
		       !num_type::lt(p.pnt[1], bx.first.pnt[1]) &&
		       !num_type::lt(p.pnt[2], bx.first.pnt[2]) &&
		       !num_type::gt(p.pnt[0], bx.second.pnt[0]) &&
		       !num_type::gt(p.pnt[1], bx.second.pnt[1]) &&
		       !num_type::gt(p.pnt[2], bx.second.pnt[2]);
	} else {
		for (dims_type i = 0; i < num_dims; ++i) {
			if (num_type::lt(p.pnt[i], bx.first.pnt[i]) ||
			    num_type::gt(p.pnt[i], bx.second.pnt[i])) {
				return false;
			}
		}
		return true;
	}
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::box_intersect_box(
	box_type const &a, box_type const &b)
{
	assert(legal_box(a) && legal_box(b));

	if constexpr (num_dims == 2) {
		return !(num_type::lt(a.second.pnt[0], b.first.pnt[0]) ||
			 num_type::gt(a.first.pnt[0], b.second.pnt[0]) ||
			 num_type::lt(a.second.pnt[1], b.first.pnt[1]) ||
			 num_type::gt(a.first.pnt[1], b.second.pnt[1]));
	} else if constexpr (num_dims == 3) {
		return !(num_type::lt(a.second.pnt[0], b.first.pnt[0]) ||
			 num_type::gt(a.first.pnt[0], b.second.pnt[0]) ||
			 num_type::lt(a.second.pnt[1], b.first.pnt[1]) ||
			 num_type::gt(a.first.pnt[1], b.second.pnt[1]) ||
			 num_type::lt(a.second.pnt[2], b.first.pnt[2]) ||
			 num_type::gt(a.first.pnt[2], b.second.pnt[2]));
	} else {
		for (dims_type i = 0; i < num_dims; ++i) {
			if (num_type::lt(a.second.pnt[i], b.first.pnt[i]) ||
			    num_type::gt(a.first.pnt[i], b.second.pnt[i])) {
				return false;
			}
		}
		return true;
	}
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::is_box_line_in_dimension(
	box_type const &box, dims_type d)
{
	assert(legal_box(box));
	return num_type::eq(box.first[d], box.second[d]);
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::vertical_line_split_box(
	coord_type const &line, box_type const &box, dims_type d)
{
	assert(legal_box(box));
	return vertical_line_intersect_box_exclude(line, box, d) ||
	       (vertical_line_on_box_right_edge(line, box, d) &&
		!is_box_line_in_dimension(box, d));
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool
base_tree<Point, DerivedTree, SkHeight,
	  ImbaRatio>::vertical_line_on_box_left_edge(coord_type const &line,
						     box_type const &box,
						     dims_type d)
{
	assert(legal_box(box));
	return num_type::eq(line, box.first.pnt[d]);
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool
base_tree<Point, DerivedTree, SkHeight,
	  ImbaRatio>::vertical_line_on_box_right_edge(coord_type const &line,
						      box_type const &box,
						      dims_type d)
{
	assert(legal_box(box));
	return num_type::eq(line, box.second.pnt[d]);
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::vertical_line_on_box_edge(
	coord_type const &line, box_type const &box, dims_type d)
{
	assert(legal_box(box));
	return vertical_line_on_box_left_edge(line, box, d) ||
	       vertical_line_on_box_right_edge(line, box, d);
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::vertical_line_intersect_box(
	coord_type const &line, box_type const &box, dims_type d)
{
	assert(legal_box(box));
	return num_type::geq(line, box.first.pnt[d]) &&
	       num_type::leq(line, box.second.pnt[d]);
}

// NOTE: if the line @line is one the boundary of the box, then it will be
// considered as not intersect
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline bool base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::
	vertical_line_intersect_box_exclude(coord_type const &line,
					    box_type const &box, dims_type d)
{
	assert(legal_box(box));
	return num_type::gt(line, box.first.pnt[d]) &&
	       num_type::lt(line, box.second.pnt[d]);
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
inline typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::box_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_empty_box()
{
	return box_type(Point(std::numeric_limits<coord_type>::max()),
			Point(std::numeric_limits<coord_type>::lowest()));
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::box_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_box(box_type const &x,
							    box_type const &y)
{
	return box_type(x.first.min_coords(y.first),
			x.second.max_coords(y.second));
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::box_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_box(slice_type V)
{
	if (V.size() == 0) {
		return get_empty_box();
	} else {
		auto minmax = [&](box_type const &x, box_type const &y) {
			return box_type(x.first.min_coords(y.first),
					x.second.max_coords(y.second));
		};
		auto boxes =
			parlay::delayed_seq<box_type>(V.size(), [&](size_t i) {
				return box_type(V[i].get_coords(),
						V[i].get_coords());
			});
		return parlay::reduce(boxes,
				      parlay::make_monoid(minmax, boxes[0]));
	}
}

// NOTE: this function omit the possibility that T contains the bounding box --
// it will always try to reduce a bounding box in the tree T
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::box_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_box(node *T)
{
	if (T->size == 0) {
		return get_empty_box();
	}
	if (T->is_leaf) {
		leaf_type *tl = static_cast<leaf_type *>(T);
		if (tl->is_dummy) {
			return box_type(tl->pts[0], tl->pts[0]);
		}
		return get_box(tl->pts.cut(0, tl->size));
		// auto nb = get_box(tl->pts.cut(0, tl->size));
		// assert(same_box(nb, tl->get_box()));
		// return nb;
	}
	interior_type *ti = static_cast<interior_type *>(T);
	if constexpr (is_binary_node<interior_type>) {
		box_type const &left_box =
			get_box<leaf_type, interior_type>(ti->left);
		box_type const &right_box =
			get_box<leaf_type, interior_type>(ti->right);
		return get_box(left_box, right_box);
		// auto nb = get_box(left_box, right_box);
		// assert(same_box(nb, ti->get_box()));
		// return nb;
	} else if constexpr (is_multi_node<interior_type>) {
		box_seq_type return_box_seq(interior_type::get_regions());
		for (size_t i = 0; i < interior_type::get_regions(); i++) {
			return_box_seq[i] = get_box<leaf_type, interior_type>(
				ti->tree_nodes[i]);
		}
		return get_box(return_box_seq);
	} else {
		static_assert(is_binary_node<interior_type> ||
				      is_multi_node<interior_type>,
			      "get_box supports only binary and multi-way "
			      "interior nodes");
		return get_empty_box();
	}
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::box_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_box(
	box_seq_type const &box_seq)
{
	box_type box = get_empty_box();
	for (auto const &b : box_seq) {
		box = get_box(box, b);
	}
	return std::move(box);
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
Point base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::get_box_center(
	box_type const &box)
{
	Point center;
	if constexpr (num_dims == 2) {
		center.pnt[0] = (box.first.pnt[0] + box.second.pnt[0]) / 2;
		center.pnt[1] = (box.first.pnt[1] + box.second.pnt[1]) / 2;
	} else if constexpr (num_dims == 3) {
		center.pnt[0] = (box.first.pnt[0] + box.second.pnt[0]) / 2;
		center.pnt[1] = (box.first.pnt[1] + box.second.pnt[1]) / 2;
		center.pnt[2] = (box.first.pnt[2] + box.second.pnt[2]) / 2;
	} else {
		for (dims_type i = 0; i < num_dims; ++i) {
			center.pnt[i] = get_box_mid(i, box);
		}
	}
	return center;
}
} // namespace psi

#endif // PSI_BASE_TREE_IMPL_BOX_OP_HPP_
