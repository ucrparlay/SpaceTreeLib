#ifndef PSI_BASE_TREE_IMPL_LEAF_OP_HPP_
#define PSI_BASE_TREE_IMPL_LEAF_OP_HPP_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#include "psi/base_tree.h"
#include "psi/dependence/concepts.h"

namespace psi
{

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename Range>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::extract_points_in_leaf(
	node *T, Range out)
{
	leaf_type *tl = static_cast<leaf_type *>(T);
	if (tl->is_dummy) {
		std::ranges::fill_n(out.begin(), tl->size, tl->pts[0]);
	} else {
		std::ranges::copy(tl->pts.begin(), tl->pts.begin() + tl->size,
				  out.begin());
	}
	return;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type>
node *base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::insert_points2_leaf(
	node *T, slice_type in)
{
	leaf_type *tl = static_cast<leaf_type *>(T);
	if (tl->is_dummy) {
		T->size += in.size();
		return T;
	}

	assert(T->size + in.size() <= leaf_capacity);
	if (tl->pts.size() == 0) {
		assert(tl->size == 0);
		tl->pts = points_type::uninitialized(leaf_capacity);
	}
	for (size_t i = 0; i < in.size(); i++) {
		tl->pts[tl->size + i] = in[i];
	}
	tl->size += in.size();

	tl->update_aug(tl->pts.cut(0, tl->size));
	return T;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename ReturnType>
ReturnType
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::delete_points4_leaf(
	node *T, slice_type in)
{
	assert(T->size >= in.size());
	leaf_type *tl = static_cast<leaf_type *>(T);

	if (tl->is_dummy) {
		assert(in.size() <= T->size);
		tl->size -= in.size(); // WARN: this assumes that in\in T
		if (tl->size == 0) {
			tl->is_dummy = false;
			tl->pts = points_type::uninitialized(leaf_capacity);
			tl->reset_aug();
		}

		if constexpr (std::same_as<ReturnType, node *>) {
			return T;
		} else if constexpr (std::same_as<ReturnType, node_box_type>) {
			if constexpr (has_box<typename leaf_type::at_type>) {
				return node_box_type(T, tl->get_box());
			} else {
				return node_box_type(
					T, T->size ? box_type(tl->pts[0],
							      tl->pts[0])
						   : get_empty_box());
			}
		} else {
			;
		}
	}

	auto it = tl->pts.begin(), end = tl->pts.begin() + tl->size;
	for (size_t i = 0; i < in.size(); i++) {
		it = std::ranges::find(tl->pts.begin(), end, in[i]);
		assert(it != end);
		std::ranges::iter_swap(it, --end);
	}

	assert(std::cmp_equal(std::distance(tl->pts.begin(), end),
			      tl->size - in.size()));
	tl->size -= in.size();
	assert(tl->size >= 0);
	tl->update_aug(tl->pts.cut(0, tl->size));
	// assert(tl->get_box() == get_box(tl->pts.cut(0, tl->size)));

	if constexpr (std::same_as<ReturnType, node *>) {
		return T;
	} else if constexpr (std::same_as<ReturnType, node_box_type>) {
		if constexpr (has_box<typename leaf_type::at_type>) {
			return node_box_type(T, tl->get_box());
			// return node_box_type(T, get_box(tl->pts.cut(0,
			// tl->size)));
		} else {
			return node_box_type(T,
					     get_box(tl->pts.cut(0, tl->size)));
		}
	} else {
	}
}

// NOTE: diff points from the leaf using std::set_difference
// {1, 2, 5, 5, 5, 9} ∖ {2, 5, 7} == {1, 5, 5, 9}
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename ReturnType>
ReturnType
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::diff_points4_leaf(
	node *T, slice_type in)
{
	leaf_type *tl = static_cast<leaf_type *>(T);

	if (tl->is_dummy) {
		size_t cnt = parlay::count(in, tl->pts[0]);
		tl->size = tl->size >= cnt ? tl->size - cnt : 0;
		assert(tl->size >= 0);
		if (tl->size == 0) { // set points to normal leaf when all
				     // points in dummy node has been deleted
			tl->is_dummy = false;
			tl->pts = points_type::uninitialized(leaf_capacity);
			tl->reset_aug();
		}

		if constexpr (std::same_as<ReturnType, node *>) {
			return T;
		} else if constexpr (std::same_as<ReturnType, node_box_type>) {
			if constexpr (has_box<typename leaf_type::at_type>) {
				return node_box_type(T, tl->get_box());
			} else {
				return node_box_type(
					T, T->size ? box_type(tl->pts[0],
							      tl->pts[0])
						   : get_empty_box());
			}
		} else {
			;
		}
	}

	// NOTE: for normal leaf, need to check whether all points_type are in
	// the leaf
	auto diff_res = std::ranges::set_difference(
		parlay::sort(tl->pts.cut(0, tl->size)), parlay::sort(in),
		tl->pts.begin(),
		[](Point const &p1, Point const &p2) { return p1 < p2; });
	tl->size = std::ranges::distance(tl->pts.begin(), diff_res.out);
	tl->update_aug(tl->pts.cut(0, tl->size));

	if constexpr (std::same_as<ReturnType, node *>) {
		return T;
	} else if constexpr (std::same_as<ReturnType, node_box_type>) {
		if constexpr (has_box<typename leaf_type::at_type>) {
			return node_box_type(T, tl->get_box());
		} else {
			return node_box_type(T,
					     get_box(tl->pts.cut(0, tl->size)));
		}
	} else {
		;
	}
}
} // namespace psi

#endif // PSI_BASE_TREE_IMPL_LEAF_OP_HPP_
