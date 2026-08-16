#ifndef PSI_BASE_TREE_IMPL_REBUILD_HPP_
#define PSI_BASE_TREE_IMPL_REBUILD_HPP_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#include "psi/base_tree.h"

namespace psi
{

/* rebuild the tree */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type, bool granularity>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::prepare_rebuild(
	node *T, points_type &wx, points_type &wo)
{
	wo = points_type::uninitialized(T->size);
	wx = points_type::uninitialized(T->size);
	flatten_rec<leaf_type, interior_type, slice_type, granularity>(
		T, parlay::make_slice(wx));
	delete_tree_recursive<leaf_type, interior_type, granularity>(T);
	return;
}

/* rebuild with new input in */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type, bool granularity>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::prepare_rebuild(
	node *T, slice_type in, points_type &wx, points_type &wo)
{
	wo = points_type::uninitialized(T->size + in.size());
	wx = points_type::uninitialized(T->size + in.size());
	parlay::parallel_for(0, in.size(), [&](size_t j) { wx[j] = in[j]; });
	flatten_rec<leaf_type, interior_type, slice_type, granularity>(
		T, wx.cut(in.size(), wx.size()));
	delete_tree_recursive<leaf_type, interior_type, granularity>(T);
	return;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type, typename PrepareFunc,
	  typename... Args>
node *base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::rebuild_with_insert(
	node *T, PrepareFunc prepare_func, slice_type in, Args &&...args)
{
	points_type w_in, w_out;
	prepare_rebuild<leaf_type, interior_type>(T, in, w_in, w_out);
	auto additional_arg = prepare_func(T, w_in, w_out);
	static_assert(
		std::is_invocable_v<decltype(&DerivedTree::build_recursive),
				    DerivedTree *, slice_type, slice_type,
				    Args &&..., box_type>);
	return static_cast<DerivedTree *>(this)->build_recursive(
		parlay::make_slice(w_in), parlay::make_slice(w_out),
		std::forward<Args>(args)..., additional_arg);
}

/*
 * PARA: when granularity set to false, it will disable the default value for
 * granularity size, which is base_type::serial_build_cutoff; instead, it will
 * check whether the force parallel flag has been enabled
 */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type, bool granularity,
	  typename... Args>
node *base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::rebuild_single_tree(
	node *T, Args &&...args)
{
	points_type wx, wo;
	prepare_rebuild<leaf_type, interior_type, granularity>(T, wx, wo);
	static_assert(
		std::is_invocable_v<decltype(&DerivedTree::build_recursive),
				    DerivedTree *, slice_type, slice_type,
				    Args &&...>);
	return static_cast<DerivedTree *>(this)->build_recursive(
		parlay::make_slice(wx), parlay::make_slice(wo),
		std::forward<Args>(args)...);
}

/*
 * traverse a tree, if it satisfy the condition, then rebuild a binary
 * tree
 * PARA: if allow_enable_rebuild enabled, this method will re-balance the
 * tree; otherwise, it flattens all sparse node into leaf nodes
 */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_binary_node interior_type, bool granularity,
	  typename PrepareFunc, typename... Args>
node *
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::rebuild_tree_recursive(
	node *T, PrepareFunc &&prepare_func, bool const allow_inba_rebuild,
	Args &&...args)
{
	if (T->is_leaf) {
		return T;
	}

	interior_type *ti = static_cast<interior_type *>(T);
	/* rebuild the tree if it is sparse or imbalance */
	if (sparse_node(0, ti->size) ||
	    (allow_inba_rebuild && imbalance_node(ti->left->size, ti->size))) {
		return rebuild_single_tree<leaf_type, interior_type,
					   granularity>(
			T, std::forward<Args>(args)...);
	}

	auto const [left_args, right_args] =
		prepare_func(T, std::forward<Args>(args)...);

	node *L, *R;
	parlay::par_do_if(
		force_parallel_recursion<interior_type, granularity>(ti),
		[&] {
			L = std::apply(
				[&](auto &&...left_args) {
					return rebuild_tree_recursive<
						leaf_type, interior_type,
						granularity>(
						ti->left, prepare_func,
						allow_inba_rebuild,
						std::forward<
							decltype(left_args)>(
							left_args)...);
				},
				left_args);
		},
		[&] {
			R = std::apply(
				[&](auto &&...right_args) {
					return rebuild_tree_recursive<
						leaf_type, interior_type,
						granularity>(
						ti->right, prepare_func,
						allow_inba_rebuild,
						std::forward<
							decltype(right_args)>(
							right_args)...);
				},
				right_args);
		});

	update_interior<interior_type>(T, L, R);
	return T;
}

/* rebuild a multi-node tree */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_multi_node interior_type, bool granularity,
	  typename PrepareFunc, typename... Args>
node *
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::rebuild_tree_recursive(
	node *T, PrepareFunc &&prepare_func,
	[[maybe_unused]] bool const allow_inba_rebuild, Args &&...args)
{
	if (T->is_leaf) {
		return T;
	}

	interior_type *ti = static_cast<interior_type *>(T);
	if (sparse_node(0, ti->size)) {
		return rebuild_single_tree<leaf_type, interior_type,
					   granularity>(
			T, std::forward<Args>(args)...);
	}

	typename interior_type::node_arr_type new_nodes;
	parlay::parallel_for(
		0, interior_type::get_regions(),
		[&](bucket_type i) {
			auto const new_args =
				prepare_func(T, i, std::forward<Args>(args)...);
			std::apply(
				[&](auto &&...new_args) {
					new_nodes[i] = rebuild_tree_recursive<
						leaf_type, interior_type,
						granularity>(
						ti->tree_nodes[i], prepare_func,
						allow_inba_rebuild,
						std::forward<
							decltype(new_args)>(
							new_args)...);
				},
				new_args);
		},
		force_parallel_recursion<interior_type, granularity>(ti)
			? 1
			: interior_type::get_regions());

	update_interior<interior_type>(T, new_nodes);
	return T;
}

} /* namespace psi */

#endif /* PSI_BASE_TREE_IMPL_REBUILD_HPP_ */
