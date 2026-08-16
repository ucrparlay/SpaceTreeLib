#ifndef PSI_BASE_TREE_IMPL_TREE_OP_FLATTEN_HPP_
#define PSI_BASE_TREE_IMPL_TREE_OP_FLATTEN_HPP_

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

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <supports_force_parallel interior_type, bool granularity>
inline bool
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::force_parallel_recursion(
	interior_type const *ti)
{
#ifndef DISABLE_BATCH_DELETE_SIZE_OPT
	return (granularity && ti->size > serial_build_cutoff) ||
	       (!granularity && ti->force_parallel());
#else
	return (granularity) || (!granularity && ti->force_parallel());
#endif // !DISABLE_BATCH_DELETE_SIZE_OPT
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_binary_node interior_type, typename Range,
	  bool granularity>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::flatten_rec(node *T,
								     Range out)
{
	assert(T->size == out.size());

	if (T->size == 0)
		return;

	if (T->is_leaf) {
		extract_points_in_leaf<leaf_type>(T, out);
		return;
	}

	interior_type *ti = static_cast<interior_type *>(T);
	assert(ti->size == ti->left->size + ti->right->size);
	parlay::par_do_if(
		// WARN: check parallelisim using node size can be biased
		force_parallel_recursion<interior_type, granularity>(ti),
		[&]() {
			flatten_rec<leaf_type, interior_type>(
				ti->left, out.cut(0, ti->left->size));
		},
		[&]() {
			flatten_rec<leaf_type, interior_type>(
				ti->right, out.cut(ti->left->size, ti->size));
		});

	return;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type, typename Range,
	  bool granularity>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::flatten_rec(node *T,
								     Range out)
	requires(!is_binary_node<interior_type>)
{
	assert(T->size == out.size());

	if (T->size == 0)
		return;

	if (T->is_leaf) {
		extract_points_in_leaf<leaf_type>(T, out);
		return;
	}

	interior_type *ti = static_cast<interior_type *>(T);

	assert(ti->size == std::accumulate(ti->tree_nodes.begin(),
					   ti->tree_nodes.end(),
					   static_cast<size_t>(0),
					   [](size_t acc, node *n) -> size_t {
						   return acc + n->size;
					   }));

	parlay::parallel_for(
		0, ti->tree_nodes.size(),
		[&](bucket_type i) {
			size_t start = 0;
			for (bucket_type j = 0; j < i; ++j) {
				start += ti->tree_nodes[j]->size;
			}
			flatten_rec<leaf_type, interior_type, Range>(
				ti->tree_nodes[i],
				out.cut(start,
					start + ti->tree_nodes[i]->size));
		},
		force_parallel_recursion<interior_type, granularity>(ti)
			? 1
			: ti->get_sub_tree_num());

	return;
}

// NOTE: for multi node @T, it only flatten the subtree with id @idx to @out
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_multi_node interior_type, typename Range,
	  bool granularity>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::partial_flatten(
	node *T, Range out, bucket_type idx)
{
	if (idx == 1) {
		assert(T->size == out.size());
		flatten_rec<leaf_type, interior_type>(T, out.cut(0, T->size));
		return;
	} else if (idx >= interior_type::get_regions()) {
		node *ns = static_cast<interior_type *>(T)
				   ->tree_nodes[idx -
						interior_type::get_regions()];
		assert(ns->size == out.size());
		flatten_rec<leaf_type, interior_type>(ns, out.cut(0, ns->size));
		return;
	}

	interior_type *ti = static_cast<interior_type *>(T);
	size_t l_size = ti->merge_size(idx << 1),
	       r_size = ti->merge_size(idx << 1 | 1);
	assert(l_size + r_size == out.size());
	parlay::par_do_if(
		force_parallel_recursion<interior_type, granularity>(
			static_cast<interior_type *>(T)),
		[&]() {
			partial_flatten<leaf_type, interior_type>(
				T, out.cut(0, l_size), idx << 1);
		},
		[&]() {
			partial_flatten<leaf_type, interior_type>(
				T, out.cut(l_size, l_size + r_size),
				idx << 1 | 1);
		});
	return;
}
} // namespace psi

#endif // PSI_BASE_TREE_IMPL_TREE_OP_FLATTEN_HPP_
