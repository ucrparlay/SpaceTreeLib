#ifndef PSI_BASE_TREE_IMPL_TREE_OP_BUILD_INNER_TREE_HPP_
#define PSI_BASE_TREE_IMPL_TREE_OP_BUILD_INNER_TREE_HPP_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#include "../../base_tree.h"

namespace psi
{
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_binary_node interior_type, typename ReturnType>
ReturnType base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::build_inner_tree(
	bucket_type idx, hyper_plane_seq_type const &pivots,
	parlay::sequence<ReturnType> const &tree_nodes)
{
	if (idx > pivot_num) {
		assert(idx - pivot_num - 1 < bucket_num);
		return tree_nodes[idx - pivot_num - 1];
	}

	ReturnType const L = build_inner_tree<leaf_type, interior_type>(
		idx << 1, pivots, tree_nodes);
	ReturnType const R = build_inner_tree<leaf_type, interior_type>(
		idx << 1 | 1, pivots, tree_nodes);

	return alloc_interior_node<interior_type>(L, R, pivots[idx]);
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <is_multi_node interior_type>
node *base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::build_inner_tree(
	bucket_type idx, hyper_plane_seq_type const &pivots,
	parlay::sequence<node *> const &tree_nodes)
{
	assert(idx < pivot_num + bucket_num + 1);

	if (idx > pivot_num) {
		return tree_nodes[idx - pivot_num - 1];
	}

	typename DerivedTree::orth_node_arr_type multi_nodes;
	typename DerivedTree::splitter_type split;
	for (dims_type i = 0; i < DerivedTree::node_regions; ++i) {
		multi_nodes[i] = build_inner_tree<interior_type>(
			idx * DerivedTree::node_regions + i, pivots,
			tree_nodes);
	}
	for (dims_type i = 0; i < DerivedTree::splitter_num; ++i) {
		split[i] = pivots[idx * (1 << i)];
		assert(i == 0 ||
		       pivots[idx * (1 << i)] == pivots[idx * (1 << i) + 1]);
	}

	return alloc_interior_node<interior_type>(multi_nodes, split);
}
} // namespace psi

#endif // PSI_BASE_TREE_IMPL_TREE_OP_BUILD_INNER_TREE_HPP_
