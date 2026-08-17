#ifndef PSI_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_
#define PSI_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_

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
/* update the info of T by new children L and R */
template <typename Traits, typename DerivedTree>
template <is_binary_node interior_type, bool UpdateParFlag>
inline void base_tree<Traits, DerivedTree>::update_interior(node *T, node *L,
							    node *R)
{
	assert(!T->is_leaf);
	interior_type *ti = static_cast<interior_type *>(T);
	if constexpr (UpdateParFlag) {
		ti->reset_parallel_flag();
	}
	ti->size = L->size + R->size;
	ti->left = L;
	ti->right = R;
	ti->update_aug(L, R);
	return;
}

template <typename Traits, typename DerivedTree>
template <is_binary_node interior_type, bool UpdateParFlag>
inline void
base_tree<Traits, DerivedTree>::update_interior(node *T, node_box_type const &L,
						node_box_type const &R)
{
	assert(!T->is_leaf);
	interior_type *ti = static_cast<interior_type *>(T);
	if constexpr (UpdateParFlag) {
		ti->reset_parallel_flag();
	}
	ti->size = L.first->size + R.first->size;
	ti->left = L.first;
	ti->right = R.first;
	ti->update_aug(L.first, R.first);
	return;
}

template <typename Traits, typename DerivedTree>
template <is_multi_node interior_type, bool UpdateParFlag>
inline void base_tree<Traits, DerivedTree>::update_interior(
	node *T, typename interior_type::node_arr_type const &new_nodes)
{
	assert(!T->is_leaf);
	interior_type *ti = static_cast<interior_type *>(T);
	if constexpr (UpdateParFlag) {
		ti->reset_parallel_flag();
	}
	/* size_t init: a plain 0 makes the accumulator int and truncates. */
	ti->size = std::accumulate(
		new_nodes.begin(), new_nodes.end(), size_t{0},
		[](size_t acc, node *n) -> size_t { return acc + n->size; });
	ti->tree_nodes = new_nodes;
	ti->update_aug(new_nodes);
	return;
}
} /* namespace psi */

#endif /* PSI_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_ */
