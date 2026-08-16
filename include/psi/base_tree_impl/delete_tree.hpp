#ifndef PSI_BASE_TREE_IMPL_DELETE_TREE_HPP_
#define PSI_BASE_TREE_IMPL_DELETE_TREE_HPP_

#include "../base_tree.h"
#include "dependence/concepts.h"

namespace psi
{

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::delete_tree_wrapper()
{
	if (this->root_ == nullptr) {
		return;
	}
	delete_tree_nodes<leaf_type, interior_type>();
	this->tree_box_ = get_empty_box();
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::delete_tree_nodes()
{
	if (this->root_ == nullptr) {
		return;
	}
	delete_tree_recursive<leaf_type, interior_type>(this->root_);
	this->root_ = nullptr;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio> //* delete tree in parallel
template <typename leaf_type, is_binary_node interior_type, bool granularity>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::delete_tree_recursive(
	node *T)
{
	if (T == nullptr) {
		return;
	}
	if (T->is_leaf) {
		free_node<leaf_type>(T);
	} else {
		interior_type *ti = static_cast<interior_type *>(T);
		// NOTE: enable granularity control by default, if it is
		// disabled, always delete in parallel
		parlay::par_do_if(
			force_parallel_recursion<interior_type, granularity>(
				ti),
			[&] {
				delete_tree_recursive<leaf_type, interior_type>(
					ti->left);
			},
			[&] {
				delete_tree_recursive<leaf_type, interior_type>(
					ti->right);
			});
		free_node<interior_type>(T);
	}
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio> //* delete tree in parallel
template <typename leaf_type, is_multi_node interior_type, bool granularity>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::delete_tree_recursive(
	node *T)
{
	if (T == nullptr) {
		return;
	}
	if (T->is_leaf) {
		free_node<leaf_type>(T);
	} else {
		interior_type *ti = static_cast<interior_type *>(T);

		// NOTE: enable granularity control by default, if it is
		// disabled, always delete in parallel
		parlay::parallel_for(
			0, ti->tree_nodes.size(),
			[&](size_t i) {
				delete_tree_recursive<leaf_type, interior_type>(
					ti->tree_nodes[i]);
			},
			force_parallel_recursion<interior_type, granularity>(ti)
				? 1
				: interior_type::get_regions());

		free_node<interior_type>(T);
	}
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio> //* delete tree in parallel
template <typename leaf_type, is_dynamic_node interior_type, bool granularity>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::delete_tree_recursive(
	node *T)
{
	if (T == nullptr) {
		return;
	}
	if (T->is_leaf) {
		free_node<leaf_type>(T);
	} else {
		interior_type *ti = static_cast<interior_type *>(T);

		// NOTE: enable granularity control by default, if it is
		// disabled, always delete in parallel
		parlay::parallel_for(
			0, ti->tree_nodes.size(),
			[&](size_t i) {
				delete_tree_recursive<leaf_type, interior_type>(
					ti->tree_nodes[i]);
			},
			force_parallel_recursion<interior_type, granularity>(ti)
				? 1
				: ti->tree_nodes.size());

		free_node<interior_type>(T);
	}
}
} // namespace psi

#endif // PSI_BASE_TREE_IMPL_DELETE_TREE_HPP_
