#ifndef PSI_BASE_TREE_IMPL_VALIDATION_HPP_
#define PSI_BASE_TREE_IMPL_VALIDATION_HPP_

#include <parlay/parallel.h>

#include <concepts>

#include "psi/base_tree.h"
#include "parlay/primitives.h"

namespace psi
{

template <typename Traits, typename DerivedTree>
template <typename leaf_type, typename interior_type>
typename base_tree<Traits, DerivedTree>::box_type
base_tree<Traits, DerivedTree>::check_box(node *T, box_type const &box)
{
	if (T->is_leaf) {
		auto const *tl = static_cast<leaf_type const *>(T);
		auto const node_box =
			Traits::template get_box<leaf_type, interior_type>(T);
		if constexpr (has_box<typename leaf_type::at_type>) {
			assert(same_box(node_box, tl->get_box()));
		} else {
			assert(within_box(node_box, box));
		}
		return node_box;
	}
	interior_type *ti = static_cast<interior_type *>(T);
	assert(!ti->get_parallel_flag_ini_status()); /* ensure that */
						     /*
						      * uninitialized force
						      * parallelism
						      */
	if constexpr (is_binary_node<interior_type> &&
		      !has_box<typename interior_type::at_type>) { /* use */
		/* hyperplane */
		/* as splitter */
		box_type lbox(box), rbox(box);
		lbox.second.pnt[ti->split.second] = ti->split.first;
		rbox.first.pnt[ti->split.second] = ti->split.first;
		box_type const left_return_box =
			check_box<leaf_type, interior_type>(ti->left, lbox);
		box_type const right_return_box =
			check_box<leaf_type, interior_type>(ti->right, rbox);
		box_type const new_box =
			get_box(left_return_box, right_return_box);

		assert(within_box(left_return_box, lbox));
		assert(within_box(right_return_box, rbox));
		assert(within_box(new_box, box));
		return new_box;
	} else if constexpr (is_binary_node<interior_type> &&
			     has_box<typename interior_type::at_type>) { /* kd
									  */
		/*
		 * with
		 * box
		 */
		auto left_box =
			retrieve_box<leaf_type, interior_type>(ti->left);
		auto right_box =
			retrieve_box<leaf_type, interior_type>(ti->right);
		box_type const left_return_box =
			check_box<leaf_type, interior_type>(ti->left, left_box);
		box_type const right_return_box =
			check_box<leaf_type, interior_type>(ti->right,
							    right_box);
		box_type const new_box =
			get_box(left_return_box, right_return_box);
		assert(same_box(
			left_return_box,
			retrieve_box<leaf_type, interior_type>(ti->left)));
		assert(same_box(
			right_return_box,
			retrieve_box<leaf_type, interior_type>(ti->right)));
		assert(same_box(new_box, ti->get_box()));
		return new_box;
	} else if constexpr (is_multi_node<interior_type> &&
			     !has_box<typename interior_type::
					      at_type>) { /* orth
							     without
							   */
		/* box */
		box_seq_type new_box(
			ti->template compute_subregions<box_seq_type>(box));
		box_seq_type return_box_seq(new_box.size());
		assert(new_box.size() == ti->tree_nodes.size());
		for (size_t i = 0; i < ti->tree_nodes.size(); i++) {
			return_box_seq[i] = check_box<leaf_type, interior_type>(
				ti->tree_nodes[i], new_box[i]);
			assert(within_box(return_box_seq[i], new_box[i]));
		}
		auto return_box = get_box(return_box_seq);
		assert(within_box(return_box, box));
		return return_box;
	} else if constexpr (is_multi_node<interior_type> &&
			     has_box<typename interior_type::at_type>) { /* orth
									  */
		/*
		 * with
		 * box
		 */
		box_seq_type new_box(
			ti->template compute_subregions<box_seq_type>(box));
		box_seq_type return_box_seq(new_box.size());
		assert(new_box.size() == ti->tree_nodes.size());
		for (size_t i = 0; i < ti->tree_nodes.size(); i++) {
			return_box_seq[i] = check_box<leaf_type, interior_type>(
				ti->tree_nodes[i], new_box[i]);
			assert(same_box(return_box_seq[i],
					retrieve_box<leaf_type, interior_type>(
						ti->tree_nodes[i])));
		}
		auto return_box = get_box(return_box_seq);
		assert(same_box(return_box, ti->get_box()));
		return return_box;
	} else {
		static_assert(is_binary_node<interior_type> ||
				      is_multi_node<interior_type>,
			      "check_box supports only binary and multi-way "
			      "interior nodes");
		return get_empty_box();
	}
}

template <typename Traits, typename DerivedTree>
template <typename leaf_type, typename interior_type>
size_t base_tree<Traits, DerivedTree>::check_size(node *T)
{
	if (T->is_leaf) {
		return T->size;
	}
	if constexpr (is_binary_node<interior_type>) {
		interior_type *ti = static_cast<interior_type *>(T);
		assert(!ti->get_parallel_flag_ini_status());
		size_t l, r;
		parlay::par_do(
			[&l, &ti] {
				l = check_size<leaf_type, interior_type>(
					ti->left);
			},
			[&r, &ti] {
				r = check_size<leaf_type, interior_type>(
					ti->right);
			});
		assert(l + r == T->size);
		return T->size;
	} else {
		interior_type *ti = static_cast<interior_type *>(T);
		assert(!ti->get_parallel_flag_ini_status());
		size_t sum = 0;
		for (bucket_type i = 0; i < ti->tree_nodes.size(); ++i) {
			sum += check_size<leaf_type, interior_type>(
				ti->tree_nodes[i]);
		}
		assert(sum == T->size);
		return T->size;
	}
}

template <typename Traits, typename DerivedTree>
template <typename leaf_type, typename interior_type>
void base_tree<Traits, DerivedTree>::check_tree_same_sequential(node *T,
								int dim)
{
	if (T->is_leaf) {
		return;
	}
	if constexpr (is_binary_node<interior_type>) {
		interior_type *ti = static_cast<interior_type *>(T);
		if (ti->split.second != dim) {
			std::cout << int(ti->split.second) << " " << int(dim)
				  << " " << ti->size;
		}
		assert(ti->split.second == dim);
		/* TODO: maybe need to add the split rule in the base tree? */
		dim = (dim + 1) % num_dims;
		parlay::par_do_if(
			T->size > 1000,
			[&]() {
				check_tree_same_sequential<leaf_type,
							   interior_type>(
					ti->left, dim);
			},
			[&]() {
				check_tree_same_sequential<leaf_type,
							   interior_type>(
					ti->right, dim);
			});
	} else {
		assert(is_multi_node<interior_type>);
		interior_type *ti = static_cast<interior_type *>(T);
		assert(std::cmp_equal((1 << ti->split.size()),
				      ti->tree_nodes.size()));
		for (size_t i = 0; i < ti->split.size(); i++) {
			assert(ti->split[i].second == dim);
			dim += 1;
		}
		assert(dim == num_dims);
		for (size_t i = 0; i < ti->tree_nodes.size(); i++) {
			check_tree_same_sequential<leaf_type, interior_type>(
				ti->tree_nodes[i], 0);
		}
	}
	return;
}

template <typename Traits, typename DerivedTree>
template <typename leaf_type, typename interior_type, typename SplitRule>
void base_tree<Traits, DerivedTree>::validate()
{
	std::cout << ">>> begin validate tree\n" << std::flush;

	check_size<leaf_type, interior_type>(this->root_);
	std::cout << "Correct size\n" << std::flush;

	/* tree property */
	if constexpr (is_binary_node<interior_type> ||
		      is_multi_node<interior_type>) {
		if (legal_box(check_box<leaf_type, interior_type>(
			    this->root_, this->tree_box_))) {
			std::cout << "Correct bounding Box\n" << std::flush;
		} else {
			std::cout << "wrong bounding Box\n" << std::flush;
			abort();
		}

		/*
		 * used to check rotate dimension
		 * For kdtree binary node, the dummy node may break the rotation
		 * manner, since if current dimension is un-splitable, one has
		 * to switch to another dimension
		 */
		if constexpr (is_rotate_dim_split<
				      typename SplitRule::dim_rule_type> &&
			      is_multi_node<interior_type>) {
			check_tree_same_sequential<leaf_type, interior_type>(
				this->root_, 0);
			std::cout << "Correct rotate dimension\n" << std::flush;
		}
	} else {
		/*
		 * p_tree lands here: its node type is CPAM's, which none of
		 * the branches above understand, and CPAM is not ours to walk.
		 * Say so rather than aborting on what reads like an
		 * unreachable case -- p_tree is covered by the black-box
		 * property tests under tests/unit instead.
		 */
		std::cout << "no structural validation for this node type; "
			     "see tests/unit\n"
			  << std::flush;
	}

	std::cout << "<<< end validate tree\n" << std::flush;
	return;
}

template <typename Traits, typename DerivedTree>
template <typename leaf_type, typename interior_type>
size_t base_tree<Traits, DerivedTree>::get_tree_height()
{
	if (this->root_ == nullptr) {
		return 0;
	}
	size_t deep = 0;
	return get_max_tree_depth<leaf_type, interior_type>(this->root_, deep);
}

template <typename Traits, typename DerivedTree>
template <typename leaf_type, typename interior_type>
size_t base_tree<Traits, DerivedTree>::get_max_tree_depth(node *T, size_t deep)
{
	if (T == nullptr) {
		return deep;
	}
	if (T->is_leaf) {
		return deep;
	}

	interior_type *ti = static_cast<interior_type *>(T);
	if constexpr (is_binary_node<interior_type>) {
		int l = get_max_tree_depth<leaf_type, interior_type>(ti->left,
								     deep + 1);
		int r = get_max_tree_depth<leaf_type, interior_type>(ti->right,
								     deep + 1);
		return std::max(l, r);
	} else if constexpr (is_multi_node<interior_type>) {
		size_t max_depth = 0;
		for (size_t i = 0; i < ti->tree_nodes.size(); i++) {
			max_depth = std::max(
				max_depth,
				get_max_tree_depth<leaf_type, interior_type>(
					ti->tree_nodes[i], deep + 1));
		}
		return max_depth;
	} else {
		return 0;
	}
}

template <typename Traits, typename DerivedTree>
template <typename leaf_type, typename interior_type>
double base_tree<Traits, DerivedTree>::get_ave_tree_height()
{
	if (this->root_ == nullptr) {
		return 0.0;
	}
	size_t leaf_num = 0, depth_sum = 0;
	count_tree_heights<leaf_type, interior_type>(this->root_, 0, leaf_num,
						     depth_sum);
	return leaf_num == 0 ? 0.0 : double(depth_sum) / double(leaf_num);
}

template <typename Traits, typename DerivedTree>
template <typename leaf_type, typename interior_type>
size_t base_tree<Traits, DerivedTree>::count_tree_nodes_num(node *T)
{
	if (T == nullptr) {
		return 0;
	}
	if (T->is_leaf) {
		return 1;
	}

	interior_type *ti = static_cast<interior_type *>(T);
	if constexpr (is_binary_node<interior_type>) {
		size_t l, r;
		parlay::par_do(
			[&]() {
				l = count_tree_nodes_num<leaf_type,
							 interior_type>(
					ti->left);
			},
			[&]() {
				r = count_tree_nodes_num<leaf_type,
							 interior_type>(
					ti->right);
			});
		return l + r + 1;
	} else {
		size_t sum = 0;
		for (int i = 0; i < ti->tree_nodes.size(); i++) {
			sum += count_tree_nodes_num<leaf_type, interior_type>(
				ti->tree_nodes[i]);
		}
		return sum + 1;
	}
}

template <typename Traits, typename DerivedTree>
template <typename leaf_type, typename interior_type>
void base_tree<Traits, DerivedTree>::count_tree_heights(node *T, size_t deep,
							size_t &leaf_num,
							size_t &depth_sum)
{
	if (T->is_leaf) {
		leaf_num++;
		depth_sum += deep;
		return;
	}

	interior_type *ti = static_cast<interior_type *>(T);
	if constexpr (is_binary_node<interior_type>) {
		count_tree_heights<leaf_type, interior_type>(
			ti->left, deep + 1, leaf_num, depth_sum);
		count_tree_heights<leaf_type, interior_type>(
			ti->right, deep + 1, leaf_num, depth_sum);
	} else if constexpr (is_multi_node<interior_type>) {
		for (size_t i = 0; i < ti->tree_nodes.size(); i++) {
			count_tree_heights<leaf_type, interior_type>(
				ti->tree_nodes[i], deep + 1, leaf_num,
				depth_sum);
		}
	}
}

} /* namespace psi */

#endif /* PSI_BASE_TREE_IMPL_VALIDATION_HPP_ */
