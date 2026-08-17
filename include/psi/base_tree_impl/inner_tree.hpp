#ifndef PSI_BASE_TREE_IMPL_INNER_TREE_HPP_
#define PSI_BASE_TREE_IMPL_INNER_TREE_HPP_

#include <utility>

#include "psi/base_tree.h"
#include "psi/dependence/concepts.h"

namespace psi
{
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename interior_type>
struct base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::inner_tree {
	using base_type = base_tree<Point, DerivedTree, SkHeight, ImbaRatio>;
	inner_tree()
	    : tags_num(0), tags(node_tag_seq_type::uninitialized(
				   pivot_num + bucket_num + 1)),
	      sums_tree(
		      parlay::sequence<balls_type>(pivot_num + bucket_num + 1)),
	      rev_tag(bucket_seq_type::uninitialized(bucket_num))
	{
	}

	/* helpers */
	bool assert_size(node *T) const
	{
		if (T->is_leaf) {
			leaf_type *ti = static_cast<leaf_type *>(T);
			assert(T->size <= ti->pts.size() &&
			       T->size <= leaf_capacity);
			return true;
		}
		interior_type *ti = static_cast<interior_type *>(T);
		if constexpr (is_binary_node<interior_type>) {
			assert(ti->size == ti->left->size + ti->right->size);
		} else if constexpr (is_multi_node<interior_type>) {
			assert(std::cmp_equal(
				ti->size,
				std::accumulate(ti->tree_nodes.begin(),
						ti->tree_nodes.end(), 0,
						[](size_t sum, node *T) {
							return sum + T->size;
						})));
		} else {
			static_assert(is_binary_node<interior_type> ||
				      is_multi_node<interior_type>);
		}
		return true;
	}

	/* cores */
	inline void reset_tags_num()
	{
		tags_num = 0;
	}

	box_type get_box_by_region_idx(int const idx, box_type const &box)
	{
		if constexpr (is_binary_node<interior_type>) {
			assert(get_depth_by_index(10) == 3);
			box_type bx(box);
			int h = get_depth_by_index(idx);
			for (int i = h - 1, new_idx = 1; i >= 0; i--) {
				int local_id = (idx >> i) & 1;
				auto ti = static_cast<interior_type *>(
					tags[new_idx].first);
				auto &target = local_id ? bx.first : bx.second;
				target.pnt[ti->split.second] = ti->split.first;
				new_idx = new_idx << 1 | local_id;
				assert(new_idx <= idx);
			}
			return std::move(bx);
		} else if constexpr (is_multi_node<interior_type>) {
			assert(get_depth_by_index(19) == 4);
			assert(num_dims ==
			       static_cast<bucket_type>(std::log2(
				       interior_type::get_regions())));

			box_type bx(box);
			bucket_type h = get_depth_by_index(idx);
			for (bucket_type i = h, new_idx = 1; i > 0;
			     i -= num_dims) {
				bucket_type local_id = (idx >> (i - num_dims)) &
						       ((1 << num_dims) - 1);
				auto ti = static_cast<interior_type *>(
					tags[new_idx].first);
				ti->modify_box_by_id(local_id, bx);
				new_idx = new_idx << num_dims | local_id;
			}
			return std::move(bx);
		} else {
			static_assert(is_binary_node<interior_type> ||
				      is_multi_node<interior_type>);
		}
	}

	/*
	 * Each node in the skeleton receives a tag
	 * A leaf_type node receives the tag < BUCKETNUM
	 * All internal node has tag == BUCKETNUM
	 */
	void assign_node_tag(node *T, bucket_type idx)
	{
		if (T->is_leaf || idx > pivot_num) {
			assert(tags_num < bucket_num);
			tags[idx] = node_tag_type(T, tags_num);
			rev_tag[tags_num++] = idx; /* cannot remove */
			return;
		}

		/* INFO: BUCKET ID in [0, bucket_num) */
		tags[idx] = node_tag_type(T, bucket_num);
		interior_type *ti = static_cast<interior_type *>(T);
		if constexpr (is_binary_node<interior_type>) {
			assign_node_tag(ti->left, idx << 1);
			assign_node_tag(ti->right, idx << 1 | 1);
		} else if constexpr (is_multi_node<interior_type>) {
			for (bucket_type i = 0;
			     i < interior_type::get_regions(); ++i) {
				assign_node_tag(
					ti->tree_nodes[i],
					idx * interior_type::get_regions() + i);
			}
		} else {
			static_assert(is_binary_node<interior_type> ||
				      is_multi_node<interior_type>);
		}
		return;
	}

	/*
	 * reduce sums is travsersal of the skeleton with counting the
	 * points sieved onto every node, it is good to determine whether we
	 * need forcing parallel in the following operations, e.g.,
	 * flatten/rebuild the tree.
	 */
	void reduce_sums(bucket_type idx)
	{
		if (idx > pivot_num || tags[idx].first->is_leaf) {
			assert(tags[idx].second < bucket_num);
			sums_tree[idx] = sums[tags[idx].second];

			/*
			 * no need to update the parallel flag here, as it
			 * is either a leaf node or it will be handled by
			 * recursive calls
			 */
			return;
		}

		if constexpr (is_binary_node<interior_type>) {
			reduce_sums(idx << 1);
			reduce_sums(idx << 1 | 1);
			sums_tree[idx] =
				sums_tree[idx << 1] + sums_tree[idx << 1 | 1];
		} else if constexpr (is_multi_node<interior_type>) {
			sums_tree[idx] = 0;
			for (bucket_type i = 0;
			     i < interior_type::get_regions(); ++i) {
				reduce_sums(idx * interior_type::get_regions() +
					    i);
				sums_tree[idx] += sums_tree
					[idx * interior_type::get_regions() +
					 i];
			}
		}

		/*
		 * Don't add force parallel here as it not precise:
		 * whether a node should be rebuilt depends on only the tomb
		 * status, rather than the points sieved to that node
		 */
		return;
	}

	static inline bucket_type get_depth_by_index(bucket_type idx)
	{
		bucket_type h = 0;
		while (idx > 1) {
			idx >>= 1;
			++h;
		}
		return h;
	}

	/*
	 * a bucket/leaf has id bucket_num+1
	 * a node needs to be rebuilt has id bucket_num+2
	 * otherwise, it has id bucket_num
	 */
	template <typename ViolateFunc>
	void pick_tag(bucket_type idx, ViolateFunc &&violate_func)
	{
		if (idx > pivot_num || tags[idx].first->is_leaf) {
			tags[idx].second = bucket_num + 1;
			rev_tag[tags_num++] = idx;
			return;
		}

		assert(tags[idx].second == bucket_num &&
		       (!tags[idx].first->is_leaf));
		if (invoke_with_optional_arg<bool>(violate_func, idx)) {
			tags[idx].second = bucket_num + 2;
			rev_tag[tags_num++] = idx;
			return;
		}

		if constexpr (is_binary_node<interior_type>) {
			pick_tag(idx << 1, violate_func);
			pick_tag(idx << 1 | 1, violate_func);
		} else if constexpr (is_multi_node<interior_type>) {
			for (bucket_type i = 0;
			     i < interior_type::get_regions(); ++i) {
				pick_tag(idx * interior_type::get_regions() + i,
					 violate_func);
			}
		}
		return;
	}

	template <typename... Args>
	void tag_imbalance_node(Args &&...args)
	{
		reduce_sums(1);
		reset_tags_num();
		pick_tag(1, std::forward<Args>(args)...);
		assert(assert_size(tags[1].first));
		return;
	}

	/*
	 * the node which needs to be rebuilt has tag bucket_num+3
	 * the *bucket* node whose ancestor has been rebuilt has tag
	 * bucket_num+2 the *bucket* node whose ancestor has not been ... has
	 * bucket_num+1 otherwise, it's bucket_num
	 */
	template <bool SetParallelFlag, typename ViolateFunc>
	void mark_tomb(bucket_type idx, bucket_type &re_num,
		       size_t &tot_re_size, box_seq_type &box_seq,
		       box_type const &box, bool has_tomb,
		       ViolateFunc &&violate_func)
	{
		if (idx > pivot_num || tags[idx].first->is_leaf) {
			assert(tags[idx].second >= 0 &&
			       tags[idx].second < bucket_num);
			if (!has_tomb) {
				tags[idx].second = bucket_num + 2;
				if constexpr (SetParallelFlag) { /* INFO: the */
								 /*
								  * sub-tree
								  * will be
								  * rebuilt in
								  * the future
								  * and need to
								  * force the
								  * parallisim
								  */
					if (!tags[idx].first->is_leaf) {
						auto ti = static_cast<
							interior_type *>(
							tags[idx].first);
						ti->set_parallel_flag(
							ti->size >
							base_type::
								serial_build_cutoff);
					}
				}
			} else {
				tags[idx].second = bucket_num + 1;
			}

			box_seq[tags_num] = box;
			rev_tag[tags_num++] = idx;
			return;
		}

		/* no need to mark the internal nodes with tag bucket_num */
		assert(tags[idx].second == bucket_num &&
		       (!tags[idx].first->is_leaf));
		interior_type *ti =
			static_cast<interior_type *>(tags[idx].first);
		if (has_tomb &&
		    invoke_with_optional_arg<bool>(violate_func, idx)) {
			tags[idx].second = bucket_num + 3;
			has_tomb = false;
			ti->set_parallel_flag(ti->size >
					      base_type::serial_build_cutoff);
			re_num++, tot_re_size += ti->size;
		}

		if constexpr (is_binary_node<interior_type>) {
			box_cut_type box_cut(box, ti->split, true);
			mark_tomb<SetParallelFlag>(idx << 1, re_num,
						   tot_re_size, box_seq,
						   box_cut.get_first_box_cut(),
						   has_tomb, violate_func);
			mark_tomb<SetParallelFlag>(idx << 1 | 1, re_num,
						   tot_re_size, box_seq,
						   box_cut.get_second_box_cut(),
						   has_tomb, violate_func);
		} else if constexpr (is_multi_node<interior_type>) {

			for (bucket_type i = 0;
			     i < interior_type::get_regions(); ++i) {
				mark_tomb<SetParallelFlag>(
					idx * interior_type::get_regions() + i,
					re_num, tot_re_size, box_seq,
					ti->get_box_by_region_id(i, box),
					has_tomb, violate_func);
			}
		}
		return;
	}

	/*
	 * the node which needs to be rebuilt has tag bucket_num+3
	 * the *bucket* node whose ancestor has been rebuilt has tag
	 * bucket_num+2 the *bucket* node whose ancestor has not been ... has
	 * bucket_num+1 otherwise, it's bucket_num
	 * TODO: maybe we can make tagImbalance node and
	 * tagImbalancedNodeDeletion together
	 */

	template <bool SetParallelFlag, typename... Args>
	auto tag_imbalance_node_deletion(Args &&...args)
	{
		reduce_sums(1);
		reset_tags_num();
		bucket_type re_num = 0;
		size_t tot_re_size = 0;
		mark_tomb<SetParallelFlag>(1, re_num, tot_re_size,
					   std::forward<Args>(args)...);
		return std::make_pair(re_num, tot_re_size);
	}

	void tag_puffy_nodes_recursive(bucket_type idx)
	{
		if (idx > pivot_num || tags[idx].first->is_leaf) {
			tags[idx].second =
				bucket_num +
				2; /* ensure the following update_interior */
				   /* meets the base case */
			if (!tags[idx].first->is_leaf) {
				interior_type *ti =
					static_cast<interior_type *>(
						tags[idx].first);
				ti->set_parallel_flag(
					ti->size >
					base_type::serial_build_cutoff);
			}
			return;
		}

		assert(!tags[idx].first->is_leaf);
		interior_type *ti =
			static_cast<interior_type *>(tags[idx].first);
		ti->set_parallel_flag(ti->size >
				      base_type::serial_build_cutoff);

		if constexpr (is_binary_node<interior_type>) {
			tag_puffy_nodes_recursive(idx << 1);
			tag_puffy_nodes_recursive(idx << 1 | 1);
		} else if constexpr (is_multi_node<interior_type>) {
			for (bucket_type i = 0;
			     i < interior_type::get_regions(); ++i) {
				tag_puffy_nodes_recursive(
					idx * interior_type::get_regions() + i);
			}
		}
		return;
	}

	void tag_puffy_nodes()
	{
		tag_puffy_nodes_recursive(1);
		return;
	}

	/*
	 * update_pointer: update the pointer only, if it contains box,
	 * return
	 * empty box
	 * update_pointer_box: update pointer and box
	 * tag_rebuild_node: update the pointer and box, meanwhile it assign the
	 * imbalance node with a new tag
	 * post_del_update: update the skeleton after rebuild (e.g., size,
	 * children),
	 * which needs to avoid touch the deleted nodes
	 */
	enum update_type {
		update_pointer,
		update_pointer_box,
		tag_rebuild_node,
		post_del_update
	};

	/* update the skeleton based on the @update_type */
	template <update_type ut, bool UpdateParFlag = true,
		  typename ReturnType, typename... Args>
		requires is_pointer_to_node<ReturnType> ||
			 is_node_box<ReturnType, Point>
	ReturnType
	update_inner_tree(parlay::sequence<ReturnType> const &tree_nodes,
			  Args &&...args)
	{
		bucket_type p = 0;
		if constexpr (ut == update_pointer ||
			      ut == update_pointer_box) { /* update the */
							  /*
							   * inner tree nodes or
							   * box
							   */
			return update_inner_tree_recursive<ut, UpdateParFlag>(
				1, tree_nodes, p, [&]() {});
		} else if constexpr (ut == tag_rebuild_node) { /* tag the */
							       /*
								* node that needs
								* to be rebuild
								*/
			this->reset_tags_num();

			auto func_2_rebuild_node =
				[&](auto const &...params) -> void {
				if constexpr (is_binary_node<interior_type>) {
					/*
					 * needs to
					 * save the
					 * box
					 */
					auto const &new_box = std::get<0>(
						std::forward_as_tuple(
							params...));
					auto const &idx = std::get<1>(
						std::forward_as_tuple(
							params...));
					rev_tag[tags_num] = idx;
					find_var<box_seq_type>(
						std::forward<Args>(
							args)...)[tags_num++] =
						new_box;
				} else if constexpr (is_multi_node<
							     interior_type>) {
					/* the */
					/*
					 * box is fixed
					 * in orth
					 * node, no
					 * need to save
					 */
					auto const &idx = std::get<0>(
						std::forward_as_tuple(
							params...));
					rev_tag[tags_num++] = idx;
				}
			};

			return update_inner_tree_recursive<ut, UpdateParFlag>(
				1, tree_nodes, p, func_2_rebuild_node);
		} else if constexpr (ut ==
				     post_del_update) { /* avoid touch the */
							/*
							 * node that has been
							 * deleted
							 */
			/*
			 * PARA: op == 0 -> toggle whether under a rebuild tree
			 * op == 1 -> query current status
			 */
			bool under_rebuild_tree = false;
			return update_inner_tree_recursive<ut, UpdateParFlag>(
				1, tree_nodes, p, [&](bool op) -> bool {
					return op == 0 ? (under_rebuild_tree =
								  !under_rebuild_tree)
						       : under_rebuild_tree;
				});
		} else {
		}
	}

	/* udpate inner tree for binary nodes */
	template <update_type ut, bool UpdateParFlag, typename ReturnType,
		  typename Func>
		requires is_binary_node<interior_type>
	ReturnType update_inner_tree_recursive(
		bucket_type idx, parlay::sequence<ReturnType> const &tree_nodes,
		bucket_type &p, Func &&func)
	{
		/* needs to ensure this success for both insert and delete */
		if (this->tags[idx].second == bucket_num + 1 ||
		    this->tags[idx].second == bucket_num + 2) {
			return tree_nodes[p++];
		}

		if constexpr (ut == post_del_update) {
			if (this->tags[idx].second == bucket_num + 3) {
				func(0); /* close the under_rebuild_tree flag */
				assert(func(1) == true);
			}
		}

		ReturnType const &left =
			update_inner_tree_recursive<ut, UpdateParFlag>(
				idx << 1, tree_nodes, p, func);
		ReturnType const &right =
			update_inner_tree_recursive<ut, UpdateParFlag>(
				idx << 1 | 1, tree_nodes, p, func);

		if constexpr (ut ==
			      update_pointer) { /* only update the pointers */
			update_interior<interior_type, UpdateParFlag>(
				this->tags[idx].first, left, right);
			if constexpr (is_pointer_to_node<ReturnType>) {
				return this->tags[idx].first;
			} else { /* if only update pointer, then avoid */
				 /* update box */
				return node_box_type(this->tags[idx].first,
						     box_type());
			}
		} else if constexpr (ut ==
				     update_pointer_box) { /* update pointer and
							    */
							   /* box */
			update_interior<interior_type, UpdateParFlag>(
				this->tags[idx].first, left, right);
			return node_box_type(
				this->tags[idx].first,
				base_type::get_box(left.second, right.second));
		} else if constexpr (ut ==
				     tag_rebuild_node) { /* retag and save */
							 /* box for rebuild */
			update_interior<interior_type, UpdateParFlag>(
				this->tags[idx].first, left, right);
			auto new_box =
				base_type::get_box(left.second, right.second);
			if (this->tags[idx].second ==
			    base_type::bucket_num + 3) {
				func(new_box, idx);
			}
			return node_box_type(this->tags[idx].first,
					     std::move(new_box));
		} else if constexpr (ut == post_del_update) { /* avoid update
								 pointers */
			if (!func(1)) { /* query whether under the rebuild_tree
					 */
				update_interior<interior_type, UpdateParFlag>(
					this->tags[idx].first, left, right);
				return node_box_type(
					this->tags[idx].first,
					box_type()); /* box computed before */
			} else if (this->tags[idx].second ==
				   bucket_num + 3) { /* recurse back */
				func(0); /* disable the under_rebuild_tree flag
					  */
				if (!this->tags[idx].first->is_leaf) {
					static_cast<interior_type *>(
						this->tags[idx].first)
						->reset_parallel_flag();
				}
				assert(func(1) == false);
				return node_box_type(this->tags[idx].first,
						     box_type());
			} else { /* the tree has been deleted */
				return node_box_type(nullptr, box_type());
			}
		} else {
		}
	}

	/* update inner tree for multi nodes */
	template <update_type ut, bool UpdateParFlag, typename ReturnType,
		  typename Func>
		requires is_multi_node<interior_type>
	ReturnType update_inner_tree_recursive(
		bucket_type idx, parlay::sequence<ReturnType> const &tree_nodes,
		bucket_type &p, Func &&func)
	{
		if (tags[idx].second == base_type::bucket_num + 1 ||
		    tags[idx].second == base_type::bucket_num + 2) {
			return tree_nodes[p++];
		}

		if constexpr (ut == post_del_update) {
			if (this->tags[idx].second == bucket_num + 3) {
				func(0); /* close the under_rebuild_tree flag */
				assert(func(1) == true);
			}
		}

		typename interior_type::node_arr_type new_nodes;
		for (bucket_type i = 0; i < interior_type::get_regions(); ++i) {
			new_nodes[i] =
				update_inner_tree_recursive<ut, UpdateParFlag>(
					idx * interior_type::get_regions() + i,
					tree_nodes, p, func);
		}

		if constexpr (ut == update_pointer) {
			base_type::template update_interior<interior_type,
							    UpdateParFlag>(
				tags[idx].first, new_nodes);
			return this->tags[idx].first;
		} else if constexpr (ut == tag_rebuild_node) {
			update_interior<interior_type, UpdateParFlag>(
				this->tags[idx].first, new_nodes);
			if (this->tags[idx].second ==
			    base_type::bucket_num + 3) {
				func(idx);
			}
			return this->tags[idx].first;
		} else if constexpr (ut == post_del_update) { /* avoid update */
							      /*
							       * pointers for deleted
							       * trees
							       */
			if (!func(1)) { /* not under rebuild tree */
				update_interior<interior_type, UpdateParFlag>(
					this->tags[idx].first, new_nodes);
				return this->tags[idx]
					.first; /* box has been computed before
						 */
			} else if (this->tags[idx].second ==
				   bucket_num + 3) { /* back */
				func(0); /* disable the under_rebuild_tree flag
					  */
				assert(func(1) == false);
				if (!this->tags[idx].first->is_leaf) {
					static_cast<interior_type *>(
						this->tags[idx].first)
						->reset_parallel_flag();
				}
				return this->tags[idx].first;
			} else { /* the tree has been deleted */
				return nullptr;
			}
		} else {
		}
	}

	/* variables */
	bucket_type tags_num;
	node_tag_seq_type tags; /* PARA: Assign each node a tag, aka skeleton */
	parlay::sequence<balls_type> sums_tree;
	bucket_seq_type
		rev_tag; /* PARA: maps tag to the position in skeleton */
	parlay::sequence<balls_type> sums;
};
}; /* namespace psi */

#endif /* PSI_BASE_TREE_IMPL_INNER_TREE_HPP_ */
