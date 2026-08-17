#ifndef PSI_KD_TREE_IMPL_KD_BATCH_DELETE_HPP
#define PSI_KD_TREE_IMPL_KD_BATCH_DELETE_HPP

#include "psi/kd_tree.h"

namespace psi
{

/* default batch delete */
template <typename Traits>
template <typename Range>
void kd_tree<Traits>::batch_delete(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { batch_delete_(A); });
	return;
}

/* assume all points_type are fully covered in the tree */
template <typename Traits>
void kd_tree<Traits>::batch_delete_(slice_type A)
{
	if (this->root_ == nullptr) {
		assert(A.size() == 0);
		return;
	}

	points_type B = points_type::uninitialized(A.size());
	node *T = this->root_;
	box_type box = this->tree_box_;
	dims_type d =
		T->is_leaf ? 0 : static_cast<interior_type *>(T)->split.second;
	std::tie(this->root_, this->tree_box_) =
		batch_delete_recursive(T, box, A, parlay::make_slice(B), d, 1);
	return;
}

/*
 * delete with rebuild, with the assumption that all points_type are in
 * the tree the param d can be only used when rotate cutting is applied
 */
template <typename Traits>
typename kd_tree<Traits>::node_box_type kd_tree<Traits>::batch_delete_recursive(
	node *T, typename kd_tree<Traits>::box_type const &box, slice_type in,
	slice_type out, dims_type d, bool has_tomb)
{
	size_t n = in.size();

	if (n == 0) {
		if constexpr (has_box<typename interior_type::at_type> &&
			      has_box<typename leaf_type::at_type>) {
			assert(base_type::same_box(
				base_type::template get_box<leaf_type,
							    interior_type>(T),
				retrieve_box<leaf_type, interior_type>(T)));
			return node_box_type(
				T, retrieve_box<leaf_type, interior_type>(T));
		} else {
			assert(base_type::within_box(
				base_type::template get_box<leaf_type,
							    interior_type>(T),
				box));
			return node_box_type(T, box);
		}
	}

#ifndef DISABLE_BATCH_DELETE_SIZE_OPT
	/* INFO: it can be used to accelerate the whole deletion process */
	if (n == T->size) {
		if (has_tomb) { /* rebuild this subtree */
			base_type::template delete_tree_recursive<
				leaf_type, interior_type>(T);
			return node_box_type(
				alloc_empty_leaf_node<slice_type, leaf_type>(),
				base_type::get_empty_box());
		}
		/* within a rebuild tree */
		if (!T->is_leaf) { /* interior */
			auto ti = static_cast<interior_type *>(T);
			ti->reset_aug(); /* needs to put before set parallel */
					 /* flag */
			/*
			 * only set the flag for root, the remaining tree
			 * is still unset
			 */
			ti->set_parallel_flag(T->size >
					      base_type::serial_build_cutoff);
		} else { /* leaf */
			auto tl = static_cast<leaf_type *>(T);
			tl->reset_aug();
		}
		T->size = 0;
		return node_box_type(T, base_type::get_empty_box());
	}
#endif

	if (T->is_leaf) {
		return base_type::template delete_points4_leaf<leaf_type,
							       node_box_type>(
			T, in);
	}

	if (in.size() <= base_type::serial_build_cutoff) {
		interior_type *ti = static_cast<interior_type *>(T);
		points_iter_type split_iter =
			std::ranges::partition(in, [&](point_type const &p) {
				return num_type::lt(p.pnt[ti->split.second],
						    ti->split.first);
			}).begin();

		/*
		 * put the tomb if the remaining points number are below
		 * sparse_leaf_threshold (to avoid next insertion exceeds the
		 * limit) or imbalance
		 */
		bool put_tomb =
			has_tomb &&
			(base_type::sparse_node(in.size(), ti->size) ||
			 (split_rule_.allow_rebuild() &&
			  base_type::imbalance_node(
				  ti->left->size - (split_iter - in.begin()),
				  ti->size - in.size())));
		has_tomb = put_tomb ? false : has_tomb;
		assert(put_tomb ? (!has_tomb) : true);

		dims_type next_dim = split_rule_.next_dimension(d);
		box_cut_type box_cut(box, ti->split, true);

		auto [L, lbox] = batch_delete_recursive(
			ti->left, box_cut.get_first_box_cut(),
			in.cut(0, split_iter - in.begin()),
			out.cut(0, split_iter - in.begin()), next_dim,
			has_tomb);
		auto [R, rbox] = batch_delete_recursive(
			ti->right, box_cut.get_second_box_cut(),
			in.cut(split_iter - in.begin(), n),
			out.cut(split_iter - in.begin(), n), next_dim,
			has_tomb);

		bool const force_parallel_flag =
			ti->size > base_type::serial_build_cutoff;

		base_type::template update_interior<interior_type>(T, L, R);
		ti = static_cast<interior_type *>(T);
		assert(T->size == L->size + R->size && ti->split.second >= 0 &&
		       ti->is_leaf == false);

		if (!has_tomb) {
			ti->set_parallel_flag(force_parallel_flag);
		}

		/* rebuild */
		if (put_tomb) {
			assert(base_type::sparse_node(0, ti->size) ||
			       (split_rule_.allow_rebuild() &&
				base_type::imbalance_node(ti->left->size,
							  ti->size)));
			auto const new_box = base_type::get_box(lbox, rbox);

			assert(base_type::within_box(
				base_type::template get_box<leaf_type,
							    interior_type>(T),
				new_box));
			return node_box_type(
				base_type::template rebuild_single_tree<
					leaf_type, interior_type, false>(
					T, d, new_box),
				new_box);
		}

		return node_box_type(T, base_type::get_box(lbox, rbox));
	}

	inner_tree IT;
	IT.assign_node_tag(T, 1);
	assert(IT.tags_num > 0 && IT.tags_num <= base_type::bucket_num);
	base_type::template sieve_points<interior_type>(in, out, n, IT.tags,
							IT.sums, IT.tags_num);

	auto tree_nodes = node_box_seq_type::uninitialized(IT.tags_num);
	auto box_seq = parlay::sequence<box_type>::uninitialized(IT.tags_num);

	/* enable the force parallel flag in batch deletion */
	/* Results are used by asserts only, but the call tags the skeleton. */
	[[maybe_unused]] auto [re_num, tot_re_size] =
		IT.template tag_imbalance_node_deletion<true>(
			box_seq, box, has_tomb, [&](bucket_type idx) -> bool {
				interior_type *ti =
					static_cast<interior_type *>(
						IT.tags[idx].first);
				return base_type::sparse_node(IT.sums_tree[idx],
							      ti->size) ||
				       (split_rule_.allow_rebuild() &&
					base_type::imbalance_node(
						ti->left->size -
							IT.sums_tree[idx << 1],
						ti->size - IT.sums_tree[idx]));
			});

	assert(re_num <= IT.tags_num);

	/* delete the points in the tree */
	parlay::parallel_for(
		0, IT.tags_num,
		[&](decltype(IT.tags_num) i) {
			size_t start = 0;
			for (decltype(IT.tags_num) j = 0; j < i; j++) {
				start += IT.sums[j];
			}

			assert(IT.sums_tree[IT.rev_tag[i]] == IT.sums[i]);
			assert(IT.tags[IT.rev_tag[i]].first->size >=
			       IT.sums[i]);
			assert(base_type::within_box(
				base_type::get_box(
					out.cut(start, start + IT.sums[i])),
				base_type::template get_box<leaf_type,
							    interior_type>(
					IT.tags[IT.rev_tag[i]].first)));

			dims_type next_dim = d,
				  depth = IT.get_depth_by_index(IT.rev_tag[i]);
			for (bucket_type i = 0; i < depth; i++) {
				next_dim = split_rule_.next_dimension(next_dim);
			}

			tree_nodes[i] = batch_delete_recursive(
				IT.tags[IT.rev_tag[i]].first, box_seq[i],
				out.cut(start, start + IT.sums[i]),
				in.cut(start, start + IT.sums[i]), next_dim,
				IT.tags[IT.rev_tag[i]].second ==
					base_type::bucket_num + 1);
		},
		1);

	/*
	 * handling of rebuild (in parallel)
	 * get new box for skeleton root and rebuild nodes
	 */
	box_type const new_box = std::get<1>(
		IT.template update_inner_tree<inner_tree::tag_rebuild_node>(
			tree_nodes, box_seq));
	assert(IT.tags_num == re_num);

	/* launch the rebuild in parallel */
	parlay::parallel_for(0, IT.tags_num, [&](size_t i) {
		assert(IT.tags[IT.rev_tag[i]].second ==
		       base_type::bucket_num + 3);

		if (IT.tags[IT.rev_tag[i]].first->size == 0) { /* empty */
#ifndef DISABLE_BATCH_DELETE_SIZE_OPT
			base_type::template delete_tree_recursive<
				leaf_type, interior_type, false>(
				IT.tags[IT.rev_tag[i]].first);
#else
      base_type::template delete_tree_recursive<leaf_type, interior_type, true>(
          IT.tags[IT.rev_tag[i]].first);
#endif
			IT.tags[IT.rev_tag[i]].first =
				alloc_empty_leaf_node<slice_type, leaf_type>();
		} else { /* rebuild */
			assert(base_type::within_box(
				base_type::template get_box<leaf_type,
							    interior_type>(
					IT.tags[IT.rev_tag[i]].first),
				box_seq[i]));

			dims_type next_dim = d,
				  depth = IT.get_depth_by_index(IT.rev_tag[i]);
			for (bucket_type i = 0; i < depth; i++) {
				next_dim = split_rule_.next_dimension(next_dim);
			}
			IT.tags[IT.rev_tag[i]].first =
				base_type::template rebuild_single_tree<
					leaf_type, interior_type, false>(
					IT.tags[IT.rev_tag[i]].first, next_dim,
					box_seq[i]);
		}
	}); /* let the parlay decide granularity */

	auto const new_root = std::get<0>(
		IT.template update_inner_tree<inner_tree::post_del_update>(
			tree_nodes));
	return node_box_type(new_root, new_box);
}

} /* namespace psi */

#endif /* PSI_KD_TREE_IMPL_KD_BATCH_DELETE_HPP */
