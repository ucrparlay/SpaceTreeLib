#ifndef PSI_ORTH_TREE_IMPL_ORTH_BATCH_DELETE_HPP_
#define PSI_ORTH_TREE_IMPL_ORTH_BATCH_DELETE_HPP_

#include "../orth_tree.h"

namespace psi
{

// NOTE: default batch delete
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::batch_delete(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { batch_delete_(A); });
	return;
}

// NOTE: assume all points_type are fully covered in the tree
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::batch_delete_(slice_type A)
{
	if (this->root_ == nullptr) {
		assert(A.size() == 0);
		return;
	}

	points_type B = points_type::uninitialized(A.size());
	this->root_ = batch_delete_recursive(
		this->root_, A, parlay::make_slice(B), this->tree_box_, 1);
	return;
}

// NOTE: delete with rebuild, with the assumption that all points_type are in
// the tree
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
node *orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
		ImbaRatio>::batch_delete_recursive(node *T, slice_type in,
						   slice_type out,
						   box_type const &box,
						   bool has_tomb)
{
	size_t n = in.size();

	if (n == 0) {
		return T;
	}

	// INFO: delete the whole tree directly if its size equals to the input
	// size, can be used to accelerate the whole deletion process
#ifndef DISABLE_BATCH_DELETE_SIZE_OPT
	if (n == T->size) {
		if (has_tomb) {
			base_type::template delete_tree_recursive<
				leaf_type, interior_type>(T);
			return alloc_empty_leaf_node<slice_type, leaf_type>();
		}
		if (!T->is_leaf) {
			auto ti = static_cast<interior_type *>(T);
			ti->set_parallel_flag(T->size >
					      base_type::serial_build_cutoff);
		}
		T->size = 0;
		return T;
	}
#endif

	if (T->is_leaf) {
		return base_type::template delete_points4_leaf<leaf_type,
							       node *>(T, in);
	}

	if (in.size() <= base_type::serial_build_cutoff) {
		parlay::sequence<balls_type> sums(node_regions, 0);
		serial_split_skeleton(T, in, 0, 1, sums);
		assert(std::cmp_equal(
			std::accumulate(sums.begin(), sums.end(), 0), n));

		bool putTomb = has_tomb &&
			       (base_type::sparse_node(in.size(), T->size));
		has_tomb = putTomb ? false : has_tomb;
		assert(putTomb ? (!has_tomb) : true);

		auto ti = static_cast<interior_type *>(T);
		orth_node_arr_type new_nodes;

		size_t start = 0;
		for (dims_type i = 0; i < node_regions; ++i) {
			new_nodes[i] = batch_delete_recursive(
				ti->tree_nodes[i],
				in.cut(start, start + sums[i]),
				out.cut(start, start + sums[i]),
				ti->get_box_by_region_id(i, box), has_tomb);
			start += sums[i];
		}

		bool const force_parallel_flag =
			ti->size > base_type::serial_build_cutoff;

		base_type::template update_interior<interior_type>(T,
								   new_nodes);
		assert(T->is_leaf == false);

		if (!has_tomb) {
			ti->set_parallel_flag(force_parallel_flag);
		}

		if (putTomb) {
			// NOTE: the box is the one that passed from the top,
			// which should be the correct one associated with this
			// node
			assert(T->size <= base_type::leaf_capacity);
			assert(base_type::within_box(
				base_type::template get_box<leaf_type,
							    interior_type>(T),
				box));
			return base_type::template rebuild_single_tree<
				leaf_type, interior_type, false>(T, box);
		}
		return T;
	}

	inner_tree IT;
	IT.assign_node_tag(T, 1);
	assert(IT.tags_num > 0 && IT.tags_num <= base_type::bucket_num);
	base_type::template sieve_points<interior_type>(in, out, n, IT.tags,
							IT.sums, IT.tags_num);

	box_seq_type box_seq(IT.tags_num); // PARA: the box for bucket nodes
	[[maybe_unused]] auto [re_num, tot_re_size] =
		IT.template tag_imbalance_node_deletion<true>(
			box_seq, box, has_tomb, [&](bucket_type idx) -> bool {
				// NOTE: only the sparse node will be rebuilt
				return base_type::sparse_node(
					IT.sums_tree[idx],
					IT.tags[idx].first->size);
			});

	assert(re_num <= IT.tags_num);

	// NOTE: continue seieving points to leaf first
	auto tree_nodes = parlay::sequence<node *>::uninitialized(IT.tags_num);
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

			tree_nodes[i] = batch_delete_recursive(
				IT.tags[IT.rev_tag[i]].first,
				out.cut(start, start + IT.sums[i]),
				in.cut(start, start + IT.sums[i]), box_seq[i],
				IT.tags[IT.rev_tag[i]].second ==
					base_type::bucket_num + 1);
		},
		1);

	// NOTE: handling of rebuild
	// WARN: the rebuild node is on top
	// NOTE: retag the inba-nodes and save the bounding boxes
	[[maybe_unused]] node *new_node =
		IT.template update_inner_tree<inner_tree::kTagRebuildNode>(
			tree_nodes);
	assert(IT.tags_num == re_num);

	parlay::parallel_for(0, IT.tags_num, [&](size_t i) {
		assert(IT.tags[IT.rev_tag[i]].second ==
		       base_type::bucket_num + 3);

		if (IT.tags[IT.rev_tag[i]].first->size == 0) { // NOTE: empty

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
		} else { // NOTE: rebuild
			assert(base_type::within_box(
				base_type::template get_box<leaf_type,
							    interior_type>(
					IT.tags[IT.rev_tag[i]].first),
				IT.get_box_by_region_idx(IT.rev_tag[i], box)));
			IT.tags[IT.rev_tag[i]].first =
				base_type::template rebuild_single_tree<
					leaf_type, interior_type, false>(
					IT.tags[IT.rev_tag[i]].first,
					IT.get_box_by_region_idx(IT.rev_tag[i],
								 box));
		}
	}); // PERF: allow the parlay decide the granularity to accelerate the
	    // small tree rebuild

	return IT.template update_inner_tree<inner_tree::kPostDelUpdate>(
		tree_nodes);
}

} // namespace psi

#endif // PSI_ORTH_TREE_IMPL_ORTH_BATCH_DELETE_HPP_
