#ifndef PSI_ORTH_TREE_IMPL_ORTH_BATCH_DIFF_HPP_
#define PSI_ORTH_TREE_IMPL_ORTH_BATCH_DIFF_HPP_

#include <tuple>

#include "psi/orth_tree.h"

namespace psi
{

// NOTE: default batch delete
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::batch_diff(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { batch_diff_(A); });
	return;
}

// NOTE: assume points are partially covered in the tree
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::batch_diff_(slice_type A)
{
	if (this->root_ == nullptr)
		return;

	// NOTE: diff points from the tree
	points_type B = points_type::uninitialized(A.size());
	this->root_ =
		batch_diff_recursive(this->root_, A, parlay::make_slice(B));

	// NOTE: launch the rebuild
	// PARA: @prepare_func: function that computes the new parameters before
	// the rebuildtree recursive
	auto prepare_func = [&]([[maybe_unused]] node *T,
				[[maybe_unused]] size_t i,
				box_type const &box) {
		auto new_box =
			static_cast<interior_type *>(T)->get_box_by_region_id(
				i, box);
		assert(base_type::within_box(new_box, box));
		return std::make_tuple(std::move(new_box));
	};
	this->root_ = base_type::template rebuild_tree_recursive<
		leaf_type, interior_type, false>(this->root_, prepare_func,
						 split_rule_.allow_rebuild(),
						 this->tree_box_);
	return;
}

// NOTE: the orth does not need box since the box will never change
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
node *orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
		ImbaRatio>::batch_diff_recursive(node *T, slice_type in,
						 slice_type out)
{
	size_t n = in.size();

	if (n == 0) {
		return T;
	}

	if (T->is_leaf) {
		return base_type::template diff_points4_leaf<leaf_type, node *>(
			T, in);
	}

	// if (in.size()) {
	if (in.size() <= base_type::serial_build_cutoff) {
		parlay::sequence<balls_type> sums(node_regions, 0);
		serial_split_skeleton(T, in, 0, 1, sums);
		assert(std::cmp_equal(
			std::accumulate(sums.begin(), sums.end(), 0), n));

		auto ti = static_cast<interior_type *>(T);
		node_arr_type new_nodes;

		size_t start = 0;
		for (dims_type i = 0; i < node_regions; ++i) {
			new_nodes[i] = batch_diff_recursive(
				ti->tree_nodes[i],
				in.cut(start, start + sums[i]),
				out.cut(start, start + sums[i]));
			start += sums[i];
		}

		bool const force_parallel_flag =
			ti->size > base_type::serial_build_cutoff;
		base_type::template update_interior<interior_type>(T,
								   new_nodes);
		assert(T->is_leaf == false);

		if (base_type::sparse_node(0, ti->size)) {
			ti->set_parallel_flag(force_parallel_flag);
		}

		return T;
	}

	inner_tree IT;
	IT.assign_node_tag(T, 1);
	assert(IT.tags_num > 0 && IT.tags_num <= base_type::bucket_num);
	base_type::template sieve_points<interior_type>(in, out, n, IT.tags,
							IT.sums, IT.tags_num);
	IT.tag_puffy_nodes();

	// PERF: no need to call tag imbalance node here, as the bounding box
	// for orth-tree is fixed

	auto tree_nodes = parlay::sequence<node *>::uninitialized(IT.tags_num);
	parlay::parallel_for(
		0, IT.tags_num,
		[&](decltype(IT.tags_num) i) {
			size_t start = 0;
			for (decltype(IT.tags_num) j = 0; j < i; j++) {
				start += IT.sums[j];
			}

			tree_nodes[i] = batch_diff_recursive(
				IT.tags[IT.rev_tag[i]].first,
				out.cut(start, start + IT.sums[i]),
				in.cut(start, start + IT.sums[i]));

			// NOTE: after pick the tag, the tag id is same as the
			// bucket id. in order to match the base case in
			// update_inner_tree, we need to manually change it to
			// bucket_num+2, i.e., none-of its ancestor has been
			// rebuilt
			IT.tags[IT.rev_tag[i]].second =
				base_type::bucket_num + 2;
		},
		1);

	return IT.template update_inner_tree<inner_tree::kUpdatePointer, false>(
		tree_nodes);
}

} // namespace psi

#endif // PSI_ORTH_TREE_IMPL_ORTH_BATCH_DIFF_HPP_
