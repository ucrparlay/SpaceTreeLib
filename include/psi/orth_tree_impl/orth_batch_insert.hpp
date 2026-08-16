#ifndef PSI_ORTH_BATCH_INSERT_HPP_
#define PSI_ORTH_BATCH_INSERT_HPP_

#include <utility>

#include "psi/orth_tree.h"

namespace psi
{

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::batch_insert(Range &&in)
{
	static_assert(base_type::build_depth_once % md == 0);
	assert(md == base_type::num_dims);
	// TODO: handling the case that insert box is no in the tree box
	assert(base_type::within_box(base_type::get_box(in), this->tree_box_));
	base_type::ingest_range(in, [&](slice_type A) { batch_insert_(A); });
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::batch_insert_(slice_type A)
{
	if (this->root_ == nullptr) { // TODO: may check using explicity tag
		return build_(A);
	}

	points_type B = points_type::uninitialized(A.size());
	node *T = this->root_;
	// this->tree_box_ = base_type::get_box(this->tree_box_,
	// base_type::get_box(A)); PERF: no need to compute bounding box here,
	// checked previously
	this->root_ = batch_insert_recursive(T, A, B.cut(0, A.size()),
					     this->tree_box_);
	assert(this->root_ != NULL);
	return;
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::serial_split_skeleton(node *T, slice_type in,
						 dims_type dim, dims_type idx,
						 parlay::sequence<balls_type>
							 &sums)
{
	// TODO: change it using the split_rule_.split()
	if (dim == base_type::num_dims) {
		sums[idx - node_regions] = in.size();
		return;
	}

	auto mid = static_cast<interior_type *>(T)->split[dim].first;
	assert(dim == static_cast<interior_type *>(T)->split[dim].second);

	points_iter_type split_iter =
		std::ranges::partition(in, [&](Point const &p) {
			return num_type::lt(p.pnt[dim], mid);
		}).begin();

	serial_split_skeleton(T, in.cut(0, split_iter - in.begin()), dim + 1,
			      2 * idx, sums);
	serial_split_skeleton(T, in.cut(split_iter - in.begin(), in.size()),
			      dim + 1, 2 * idx + 1, sums);
	return;
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
node *orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
		ImbaRatio>::batch_insert_recursive(node *T, slice_type in,
						   slice_type out,
						   box_type const &box)
{
	size_t n = in.size();

	if (n == 0)
		return T;

	if (T->is_leaf) {
		leaf_type *tl = static_cast<leaf_type *>(T);
		// NOTE: insert the points to normal leaf if the capacity
		// allows; or check if the leaf is dummy and contains same
		// points as inputs
		if ((!tl->is_dummy &&
		     n + T->size <= base_type::leaf_capacity) ||
		    (tl->is_dummy && parlay::all_of(in, [&](Point const &p) {
			     return p == tl->pts[0];
		     }))) {
			return base_type::template insert_points2_leaf<
				leaf_type>(T, in);
		} else {
			return base_type::template rebuild_with_insert<
				leaf_type, interior_type>(
				T,
				[&](node *, points_type &, points_type &) {
					return box;
				},
				in);
		}
	}

	// if (n) {
	if (n <= base_type::serial_build_cutoff) {
		parlay::sequence<balls_type> sums(node_regions, 0);
		serial_split_skeleton(T, in, 0, 1, sums);
		assert(std::cmp_equal(
			std::accumulate(sums.begin(), sums.end(), 0), n));

		auto ti = static_cast<interior_type *>(T);
		node_arr_type new_nodes;
		size_t start = 0;
		for (dims_type i = 0; i < node_regions; ++i) {
			new_nodes[i] = batch_insert_recursive(
				ti->tree_nodes[i],
				in.cut(start, start + sums[i]),
				out.cut(start, start + sums[i]),
				ti->get_box_by_region_id(i, box));
			start += sums[i];
		}
		base_type::template update_interior<interior_type>(T,
								   new_nodes);
		assert(T->is_leaf == false);
		return T;
	}

	// NOTE: assign each node a tag
	inner_tree IT;
	IT.assign_node_tag(T, 1);
	assert(IT.tags_num > 0 && IT.tags_num <= base_type::bucket_num);

	base_type::template sieve_points<interior_type>(in, out, n, IT.tags,
							IT.sums, IT.tags_num);

	// NOTE: no need to tag imbalance node in orth tree as it never rebuilds
	// the tree, used to remap the bucket node tag to bucket_num+1 and
	// compute the bounding boxes NOTE: we pass has_tomb as true, to make
	// the leaf set to bucket_num+1 IT.tag_imbalance_node([]() -> bool {
	// return false; });
	box_seq_type box_seq(IT.tags_num); // PARA: the box for bucket nodes
	[[maybe_unused]] auto [re_num, tot_re_size] =
		IT.template tag_imbalance_node_deletion<false>(
			box_seq, box, true,
			[&](bucket_type) -> bool { return false; });

	auto tree_nodes = parlay::sequence<node *>::uninitialized(IT.tags_num);

	parlay::parallel_for(
		0, IT.tags_num,
		[&](decltype(IT.tags_num) i) {
			size_t s = 0;
			for (decltype(IT.tags_num) j = 0; j < i; j++) {
				s += IT.sums_tree[IT.rev_tag[j]];
			}

			assert(IT.tags[IT.rev_tag[i]].second ==
			       base_type::bucket_num + 1);
			tree_nodes[i] = batch_insert_recursive(
				IT.tags[IT.rev_tag[i]].first,
				out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]),
				in.cut(s, s + IT.sums_tree[IT.rev_tag[i]]),
				box_seq[i]);
		},
		1);

	return IT.template update_inner_tree<inner_tree::kUpdatePointer>(
		tree_nodes);
}
} // namespace psi

#endif // PSI_ORTH_BATCH_INSERT_HPP_
