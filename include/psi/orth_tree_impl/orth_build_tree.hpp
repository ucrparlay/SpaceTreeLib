#ifndef PSI_ORTH_TREE_IMPL_ORTH_BUILD_TREE_HPP_
#define PSI_ORTH_TREE_IMPL_ORTH_BUILD_TREE_HPP_

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include <algorithm>

#include "psi/dependence/tree_node.h"
#include "psi/orth_tree.h"

namespace psi
{
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range, typename... Args>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::build(Range &&in, Args &&...args)
{
	static_assert(base_type::build_depth_once % md == 0);
	base_type::ingest_range(in, [&](slice_type A) {
		build_(A, std::forward<Args>(args)...);
	});
}

// TODO: maybe we don't need this function, it can be directly computed by value
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::divide_rotate(hyper_plane_seq_type &pivots,
					 dims_type dim, bucket_type idx,
					 box_seq_type &box_seq,
					 box_type const &box)
{
	if (idx > base_type::pivot_num) {
		// WARN: sometimes cut dimension can be -1, never use
		// pivots[idx].first == -1 to check whether it is in bucket;
		// instead, use idx > PIVOT_NUM
		box_seq[idx - base_type::bucket_num] = box;
		pivots[idx] = hyper_plane_type(0, idx - base_type::bucket_num);
		return;
	}

	pivots[idx] = split_rule_.split_sample(slice_type(nullptr, nullptr),
					       dim, box);

	box_cut_type box_cut(box, pivots[idx], true);
	// dim = (dim + 1) % base_type::num_dims;
	dim = split_rule_.next_dimension(dim);
	divide_rotate(pivots, dim, 2 * idx, box_seq,
		      box_cut.get_first_box_cut());
	divide_rotate(pivots, dim, 2 * idx + 1, box_seq,
		      box_cut.get_second_box_cut());

	return;
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::serial_split(slice_type in, dims_type dim,
					dims_type idx, box_type const &box,
					parlay::sequence<balls_type> &sums)
{
	assert(dim <= base_type::num_dims);

	if (dim == base_type::num_dims) {
		sums[idx - node_regions] = in.size();
		return;
	}

	auto [split_iter, split] = split_rule_.split_input(in, dim, box);

	serial_split(in.cut(0, split_iter - in.begin()), dim + 1, idx << 1, box,
		     sums);
	serial_split(in.cut(split_iter - in.begin(), in.size()), dim + 1,
		     idx << 1 | 1, box, sums);
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
node *orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
		ImbaRatio>::serial_build_recursive(slice_type in,
						   slice_type out,
						   box_type const &box,
						   bool checked_duplicate)
{
	assert(in.size() == 0 ||
	       base_type::within_box(base_type::get_box(in), box));
	size_t n = in.size();

	if (n == 0) {
		return alloc_empty_leaf_node<slice_type, leaf_type>();
	}

	if (n <= base_type::leaf_capacity) {
		return alloc_normal_leaf_node<slice_type, leaf_type>(in);
	}

	assert(splitter_num == base_type::num_dims);

	dims_type dim = 0;
	parlay::sequence<balls_type> sums(node_regions, 0);
	auto splitter = interior_type::compute_splitter(box);
	serial_split(in, dim, 1, box, sums);
	assert(std::cmp_equal(std::accumulate(sums.begin(), sums.end(), 0), n));

	if (std::ranges::count(sums, 0) == node_regions - 1) { // split fails
		if (std::ranges::find_if_not(
			    in, [&](Point const &p) { // early return
				    return p.same_dimension(in[0]);
			    }) == in.end()) {
			if constexpr (is_aug_point<Point>) {
				if constexpr (
					Point::is_non_trivial_augmentation()) {
					return alloc_fix_size_leaf_node<
						slice_type, leaf_type>(
						in,
						std::max(
							in.size(),
							static_cast<size_t>(
								base_type::
									leaf_capacity)));
				} else {
					return alloc_dummy_leaf_node<slice_type,
								     leaf_type>(
						in);
				}
			} else {
				// WARN: Need to pass full range, since it needs
				// to compute the size
				return alloc_dummy_leaf_node<slice_type,
							     leaf_type>(in);
			}
		} else {
			return split_rule_.handle_undivided(*this, in, out,
							    box);
		}
	}

	orth_node_arr_type tree_nodes;
	size_t start = 0;
	for (dims_type i = 0; i < node_regions; ++i) {
		// NOTE: iterate through non-empty partitions, put them into the
		// position identified by non_empty_node
		tree_nodes[i] = serial_build_recursive(
			in.cut(start, start + sums[i]),
			out.cut(start, start + sums[i]),
			interior_type::get_box_by_region_id(i, splitter, box),
			checked_duplicate);
		start += sums[i];
	}

	return alloc_interior_node<interior_type>(tree_nodes, splitter);
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
node *orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
		ImbaRatio>::build_recursive(slice_type in, slice_type out,
					    box_type const &box)
{
	// TODO: may ensure the bucket is corresponding the the splitter
	assert(in.size() == 0 ||
	       base_type::within_box(base_type::get_box(in), box));
	size_t n = in.size();

	// if (in.size()) {
	if (n <= base_type::serial_build_cutoff) {
		return serial_build_recursive(in, out, box, false);
	}

	auto pivots = hyper_plane_seq_type::uninitialized(
		base_type::pivot_num + base_type::bucket_num + 1);
	auto box_seq = box_seq_type::uninitialized(base_type::bucket_num);
	parlay::sequence<balls_type> sums;

	divide_rotate(pivots, 0, 1, box_seq, box);
	base_type::partition(in, out, in.size(), pivots, sums);

	auto tree_nodes =
		parlay::sequence<node *>::uninitialized(base_type::bucket_num);
	auto nodes_map = bucket_seq_type::uninitialized(base_type::bucket_num);
	bucket_type zeros = std::ranges::count(sums, 0), cnt = 0;

	if (zeros == base_type::bucket_num - 1) { // NOTE: switch to seral
		// TODO: add parallelsim within this call
		// see parallel kth element
		return serial_build_recursive(in, out, box, false);
	}

	for (bucket_type i = 0; i < base_type::bucket_num; ++i) {
		if (sums[i] == 0) {
			tree_nodes[i] =
				alloc_empty_leaf_node<slice_type, leaf_type>();
		} else {
			nodes_map[cnt++] = i;
		}
	}

	parlay::parallel_for(
		0, base_type::bucket_num - zeros,
		[&](bucket_type i) {
			size_t start = 0;
			for (bucket_type j = 0; j < nodes_map[i]; ++j) {
				start += sums[j];
			}

			tree_nodes[nodes_map[i]] = build_recursive(
				out.cut(start, start + sums[nodes_map[i]]),
				in.cut(start, start + sums[nodes_map[i]]),
				box_seq[nodes_map[i]]);
		},
		1);

	return base_type::template build_inner_tree<interior_type>(1, pivots,
								   tree_nodes);
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::build_(slice_type A)
{
	base_type::template delete_tree_nodes<leaf_type, interior_type>();
	points_type B = points_type::uninitialized(A.size());
	if (!fixed_box) {
		this->tree_box_ = base_type::get_box(A);
	}
	this->root_ = build_recursive(A, B.cut(0, A.size()), this->tree_box_);
	assert(this->root_ != nullptr);
	return;
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::build_(slice_type A, box_type const &box)
{
	assert(base_type::within_box(base_type::get_box(A), box));

	base_type::template delete_tree_nodes<leaf_type, interior_type>();
	points_type B = points_type::uninitialized(A.size());
	if (!fixed_box) {
		this->tree_box_ = box;
	}
	this->root_ = build_recursive(A, B.cut(0, A.size()), this->tree_box_);
	assert(this->root_ != nullptr);
	return;
}

} // namespace psi

#endif // PSI_ORTH_TREE_IMPL_ORTH_BUILD_TREE_HPP_
