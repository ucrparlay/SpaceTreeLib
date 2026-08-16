#ifndef PSI_KD_TREE_IMPL_KD_BUILD_TREE_HPP_
#define PSI_KD_TREE_IMPL_KD_BUILD_TREE_HPP_

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include "psi/kd_tree.h"
#include "psi/dependence/tree_node.h"

namespace psi
{
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
void kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	     ImbaRatio>::build(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { build_(A); });
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	     ImbaRatio>::divide_rotate(slice_type in, splitter_seq_type &pivots,
				       dims_type dim, bucket_type idx,
				       box_seq_type &box_seq,
				       box_type const &box)
{
	if (idx > base_type::pivot_num) {
		/*
		 * sometimes cut dimension can be -1
		 * never use pivots[idx].first to check whether it is in
		 * bucket; instead, use idx > PIVOT_NUM
		 */
		box_seq[idx - base_type::bucket_num] = box;
		pivots[idx] = splitter_type(0, idx - base_type::bucket_num);
		return;
	}
	size_t n = in.size();
	dims_type cutting_dim = split_rule_.find_cutting_dimension(box, dim);
	assert(cutting_dim < base_type::num_dims);

	pivots[idx] = split_rule_.split_sample(in, cutting_dim, box);

	box_cut_type box_cut(box, pivots[idx], true);

	cutting_dim = split_rule_.next_dimension(cutting_dim);
	divide_rotate(in.cut(0, n / 2), pivots, cutting_dim, 2 * idx, box_seq,
		      box_cut.get_first_box_cut());
	divide_rotate(in.cut(n / 2, n), pivots, cutting_dim, 2 * idx + 1,
		      box_seq, box_cut.get_second_box_cut());
	return;
}

/* starting at dimesion dim and pick pivots in a rotation manner */
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	     ImbaRatio>::pick_pivots(slice_type in, size_t const &n,
				     splitter_seq_type &pivots,
				     dims_type const dim, box_seq_type &box_seq,
				     box_type const &bx)
{
	size_t size =
		std::min(n, static_cast<size_t>(32 * base_type::bucket_num));
	assert(size <= n);

	points_type arr = points_type::uninitialized(size);
	base_type::sample_points(in, arr);

	/* pick pivots */
	divide_rotate(arr.cut(0, size), pivots, dim, 1, box_seq, bx);
	return;
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
node *kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	      ImbaRatio>::serial_build_recursive(slice_type in, slice_type out,
						 dims_type dim,
						 box_type const &box)
{
	size_t n = in.size();

	if (n == 0)
		return alloc_empty_leaf_node<slice_type, leaf_type>();

	if (n <= base_type::leaf_capacity)
		return alloc_normal_leaf_node<slice_type, leaf_type>(in);

	dims_type d = split_rule_.find_cutting_dimension(box, dim);
	auto [split_iter, split] = split_rule_.split_input(in, d, box);

	if (!split.has_value()) { /* split fails */
		if (in.end() ==
		    std::ranges::find_if_not(in, [&](Point const &p) {
			    return p.same_dimension(in[0]);
		    })) { /* check whether all elements are identical */
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
				return alloc_dummy_leaf_node<slice_type,
							     leaf_type>(in);
			}
		} else { /* current dim d is same but other dims are not */
			return split_rule_.handle_undivided(*this, in, out, box,
							    dim);
		}
	}

	assert(std::ranges::all_of(in.begin(), split_iter, [&](Point &p) {
		return num_type::lt(p.pnt[split.value().second],
				    split.value().first);
	}));
	assert(std::ranges::all_of(split_iter, in.end(), [&](Point &p) {
		return num_type::geq(p.pnt[split.value().second],
				     split.value().first);
	}));

	box_cut_type box_cut(box, split.value(), true);

	d = split_rule_.next_dimension(d);
	node *L, *R;

	L = serial_build_recursive(in.cut(0, split_iter - in.begin()),
				   out.cut(0, split_iter - in.begin()), d,
				   box_cut.get_first_box_cut());
	R = serial_build_recursive(in.cut(split_iter - in.begin(), n),
				   out.cut(split_iter - in.begin(), n), d,
				   box_cut.get_second_box_cut());
	return alloc_interior_node<interior_type>(L, R, split.value());
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
node *kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	      ImbaRatio>::build_recursive(slice_type in, slice_type out,
					  dims_type dim, box_type const &bx)
{
	assert(in.size() == 0 ||
	       base_type::within_box(base_type::get_box(in), bx));

	if (in.size() <= base_type::serial_build_cutoff) {
		return serial_build_recursive(in, out, dim, bx);
	}

	auto pivots = splitter_seq_type::uninitialized(
		base_type::pivot_num + base_type::bucket_num + 1);
	auto box_seq = box_seq_type::uninitialized(base_type::bucket_num);
	parlay::sequence<balls_type> sums;

	pick_pivots(in, in.size(), pivots, dim, box_seq, bx);
	base_type::partition(in, out, in.size(), pivots, sums);

	/*
	 * if random sampling failed to split points, re-partitions using
	 * serail approach
	 */
	auto tree_nodes =
		parlay::sequence<node *>::uninitialized(base_type::bucket_num);
	auto nodes_map = bucket_seq_type::uninitialized(base_type::bucket_num);
	bucket_type zeros = std::ranges::count(sums, 0), cnt = 0;

	if (zeros == base_type::bucket_num - 1) { /* switch to seral */
		/*
		 * TODO: add parallelsim within this call
		 * see parallel kth element
		 */
		return serial_build_recursive(in, out, dim, bx);
	}

	/* alloc empty leaf beforehand to avoid spawn threads */
	for (bucket_type i = 0; i < base_type::bucket_num; ++i) {
		if (!sums[i]) {
			tree_nodes[i] =
				alloc_empty_leaf_node<slice_type, leaf_type>();
		} else {
			nodes_map[cnt++] = i;
		}
	}

	for (bucket_type i = 0; i < base_type::build_depth_once; ++i) {
		dim = split_rule_.next_dimension(dim);
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
				in.cut(start, start + sums[nodes_map[i]]), dim,
				box_seq[nodes_map[i]]);
		},
		1);

	return base_type::template build_inner_tree<leaf_type, interior_type>(
		1, pivots, tree_nodes);
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	     ImbaRatio>::build_(slice_type A)
{
	base_type::template delete_tree_nodes<leaf_type, interior_type>();
	points_type B = points_type::uninitialized(A.size());
	this->tree_box_ = base_type::get_box(A);
	this->root_ =
		build_recursive(A, B.cut(0, A.size()), 0, this->tree_box_);
	assert(this->root_ != nullptr);
	return;
}

} /* namespace psi */

#endif /* PSI_KD_TREE_IMPL_KD_BUILD_TREE_HPP_ */
