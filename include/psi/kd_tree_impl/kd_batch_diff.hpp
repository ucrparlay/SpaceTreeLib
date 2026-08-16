#ifndef PSI_KD_TREE_IMPL_KD_BATCH_DIFF_HPP_
#define PSI_KD_TREE_IMPL_KD_BATCH_DIFF_HPP_

#include "../kd_tree.h"

namespace psi
{
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
void kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	     ImbaRatio>::batch_diff(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { batch_diff_(A); });
	return;
}

// NOTE: batch delete suitable for points_type that are pratially covered in the
// tree
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	     ImbaRatio>::batch_diff_(slice_type A)
{
	if (this->root_ == nullptr)
		return;

	points_type B = points_type::uninitialized(A.size());
	node *T = this->root_;
	box_type box = this->tree_box_;

	// NOTE: diff points from the tree
	dims_type d =
		T->is_leaf ? 0 : static_cast<interior_type *>(T)->split.second;
	std::tie(T, this->tree_box_) =
		batch_diff_recursive(T, box, A, parlay::make_slice(B), d);

	// NOTE: launch rebuild to either: rebuild the imbalance tree or remove
	// the sparse node
	d = T->is_leaf ? 0 : static_cast<interior_type *>(T)->split.second;
	auto prepare_rebuild_func = [&](node *T, dims_type d,
					box_type const &box) {
		dims_type new_dim = split_rule_.next_dimension(d);
		box_cut_type box_cut(
			box, static_cast<interior_type *>(T)->split, true);
		auto left_args =
			std::make_pair(new_dim, box_cut.get_first_box_cut());
		auto right_args =
			std::make_pair(new_dim, box_cut.get_second_box_cut());
		return std::make_pair(std::move(left_args),
				      std::move(right_args));
	};

	// PERF: in the batch diff, there is no need to set the force parallel
	// flag, as the size of the tree has been updated in the first time
	// traversal, the second time only need to follow the size of the
	// current tree
	this->root_ = base_type::template rebuild_tree_recursive<
		leaf_type, interior_type, false>(
		T, prepare_rebuild_func, this->split_rule_.allow_rebuild(), d,
		this->tree_box_);

	return;
}

// NOTE: only sieve the points_type, without rebuilding the tree
// NOTE: the kdtree needs box since the box will be changed in batch diff
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
typename kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
		 ImbaRatio>::node_box_type
kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight, ImbaRatio>::
	batch_diff_recursive(
		node *T,
		typename kd_tree<Point, SplitRule, LeafAugType, InteriorAugType,
				 SkHeight, ImbaRatio>::box_type const &box,
		slice_type in, slice_type out, dims_type d)
{
	size_t n = in.size();

	if (n == 0)
		return node_box_type(T, box);

	if (T->is_leaf) {
		return base_type::template diff_points4_leaf<leaf_type,
							     node_box_type>(T,
									    in);
	}

	if (in.size() <= base_type::serial_build_cutoff) {
		// if (in.size()) {
		interior_type *ti = static_cast<interior_type *>(T);
		points_iter_type split_iter =
			std::ranges::partition(in, [&](Point const &p) {
				return num_type::lt(p.pnt[ti->split.second],
						    ti->split.first);
			}).begin();

		dims_type next_dim = split_rule_.next_dimension(d);

		box_cut_type box_cut(box, ti->split, true);
		auto [L, lbox] = batch_diff_recursive(
			ti->left, box_cut.get_first_box_cut(),
			in.cut(0, split_iter - in.begin()),
			out.cut(0, split_iter - in.begin()), next_dim);
		auto [R, rbox] = batch_diff_recursive(
			ti->right, box_cut.get_second_box_cut(),
			in.cut(split_iter - in.begin(), n),
			out.cut(split_iter - in.begin(), n), next_dim);

		bool const force_parallel_flag =
			ti->size > base_type::serial_build_cutoff;
		base_type::template update_interior<interior_type>(T, L, R);
		assert(T->size == L->size + R->size && ti->split.second >= 0 &&
		       ti->is_leaf == false);

		// TODO: replace this one by a lambda that can be pssed to
		// rebuild function as well
		if (base_type::sparse_node(0, ti->size) ||
		    (split_rule_.allow_rebuild() &&
		     base_type::imbalance_node(ti->left->size, ti->size))) {
			ti->set_parallel_flag(force_parallel_flag);
		}

		return node_box_type(T, base_type::get_box(lbox, rbox));
	}

	inner_tree IT;
	IT.assign_node_tag(T, 1);
	assert(IT.tags_num > 0 && IT.tags_num <= base_type::bucket_num);
	base_type::template sieve_points<interior_type>(in, out, n, IT.tags,
							IT.sums, IT.tags_num);

	auto tree_nodes = node_box_seq_type::uninitialized(IT.tags_num);

	IT.tag_puffy_nodes();

	parlay::parallel_for(
		0, IT.tags_num,
		// NOTE: i is the index of the tags
		[&](bucket_type i) {
			size_t start = 0;
			for (bucket_type j = 0; j < i; j++) {
				// NOTE: should have same effect as using
				// sums_tree if using sums_tree then it should
				// be sums_tree[rev_tag[j]]
				start += IT.sums[j];
			}

			dims_type next_dim = d,
				  depth = IT.get_depth_by_index(IT.rev_tag[i]);
			for (bucket_type i = 0; i < depth; i++) {
				next_dim = split_rule_.next_dimension(next_dim);
			}

			assert(base_type::within_box(
				base_type::template get_box<leaf_type,
							    interior_type>(
					IT.tags[IT.rev_tag[i]].first),
				IT.get_box_by_region_idx(IT.rev_tag[i], box)));

			tree_nodes[i] = batch_diff_recursive(
				IT.tags[IT.rev_tag[i]].first,
				IT.get_box_by_region_idx(IT.rev_tag[i], box),
				out.cut(start, start + IT.sums[i]),
				in.cut(start, start + IT.sums[i]), next_dim);
		},
		1);

	return IT.template update_inner_tree<inner_tree::kUpdatePointerBox,
					     false>(tree_nodes);
}

} // namespace psi

#endif // PSI_KD_TREE_IMPL_KD_BATCH_DIFF_HPP_
