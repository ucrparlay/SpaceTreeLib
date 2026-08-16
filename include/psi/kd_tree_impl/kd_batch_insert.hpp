#ifndef PSI_KD_TREE_IMPL_KD_BATCH_INSERT_HPP_
#define PSI_KD_TREE_IMPL_KD_BATCH_INSERT_HPP_

#include "../kd_tree.h"
#include "parlay/slice.h"

namespace psi
{
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
void kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	     ImbaRatio>::batch_insert(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { batch_insert_(A); });
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	     ImbaRatio>::batch_insert_(slice_type A)
{
	if (this->root_ == nullptr) { // TODO: may check using explicity tag
		return build_(A);
	}

	points_type B = points_type::uninitialized(A.size());
	node *T = this->root_;
	dims_type d =
		T->is_leaf ? 0 : static_cast<interior_type *>(T)->split.second;
	this->root_ = batch_insert_recursive(T, A, B.cut(0, A.size()), d);
	assert(this->root_ != nullptr);
	if constexpr (has_box<typename interior_type::at_type>) {
		this->tree_box_ = base_type::template retrieve_box<
			leaf_type, interior_type>(this->root_);
	} else {
		this->tree_box_ = base_type::get_box(this->tree_box_,
						     base_type::get_box(A));
	}
	return;
}

// NOTE: return the updated node
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
node *kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	      ImbaRatio>::batch_insert_recursive(node *T, slice_type in,
						 slice_type out, dims_type d)
{
	size_t n = in.size();

	if (n == 0)
		return T;

	auto prepare_func = [](node *, points_type &wp, points_type &) {
		return base_type::get_box(parlay::make_slice(wp));
	};

	if (T->is_leaf) {
		leaf_type *tl = static_cast<leaf_type *>(T);
		if ((!tl->is_dummy &&
		     n + T->size <= base_type::leaf_capacity) ||
		    (tl->is_dummy && parlay::all_of(in, [&](Point const &p) {
			     return p == tl->pts[0];
		     }))) {
			return base_type::template insert_points2_leaf<
				leaf_type>(T, in);
			// auto o = base_type::template
			// insert_points2_leaf<leaf_type>(T, in);
			// assert(base_type::same_box(base_type::template
			// get_box<leaf_type, interior_type>(o),
			//                    base_type::template
			//                    retrieve_box<leaf_type,
			//                    interior_type>(o)));
			// return o;
		} else { // PERF: if a nomarl leaf tl cannot handle more
			 // duplicates, leave them here
			return base_type::template rebuild_with_insert<
				leaf_type, interior_type>(T, prepare_func, in,
							  d);
		}
	}

	if (n <= base_type::serial_build_cutoff) {
		interior_type *ti = static_cast<interior_type *>(T);
		std::ranges::subrange _2ndGroup =
			std::ranges::partition(in, [&](Point const &p) {
				return num_type::lt(p.pnt[ti->split.second],
						    ti->split.first);
			});
		size_t split_pos =
			static_cast<points_iter_type>(_2ndGroup.begin()) -
			static_cast<points_iter_type>(in.begin());

		// NOTE: rebuild
		if (split_rule_.allow_rebuild() &&
		    base_type::imbalance_node(ti->left->size + split_pos,
					      ti->size + n)) {
			return base_type::template rebuild_with_insert<
				leaf_type, interior_type>(T, prepare_func, in,
							  d);
		}

		// NOTE: continue
		node *L, *R;
		d = split_rule_.next_dimension(d);
		L = batch_insert_recursive(ti->left, in.cut(0, split_pos),
					   out.cut(0, split_pos), d);
		R = batch_insert_recursive(ti->right, in.cut(split_pos, n),
					   out.cut(split_pos, n), d);
		base_type::template update_interior<interior_type>(T, L, R);
		assert(T->size == L->size + R->size && ti->split.second >= 0 &&
		       ti->is_leaf == false);
		return T;
	}

	// NOTE: assign each node a tag
	inner_tree IT;
	assert(IT.rev_tag.size() == base_type::bucket_num);
	IT.assign_node_tag(T, 1);
	assert(IT.tags_num > 0 && IT.tags_num <= base_type::bucket_num);

	base_type::template sieve_points<interior_type>(in, out, n, IT.tags,
							IT.sums, IT.tags_num);

	IT.tag_imbalance_node([&]([[maybe_unused]] bucket_type idx) -> bool {
		auto const ti =
			static_cast<interior_type *>(IT.tags[idx].first);
		return split_rule_.allow_rebuild() &&
		       base_type::imbalance_node(ti->left->size +
							 IT.sums_tree[idx << 1],
						 ti->size + IT.sums_tree[idx]);
	});
	assert(IT.tags_num > 0 && IT.tags_num <= base_type::bucket_num);
	auto tree_nodes = parlay::sequence<node *>::uninitialized(IT.tags_num);

	parlay::parallel_for(
		0, IT.tags_num,
		[&](decltype(IT.tags_num) i) {
			size_t s = 0;
			for (decltype(IT.tags_num) j = 0; j < i; j++) {
				s += IT.sums_tree[IT.rev_tag[j]];
			}

			dims_type next_dim = d,
				  depth = IT.get_depth_by_index(IT.rev_tag[i]);
			for (bucket_type i = 0; i < depth; i++) {
				next_dim = split_rule_.next_dimension(next_dim);
			}

			if (IT.tags[IT.rev_tag[i]].second ==
			    base_type::bucket_num + 1) {
				// NOTE: continue sieve
				tree_nodes[i] = batch_insert_recursive(
					IT.tags[IT.rev_tag[i]].first,
					out.cut(s,
						s + IT.sums_tree[IT.rev_tag
									 [i]]),
					in.cut(s,
					       s + IT.sums_tree[IT.rev_tag[i]]),
					next_dim);
			} else { // NOTE: launch rebuild subtree
				assert(IT.tags[IT.rev_tag[i]].second ==
				       base_type::bucket_num + 2);
				assert(IT.tags[IT.rev_tag[i]].first->size +
					       IT.sums_tree[IT.rev_tag[i]] >=
				       0);

				tree_nodes[i] = base_type::template rebuild_with_insert<
					leaf_type, interior_type>(
					IT.tags[IT.rev_tag[i]].first,
					prepare_func,
					out.cut(s,
						s + IT.sums_tree[IT.rev_tag
									 [i]]),
					next_dim);
			}
		},
		1);

	return IT.template update_inner_tree<inner_tree::kUpdatePointer>(
		tree_nodes);
}

} // namespace psi

#endif // PSI_KD_TREE_IMPL_KD_BATCH_INSERT_HPP_
