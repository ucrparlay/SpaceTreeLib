#ifndef PSI_BASE_TREE_IMPL_RANGE_QUERY_HPP_
#define PSI_BASE_TREE_IMPL_RANGE_QUERY_HPP_

#include <algorithm>
#include <utility>

#include "psi/base_tree.h"
#include "psi/dependence/concepts.h"

namespace psi
{
// NOTE: orthogonal range count in leaf
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type>
size_t
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::range_count_rectangle_leaf(
	node *T, box_type const &query_box)
{
	assert(T->is_leaf);

	leaf_type *tl = static_cast<leaf_type *>(T);
	if (tl->is_dummy) {
		return within_box(tl->pts[0], query_box) ? tl->size : 0;
	} else {
		return std::ranges::count_if(
			tl->pts.begin(), tl->pts.begin() + tl->size,
			[&](auto const &p) {
				return within_box(p, query_box);
			});
	}
}

// NOTE: recursive orthogonal range count for bianry node
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_binary_node interior_type>
size_t
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::range_count_rectangle(
	node *T, box_type const &query_box, box_type const &node_box,
	range_query_logger &logger)
{
	if (T->is_leaf) {
		// PERF: no need to check box because we already check it before
		// we recursion
		logger.vis_leaf_num++;
		return range_count_rectangle_leaf<leaf_type>(T, query_box);
	}

	logger.vis_interior_num++;
	interior_type *ti = static_cast<interior_type *>(T);

	size_t left_cnt, right_cnt;
	auto recurse = [&](node *t, box_type const &box,
			   size_t &counter) -> void {
		if (!box_intersect_box(box, query_box)) {
			logger.skip_box_num++;
			counter = 0;
		} else if (within_box(box, query_box)) {
			logger.full_box_num++;
			counter = t->size;
		} else {
			counter =
				range_count_rectangle<leaf_type, interior_type>(
					t, query_box, box, logger);
		}
	};

	if constexpr (!has_box<typename interior_type::at_type>) { // use
								   // hyperplane
		box_cut_type box_cut(node_box, ti->split, true);
		logger.generate_box_num++;
		recurse(ti->left, box_cut.get_first_box_cut(), left_cnt);
		recurse(ti->right, box_cut.get_second_box_cut(), right_cnt);
	} else if constexpr (has_box<typename interior_type::
					     at_type>) { // use bounding
							 // box
		recurse(ti->left,
			retrieve_box<leaf_type, interior_type>(ti->left),
			left_cnt);
		recurse(ti->right,
			retrieve_box<leaf_type, interior_type>(ti->right),
			right_cnt);
	} else {
		assert(0);
	}

	return left_cnt + right_cnt;
}

// NOTE: orthogonal range count for multi node
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_multi_node interior_type>
size_t
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::range_count_rectangle(
	node *T, box_type const &query_box, box_type const &node_box,
	dims_type dim, bucket_type idx, range_query_logger &logger)
{
	if (T->is_leaf) {
		logger.vis_leaf_num++;
		return range_count_rectangle_leaf<leaf_type>(T, query_box);
	}

	logger.vis_interior_num++;
	interior_type *ti = static_cast<interior_type *>(T);

	auto recurse = [&query_box, &logger](
			       node *t, box_type const &box, size_t &counter,
			       dims_type next_dim, dims_type next_idx) -> void {
		if (!box_intersect_box(box, query_box)) {
			logger.skip_box_num++;
			counter = 0;
		} else if (within_box(box, query_box)) {
			logger.full_box_num++;
			// NOTE: when reach leaf, t is the children
			counter = static_cast<interior_type *>(t)->merge_size(
				next_idx);
		} else {
			counter =
				range_count_rectangle<leaf_type, interior_type>(
					t, query_box, box, next_dim, next_idx,
					logger);
		}
	};

	size_t l, r;
	idx <<= 1;
	bool reach_leaf = idx >= interior_type::get_regions();
	box_cut_type box_cut(node_box, ti->split[dim], true);
	logger.generate_box_num++;

	// NOTE: visit left half
	assert(ti->split[dim].second == dim);

	recurse(reach_leaf ? ti->tree_nodes[idx - interior_type::get_regions()]
			   : T,
		box_cut.get_first_box_cut(), l, (dim + 1) % num_dims,
		reach_leaf ? 1 : idx);

	idx |= 1;
	recurse(reach_leaf ? ti->tree_nodes[idx - interior_type::get_regions()]
			   : T,
		box_cut.get_second_box_cut(), r, (dim + 1) % num_dims,
		reach_leaf ? 1 : idx);

	return l + r;
}

// TODO: as range_count_rectangle
// template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
//           uint_fast8_t ImbaRatio>
// template <typename leaf_type, is_binary_node interior_type>
// size_t base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::RangeCountRadius(
//     node* T, circle_type const& cl, box_type const& node_box) {
//   if (!circle_intersect_box(cl, node_box)) return 0;
//   if (within_circle(node_box, cl)) return T->size;
//
//   if (T->is_leaf) {
//     size_t cnt = 0;
//     leaf_type* tl = static_cast<leaf_type*>(T);
//     if (tl->is_dummy) {
//       if (within_circle(tl->pts[0], cl)) {
//         cnt += tl->size;
//       }
//     } else {
//       std::ranges::for_each(tl->pts.begin(), tl->pts.begin() + tl->size,
//                             [&](auto const& p) {
//                               if (within_circle(p, cl)) {
//                                 cnt++;
//                               }
//                             });
//     }
//     return cnt;
//   }
//
//   interior_type* ti = static_cast<interior_type*>(T);
//   box_type lbox(node_box), rbox(node_box);
//   lbox.second.pnt[ti->split.second] = ti->split.first;  //* loose
//   rbox.first.pnt[ti->split.second] = ti->split.first;
//
//   size_t l, r;
//   parlay::par_do_if(
//       ti->size >= serial_build_cutoff,
//       [&] { l = RangeCountRadius(ti->left, cl, lbox); },
//       [&] { r = RangeCountRadius(ti->right, cl, rbox); });
//
//   return l + r;
// };

// NOTE: orthogonal range query in leaf
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, typename Range>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::range_query_leaf(
	node *T, Range out, size_t &s, box_type const &query_box)
{
	assert(T->is_leaf);

	if (s >= out.size()) {
		return;
	}

	leaf_type *tl = static_cast<leaf_type *>(T);
	/*
	 * out is the caller's buffer and the answer size is not known until
	 * the query runs, so every write is bounded by it. When the buffer has
	 * room for the whole leaf -- the ordinary case -- the bulk paths below
	 * run untouched; otherwise the answer is truncated and the caller sees
	 * it in the returned count.
	 */
	if (tl->is_dummy) {
		if (within_box(tl->pts[0], query_box)) {
			size_t const n =
				std::min<size_t>(tl->size, out.size() - s);
			std::fill_n(out.begin() + s, n, tl->pts[0]);
			s += n;
		}
	} else if (s + tl->size <= out.size()) {
		auto result = std::ranges::copy_if(
			tl->pts.begin(), tl->pts.begin() + tl->size,
			out.begin() + s, [&](auto const &p) {
				return within_box(p, query_box);
			});
		s += std::ranges::distance(out.begin() + s, result.out);
	} else {
		for (size_t i = 0; i < tl->size && s < out.size(); i++) {
			if (within_box(tl->pts[i], query_box)) {
				out[s++] = tl->pts[i];
			}
		}
	}
	return;
}

// NOTE: orthogonal range query recusively
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_binary_node interior_type, typename Range>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::
	range_query_serial_recursive(node *T, Range out, size_t &s,
				     box_type const &query_box,
				     box_type const &node_box,
				     range_query_logger &logger)
{
	if (T->is_leaf) {
		logger.vis_leaf_num++;
		range_query_leaf<leaf_type>(T, out, s, query_box);
		return;
	}

	logger.vis_interior_num++;
	interior_type *ti = static_cast<interior_type *>(T);

	auto recurse = [&](node *t, box_type const &box) -> void {
		if (!box_intersect_box(box, query_box)) {
			logger.skip_box_num++;
			return;
		} else if (within_box(box, query_box) &&
			   s + t->size <= out.size()) {
			logger.full_box_num++;
			flatten_rec<leaf_type, interior_type>(
				t, out.cut(s, s + t->size));
			s += t->size;
			return;
		} else {
			range_query_serial_recursive<leaf_type, interior_type>(
				t, out, s, query_box, box, logger);
			return;
		}
	};

	if constexpr (!has_box<typename interior_type::at_type>) { // use
								   // hyperplane
		box_cut_type box_cut(node_box, ti->split, true);
		logger.generate_box_num++;
		recurse(ti->left, box_cut.get_first_box_cut());
		recurse(ti->right, box_cut.get_second_box_cut());
	} else if constexpr (has_box<typename interior_type::
					     at_type>) { // use bounding
							 // box
		recurse(ti->left,
			retrieve_box<leaf_type, interior_type>(ti->left));
		recurse(ti->right,
			retrieve_box<leaf_type, interior_type>(ti->right));
	} else {
		assert(0);
	}

	return;
}

template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename leaf_type, is_multi_node interior_type, typename Range>
void base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::
	range_query_serial_recursive(node *T, Range out, size_t &s,
				     box_type const &query_box,
				     box_type const &node_box, dims_type dim,
				     bucket_type idx,
				     range_query_logger &logger)
{
	if (T->is_leaf) {
		logger.vis_leaf_num++;
		range_query_leaf<leaf_type>(T, out, s, query_box);
		return;
	}

	logger.vis_interior_num++;
	interior_type *ti = static_cast<interior_type *>(T);

	auto recurse = [&query_box, &s, &out, &logger](
			       node *t, box_type const &box, dims_type next_dim,
			       bucket_type next_idx) -> void {
		if (!box_intersect_box(box, query_box)) {
			logger.skip_box_num++;
			return;
		} else if (within_box(box, query_box)) {
			size_t const candidate_size =
				static_cast<interior_type *>(t)->merge_size(
					next_idx);
			if (s + candidate_size > out.size()) {
				/* No room for the O(1) copy; the walk below
				 * clamps leaf by leaf instead. */
				range_query_serial_recursive<leaf_type,
							     interior_type>(
					t, out, s, query_box, box, next_dim,
					next_idx, logger);
				return;
			}
			logger.full_box_num++;
			partial_flatten<leaf_type, interior_type>(
				t, out.cut(s, s + candidate_size), next_idx);
			s += candidate_size;
			return;
		} else {
			range_query_serial_recursive<leaf_type, interior_type>(
				t, out, s, query_box, box, next_dim, next_idx,
				logger);
			return;
		}
	};

	idx <<= 1;
	bool reach_leaf = idx >= interior_type::get_regions();
	logger.generate_box_num++;
	box_cut_type box_cut(node_box, ti->split[dim], true);

	recurse(reach_leaf ? ti->tree_nodes[idx - interior_type::get_regions()]
			   : T,
		box_cut.get_first_box_cut(), (dim + 1) % num_dims,
		reach_leaf ? 1 : idx);

	// NOTE: visit right
	idx |= 1;
	recurse(reach_leaf ? ti->tree_nodes[idx - interior_type::get_regions()]
			   : T,
		box_cut.get_second_box_cut(), (dim + 1) % num_dims,
		reach_leaf ? 1 : idx);

	return;
}
/*
 * Convenience range_query: the answer size is not known in advance, so this
 * counts first and then reports into a buffer it owns.
 */
template <typename Point, typename DerivedTree, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
typename base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::points_type
base_tree<Point, DerivedTree, SkHeight, ImbaRatio>::range_query(
	box_type const &query_box) const
{
	auto const *tree = static_cast<DerivedTree const *>(this);
	size_t const n = tree->range_count(query_box).first;
	points_type out(n);
	if (n != 0) {
		tree->range_query(query_box, parlay::make_slice(out));
	}
	return out;
}

} // namespace psi

#endif // PSI_BASE_TREE_IMPL_RANGE_QUERY_HPP_
