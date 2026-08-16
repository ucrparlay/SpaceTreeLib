#ifndef PSI_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_
#define PSI_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_

#include <tuple>
#include <type_traits>
#include <utility>

#include "psi/orth_tree.h"
#include "psi/dependence/concepts.h"
#include "psi/dependence/tree_node.h"

namespace psi
{
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
auto orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::knn(Point const &q,
			       bounded_queue<Point, Range> &bq) const
{
	knn_logger logger;
	/* An empty queue has no slot to write a neighbour into. */
	if (this->root_ == nullptr || bq.max_size() == 0)
		return logger;
	// base_type::template knn_multi<leaf_type, interior_type>(this->root_,
	// q, DIM, bq, this->tree_box_,
	//                                       vis_node_num, generate_box_num,
	//                                       check_box_num);
	// base_type::template knn_binary<leaf_type, kd_interior_node>(T, q,
	// DIM, bq, this->tree_box_,
	//                                              logger);
	if constexpr (has_box<typename interior_type::at_type>) {
		base_type::template knn_multi<leaf_type, interior_type>(
			this->root_, q, bq, logger);
		// base_type::template knn_multi_expand_box<leaf_type,
		// interior_type>(T, q, 0, 1, bq, logger);
	} else {
		base_type::template knn_multi_expand<leaf_type, interior_type>(
			this->root_, q, 0, 1, bq, this->tree_box_, logger);
	}
	return logger;
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
size_t orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
		 ImbaRatio>::flatten(Range &&out) const
{
	size_t const n = this->get_size();
	if (this->root_ == nullptr || out.size() != n) {
		return 0;
	}
	base_type::template flatten_rec<leaf_type, interior_type>(
		this->root_, parlay::make_slice(out));
	return n;
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
auto orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::range_count(box_type const &bx) const
{
	range_query_logger logger;
	if (this->root_ == nullptr || base_type::box_is_empty(bx))
		return std::make_pair(size_t{0}, logger);
	size_t s = base_type::template range_count_rectangle<leaf_type,
							     interior_type>(
		this->root_, bx, this->tree_box_, 0, 1, logger);
	return std::make_pair(s, logger);
}

// template <typename Point, typename SplitRule, typename LeafAugType, typename
// InteriorAugType, uint_fast8_t md,
//           uint_fast8_t SkHeight, uint_fast8_t ImbaRatio>
// auto orth_tree<Point, SplitRule, md, SkHeight, ImbaRatio>::range_count(
//     circle_type const& cl) {
//   return base_type::template RangeCountRadius<leaf_type,
//   interior_type>(this->root_, cl,
//                                                        this->tree_box_);
// }

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
auto orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::range_query(box_type const &query_box,
				       Range &&out) const
{
	range_query_logger logger;
	size_t s = 0;
	if (this->root_ == nullptr || base_type::box_is_empty(query_box))
		return std::make_pair(s, logger);
	base_type::template range_query_serial_recursive<leaf_type,
							 interior_type>(
		this->root_, parlay::make_slice(out), s, query_box,
		this->tree_box_, 0, 1, logger);
	return std::make_pair(s, logger);
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
constexpr void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md,
			 SkHeight, ImbaRatio>::delete_tree()
{
	base_type::template delete_tree_wrapper<leaf_type, interior_type>();
	this->fixed_box = false;
}

} // namespace psi

#endif // PSI_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_
