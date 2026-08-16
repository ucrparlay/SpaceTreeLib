#ifndef PSI_KD_TREE_IMPL_KD_OVERRIDE_HPP_
#define PSI_KD_TREE_IMPL_KD_OVERRIDE_HPP_

#include <type_traits>
#include <utility>

#include "psi/kd_tree.h"
#include "psi/dependence/loggers.h"
#include "psi/dependence/tree_node.h"

namespace psi
{
template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
auto kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	     ImbaRatio>::knn(Point const &q,
			     bounded_queue<Point, Range> &bq) const
{
	knn_logger logger;
	/* An empty queue has no slot to write a neighbour into. */
	if (this->root_ == nullptr || bq.max_size() == 0)
		return logger;
	if constexpr (has_box<typename interior_type::at_type>) {
		base_type::template knn_binary_box<leaf_type, interior_type>(
			this->root_, q, bq, logger);
	} else {
		base_type::template knn_binary<leaf_type, interior_type>(
			this->root_, q, bq, this->tree_box_, logger);
	}
	return logger;
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
size_t kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
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
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
auto kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
	     ImbaRatio>::range_count(box_type const &bx) const
{
	range_query_logger logger;
	if (this->root_ == nullptr || base_type::box_is_empty(bx))
		return std::make_pair(size_t{0}, logger);
	size_t size = base_type::template range_count_rectangle<leaf_type,
								interior_type>(
		this->root_, bx, this->tree_box_, logger);
	return std::make_pair(size, logger);
}

// template <typename Point, typename SplitRule, typename LeafAugType, typename
// InteriorAugType,  uint_fast8_t SkHeight,
//           uint_fast8_t ImbaRatio>
// auto kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight ,
// ImbaRatio>::range_count(
//     circle_type const& cl) {
//   return base_type::template RangeCountRadius<leaf_type,
//   interior_type>(this->root_, cl,
//                                                        this->tree_box_);
// }

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
auto kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
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
		this->tree_box_, logger);
	return std::make_pair(s, logger);
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
constexpr void kd_tree<Point, SplitRule, LeafAugType, InteriorAugType, SkHeight,
		       ImbaRatio>::delete_tree()
{
	base_type::template delete_tree_wrapper<leaf_type, interior_type>();
}

} // namespace psi

#endif // PSI_KD_TREE_IMPL_KD_OVERRIDE_HPP_
