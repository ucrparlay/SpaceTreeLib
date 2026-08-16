#ifndef PSI_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_
#define PSI_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_

#include <tuple>
#include <type_traits>
#include <utility>

#include "../orth_tree.h"
#include "dependence/concepts.h"
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
	if (this->root_ == nullptr)
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
void orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::flatten(Range &&out) const
{
	if (this->root_ == nullptr) {
		assert(out.size() == 0);
		return;
	}
	base_type::template flatten_rec<leaf_type, interior_type>(
		this->root_, parlay::make_slice(out));
}

template <typename Point, typename SplitRule, typename LeafAugType,
	  typename InteriorAugType, uint_fast8_t md, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
auto orth_tree<Point, SplitRule, LeafAugType, InteriorAugType, md, SkHeight,
	       ImbaRatio>::range_count(box_type const &bx) const
{
	range_query_logger logger;
	if (this->root_ == nullptr)
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
	if (this->root_ == nullptr)
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
