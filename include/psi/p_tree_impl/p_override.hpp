#ifndef PSI_P_TREE_IMPL_P_OVERRIDE_HPP_
#define PSI_P_TREE_IMPL_P_OVERRIDE_HPP_

#include <type_traits>
#include <utility>

#include "psi/p_tree.h"
#include "psi/dependence/loggers.h"
#include "psi/dependence/tree_node.h"

namespace psi
{
template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
auto p_tree<Point, SplitRule, SkHeight, ImbaRatio>::knn(
	Point const &q, bounded_queue<Point, Range> &bq) const
{
	knn_logger logger;
	/* An empty queue has no slot to write a neighbour into. */
	if (bq.max_size() == 0)
		return logger;
	cpam_aug_map_type::template knn<base_type>(this->cpam_aug_map_, q, bq,
						   logger);
	return logger;
}

template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
size_t p_tree<Point, SplitRule, SkHeight, ImbaRatio>::flatten(Range &&out) const
{
	size_t const n = this->cpam_aug_map_.size();
	if (out.size() != n) {
		return 0;
	}
	cpam_aug_map_type::entries(this->cpam_aug_map_,
				   parlay::make_slice(out).begin());
	return n;
}

template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
auto p_tree<Point, SplitRule, SkHeight, ImbaRatio>::range_count(
	box_type const &query_box) const
{
	range_query_logger logger;
	if (base_type::box_is_empty(query_box)) {
		return std::make_pair(size_t{0}, logger);
	}

	auto size = cpam_aug_map_type::template range_count_filter2<base_type>(
		this->cpam_aug_map_, query_box, logger);

	return std::make_pair(size, logger);
}

template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
auto p_tree<Point, SplitRule, SkHeight, ImbaRatio>::range_query(
	box_type const &query_box, Range &&out) const
{
	range_query_logger logger;
	size_t cnt = 0;
	if (base_type::box_is_empty(query_box)) {
		return std::make_pair(cnt, logger);
	}
	cpam_aug_map_type::template range_report_filter2<base_type>(
		this->cpam_aug_map_, query_box, cnt, parlay::make_slice(out),
		logger);
	return std::make_pair(cnt, logger);
}

template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
constexpr void p_tree<Point, SplitRule, SkHeight, ImbaRatio>::delete_tree()
{
	cpam_aug_map_.clear();
	cpam_aug_map_.root = nullptr;
	return;
}

} // namespace psi

#endif // PSI_P_TREE_IMPL_P_OVERRIDE_HPP_
