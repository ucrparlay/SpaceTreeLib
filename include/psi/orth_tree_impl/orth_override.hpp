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
template <typename Traits>
template <typename Range>
auto orth_tree<Traits>::knn(point_type const &q,
			    bounded_queue<point_type, Range> &bq) const
{
	knn_logger logger;
	/* An empty queue has no slot to write a neighbour into. */
	if (this->root_ == nullptr || bq.max_size() == 0)
		return logger;
	if constexpr (has_box<typename interior_type::at_type>) {
		base_type::template knn_multi<leaf_type, interior_type>(
			this->root_, q, bq, logger);
	} else {
		base_type::template knn_multi_expand<leaf_type, interior_type>(
			this->root_, q, 0, 1, bq, this->tree_box_, logger);
	}
	return logger;
}

template <typename Traits>
template <typename Range>
size_t orth_tree<Traits>::flatten(Range &&out) const
{
	size_t const n = this->get_size();
	if (this->root_ == nullptr || out.size() != n) {
		return 0;
	}
	base_type::template flatten_rec<leaf_type, interior_type>(
		this->root_, parlay::make_slice(out));
	return n;
}

template <typename Traits>
auto orth_tree<Traits>::range_count(box_type const &bx) const
{
	range_query_logger logger;
	if (this->root_ == nullptr || base_type::box_is_empty(bx))
		return std::make_pair(size_t{0}, logger);
	size_t s = base_type::template range_count_rectangle<leaf_type,
							     interior_type>(
		this->root_, bx, this->tree_box_, 0, 1, logger);
	return std::make_pair(s, logger);
}

template <typename Traits>
template <typename Range>
auto orth_tree<Traits>::range_query(box_type const &query_box,
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

template <typename Traits>
constexpr void orth_tree<Traits>::delete_tree()
{
	base_type::template delete_tree_wrapper<leaf_type, interior_type>();
	this->fixed_box = false;
}

} /* namespace psi */

#endif /* PSI_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_ */
