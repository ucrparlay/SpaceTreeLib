#ifndef PSI_P_TREE_IMPL_P_BUILD_TREE_HPP_
#define PSI_P_TREE_IMPL_P_BUILD_TREE_HPP_

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include "psi/p_tree.h"
#include "parlay/utilities.h"
#include "psi/dependence/tree_node.h"

namespace psi
{
template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
void p_tree<Point, SplitRule, SkHeight, ImbaRatio>::build(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { build_(A); });
}

template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void p_tree<Point, SplitRule, SkHeight, ImbaRatio>::build_(slice_type A)
{

	this->cpam_aug_map_ = std::move(cpam_aug_map_type(A));

	return;
}

} /* namespace psi */

#endif /* PSI_P_TREE_IMPL_P_BUILD_TREE_HPP_ */
