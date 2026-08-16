#ifndef PSI_P_TREE_IMPL_P_BATCH_DIFF_HPP_
#define PSI_P_TREE_IMPL_P_BATCH_DIFF_HPP_

#include "psi/p_tree.h"

namespace psi
{
template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
void p_tree<Point, SplitRule, SkHeight, ImbaRatio>::batch_diff(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { batch_diff_(A); });
	return;
}

/*
 * batch delete suitable for points_type that are pratially covered in the
 * tree
 */
template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void p_tree<Point, SplitRule, SkHeight, ImbaRatio>::batch_diff_(slice_type A)
{
	if (!this->cpam_aug_map_.root) {
		return;
	}
	this->cpam_aug_map_ = cpam_aug_map_type::multi_diff(
		std::move(this->cpam_aug_map_), A);
	return;
}

} /* namespace psi */

#endif /* PSI_P_TREE_IMPL_P_BATCH_DIFF_HPP_ */
