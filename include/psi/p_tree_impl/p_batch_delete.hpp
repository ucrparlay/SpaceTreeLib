#ifndef PSI_P_TREE_IMPL_P_BATCH_DELETE_HPP
#define PSI_P_TREE_IMPL_P_BATCH_DELETE_HPP

#include "psi/p_tree.h"

namespace psi
{

// NOTE: default batch delete
template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
template <typename Range>
void p_tree<Point, SplitRule, SkHeight, ImbaRatio>::batch_delete(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { batch_delete_(A); });
	return;
}

// NOTE: assume all points_type are fully covered in the tree
template <typename Point, typename SplitRule, uint_fast8_t SkHeight,
	  uint_fast8_t ImbaRatio>
void p_tree<Point, SplitRule, SkHeight, ImbaRatio>::batch_delete_(slice_type A)
{
	if (!this->cpam_aug_map_.root ||
	    !cpam_aug_map_type::size(this->cpam_aug_map_.root)) {
		return;
	}
	this->cpam_aug_map_ = cpam_aug_map_type::multi_delete(
		std::move(this->cpam_aug_map_), A);
	return;
}

} // namespace psi

#endif // PSI_P_TREE_IMPL_P_BATCH_DELETE_HPP
