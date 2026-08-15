#ifndef PSI_P_TREE_IMPL_P_BATCH_DELETE_HPP
#define PSI_P_TREE_IMPL_P_BATCH_DELETE_HPP

#include "../p_tree.h"

namespace psi
{

// NOTE: default batch delete
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
	  uint_fast8_t kImbaRatio>
template <typename Range>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDelete(Range &&In)
{
	BT::IngestRange(In, [&](Slice A) { BatchDelete_(A); });
	return;
}

// NOTE: assume all Points are fully covered in the tree
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
	  uint_fast8_t kImbaRatio>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDelete_(Slice A)
{
	if (!this->cpam_aug_map_.root ||
	    !CpamAugMap::size(this->cpam_aug_map_.root)) {
		return;
	}
	this->cpam_aug_map_ =
		CpamAugMap::multi_delete(std::move(this->cpam_aug_map_), A);
	return;
}

} // namespace psi

#endif // PSI_P_TREE_IMPL_P_BATCH_DELETE_HPP
