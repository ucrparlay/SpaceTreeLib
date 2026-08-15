#ifndef PSI_P_TREE_IMPL_P_BATCH_DIFF_HPP_
#define PSI_P_TREE_IMPL_P_BATCH_DIFF_HPP_

#include "../p_tree.h"

namespace psi
{
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
	  uint_fast8_t kImbaRatio>
template <typename Range>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDiff(Range &&In)
{
	BT::IngestRange(In, [&](Slice A) { BatchDiff_(A); });
	return;
}

// NOTE: batch delete suitable for Points that are pratially covered in the tree
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
	  uint_fast8_t kImbaRatio>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDiff_(Slice A)
{
	if (!this->cpam_aug_map_.root) {
		return;
	}
	this->cpam_aug_map_ =
		CpamAugMap::multi_diff(std::move(this->cpam_aug_map_), A);
	return;
}

} // namespace psi

#endif // PSI_P_TREE_IMPL_P_BATCH_DIFF_HPP_
