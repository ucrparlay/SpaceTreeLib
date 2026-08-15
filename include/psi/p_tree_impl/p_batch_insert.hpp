#ifndef PSI_P_TREE_IMPL_P_BATCH_INSERT_HPP_
#define PSI_P_TREE_IMPL_P_BATCH_INSERT_HPP_

#include "../p_tree.h"
#include "parlay/slice.h"

namespace psi
{
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
	  uint_fast8_t kImbaRatio>
template <typename Range>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchInsert(Range &&In)
{
	BT::IngestRange(In, [&](Slice A) { BatchInsert_(A); });
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
	  uint_fast8_t kImbaRatio>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchInsert_(Slice A)
{
	if (!this->cpam_aug_map_.root ||
	    !CpamAugMap::size(this->cpam_aug_map_.root)) {
		/* Build_ , not Build: A is already tree owned storage. */
		return Build_(A);
	}
	this->cpam_aug_map_ =
		CpamAugMap::multi_insert(std::move(this->cpam_aug_map_), A);
	return;
}

} // namespace psi

#endif // PSI_P_TREE_IMPL_P_BATCH_INSERT_HPP_
