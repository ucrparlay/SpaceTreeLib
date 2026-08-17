#ifndef PSI_P_TREE_IMPL_P_BATCH_DELETE_HPP
#define PSI_P_TREE_IMPL_P_BATCH_DELETE_HPP

#include "psi/p_tree.h"

namespace psi
{

/* default batch delete */
template <typename Traits>
template <typename Range>
void p_tree<Traits>::batch_delete(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { batch_delete_(A); });
	return;
}

/* assume all points_type are fully covered in the tree */
template <typename Traits>
void p_tree<Traits>::batch_delete_(slice_type A)
{
	if (!this->cpam_aug_map_.root ||
	    !cpam_aug_map_type::size(this->cpam_aug_map_.root)) {
		return;
	}
	this->cpam_aug_map_ = cpam_aug_map_type::multi_delete(
		std::move(this->cpam_aug_map_), A);
	return;
}

} /* namespace psi */

#endif /* PSI_P_TREE_IMPL_P_BATCH_DELETE_HPP */
