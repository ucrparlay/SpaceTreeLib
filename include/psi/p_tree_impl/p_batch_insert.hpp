#ifndef PSI_P_TREE_IMPL_P_BATCH_INSERT_HPP_
#define PSI_P_TREE_IMPL_P_BATCH_INSERT_HPP_

#include "psi/p_tree.h"
#include "parlay/slice.h"

namespace psi
{
template <typename Traits>
template <typename Range>
void p_tree<Traits>::batch_insert(Range &&in)
{
	base_type::ingest_range(in, [&](slice_type A) { batch_insert_(A); });
}

template <typename Traits>
void p_tree<Traits>::batch_insert_(slice_type A)
{
	if (!this->cpam_aug_map_.root ||
	    !cpam_aug_map_type::size(this->cpam_aug_map_.root)) {
		/* build_ , not build: A is already tree owned storage. */
		return build_(A);
	}
	this->cpam_aug_map_ = cpam_aug_map_type::multi_insert(
		std::move(this->cpam_aug_map_), A);
	return;
}

} /* namespace psi */

#endif /* PSI_P_TREE_IMPL_P_BATCH_INSERT_HPP_ */
