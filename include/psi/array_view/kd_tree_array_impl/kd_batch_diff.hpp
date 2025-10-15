#ifndef PSI_ARRAY_VIEW_KD_TREE_ARRAY_IMPL_KD_BATCH_DIFF_HPP_
#define PSI_ARRAY_VIEW_KD_TREE_ARRAY_IMPL_KD_BATCH_DIFF_HPP_

#include "array_view/kd_tree_array.h"
namespace psi {
namespace array_view {

//==============================================================================
// BATCH DIFF OPERATIONS
//==============================================================================

template<typename TypeTrait>
template <typename Range>
void KdTreeArray<TypeTrait>::BatchDiff(Range&& In) {
  // TODO: Implement
  // Similar to BatchDelete but handles points not in tree
  throw std::runtime_error("KdTreeArray::BatchDiff not yet implemented");
}

}  // namespace array_view
}  // namespace psi

#endif  // PSI_ARRAY_VIEW_KD_TREE_ARRAY_IMPL_KD_BATCH_DIFF_HPP_
