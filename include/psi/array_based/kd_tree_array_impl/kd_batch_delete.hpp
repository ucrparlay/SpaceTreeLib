#ifndef PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_DELETE_HPP_
#define PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_DELETE_HPP_

#include "array_based/kd_tree_array.h"
namespace psi {
namespace array_based {

//==============================================================================
// BATCH DELETE OPERATIONS
//==============================================================================

template<typename TypeTrait>
template <typename Range>
void KdTreeArray<TypeTrait>::BatchDelete(Range&& In) {
  // TODO: Implement
  // Strategy: Rebuild tree without deleted points
  throw std::runtime_error("KdTreeArray::BatchDelete not yet implemented");
}

}  // namespace array_based
}  // namespace psi

#endif  // PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_DELETE_HPP_
