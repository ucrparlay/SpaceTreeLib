#ifndef PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_INSERT_HPP_
#define PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_INSERT_HPP_

#include "array_based/kd_tree_array.h"

namespace psi {
namespace array_based {

//==============================================================================
// BATCH INSERT OPERATIONS
//==============================================================================

template <typename TypeTrait>
void KdTreeArray<TypeTrait>::BatchInsert(Slice In) {
  // TODO: Implement
  // Strategy: Rebuild tree with new points
  // 1. Flatten existing tree
  // 2. Merge with new points
  // 3. Rebuild from scratch
  throw std::runtime_error("KdTreeArray::BatchInsert not yet implemented");
}

}  // namespace array_based
}  // namespace psi

#endif  // PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_INSERT_HPP_
