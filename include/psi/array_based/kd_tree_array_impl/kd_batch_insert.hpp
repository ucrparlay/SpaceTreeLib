#ifndef PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_INSERT_HPP_
#define PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_INSERT_HPP_

#include "array_based/kd_tree_array.h"

namespace psi {
namespace array_based {

#define KDTREEARRAY_TEMPLATE                                          \
  template <typename Point, typename SplitRule, typename NodeAugType, \
            uint_fast8_t kInnerNodeLevels, uint_fast8_t kSkHeight,    \
            uint_fast8_t kImbaRatio>

#define KDTREEARRAY_CLASS                                                 \
  KdTreeArray<Point, SplitRule, NodeAugType, kInnerNodeLevels, kSkHeight, \
              kImbaRatio>

//==============================================================================
// BATCH INSERT OPERATIONS
//==============================================================================

KDTREEARRAY_TEMPLATE
void KDTREEARRAY_CLASS::BatchInsert(Slice In) {
  // TODO: Implement
  // Strategy: Rebuild tree with new points
  // 1. Flatten existing tree
  // 2. Merge with new points
  // 3. Rebuild from scratch
  throw std::runtime_error("KdTreeArray::BatchInsert not yet implemented");
}

#undef KDTREEARRAY_TEMPLATE
#undef KDTREEARRAY_CLASS

}  // namespace array_based
}  // namespace psi

#endif  // PSI_ARRAY_BASED_KD_TREE_ARRAY_IMPL_KD_BATCH_INSERT_HPP_
