#ifndef PSPT_P_TREE_IMPL_P_BATCH_INSERT_HPP_
#define PSPT_P_TREE_IMPL_P_BATCH_INSERT_HPP_

#include "../p_tree.h"
#include "parlay/slice.h"

namespace pspt {
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchInsert(Slice A) {
  return;
}

}  // namespace pspt

#endif  // PSPT_P_TREE_IMPL_P_BATCH_INSERT_HPP_
