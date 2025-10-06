#ifndef PSI_P_TREE_IMPL_P_BATCH_INSERT_HPP_
#define PSI_P_TREE_IMPL_P_BATCH_INSERT_HPP_

#define PTREE_TEMPLATE template <typename Point, typename SplitRule, \
    uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
#define PTREE_CLASS PTree<Point, SplitRule, kSkHeight, kImbaRatio>

#include "../p_tree.h"
#include "parlay/slice.h"

namespace psi {
PTREE_TEMPLATE
void PTREE_CLASS::BatchInsert(Slice A) {
  if (!this->cpam_aug_map_.root ||
      !CpamAugMap::size(this->cpam_aug_map_.root)) {
    return Build(std::forward<Slice>(A));
  }
  this->cpam_aug_map_ =
      CpamAugMap::multi_insert(std::move(this->cpam_aug_map_), A);
  return;
}

}  // namespace psi

#undef PTREE_TEMPLATE
#undef PTREE_CLASS

#endif
