#ifndef PSI_POINTER_VIEW_P_TREE_IMPL_P_BATCH_INSERT_HPP_
#define PSI_POINTER_VIEW_P_TREE_IMPL_P_BATCH_INSERT_HPP_

#include "../p_tree.h"
#include "parlay/slice.h"

namespace psi {
namespace pointer_view {
template <typename TypeTrait>
void PTree<TypeTrait>::BatchInsert(Slice A) {
  if (!this->cpam_aug_map_.root ||
      !CpamAugMap::size(this->cpam_aug_map_.root)) {
    return Build(std::forward<Slice>(A));
  }
  this->cpam_aug_map_ =
      CpamAugMap::multi_insert(std::move(this->cpam_aug_map_), A);
  return;
}

}  // namespace pointer_view
}  // namespace psi

 
 

#endif
