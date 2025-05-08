#ifndef PSPT_P_TREE_IMPL_P_BATCH_INSERT_HPP_
#define PSPT_P_TREE_IMPL_P_BATCH_INSERT_HPP_

#include "../p_tree.h"
#include "parlay/slice.h"

namespace pspt {
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchInsert(Slice A) {
  if (!this->cpam_aug_map_.root ||
      !CpamAugMap::size(this->cpam_aug_map_.root)) {
    return Build(std::forward<Slice>(A));
  }

  parlay::internal::timer t("BatchInsert");
  // auto insert_map = CpamAugMap(A);
  // t.next("build insert map");
  // this->cpam_aug_map_ =
  //     CpamAugMap::map_union(std::move(this->cpam_aug_map_),
  //                           std::move(insert_map));
  this->cpam_aug_map_ =
      CpamAugMap::multi_insert(std::move(this->cpam_aug_map_), A);
  t.next("insert map union");
  puts(">>>>>>>>>>>>>>>>>>>");
  return;
}

}  // namespace pspt

#endif  // PSPT_P_TREE_IMPL_P_BATCH_INSERT_HPP_
