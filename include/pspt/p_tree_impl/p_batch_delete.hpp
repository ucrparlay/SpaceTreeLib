#ifndef PSPT_P_TREE_IMPL_P_BATCH_DELETE_HPP
#define PSPT_P_TREE_IMPL_P_BATCH_DELETE_HPP

#include "../p_tree.h"

namespace pspt {

// NOTE: default batch delete
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDelete(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  BatchDelete_(A);
  return;
}

// NOTE: assume all Points are fully covered in the tree
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDelete_(Slice A) {
  if (!this->cpam_aug_map_.root ||
      !CpamAugMap::size(this->cpam_aug_map_.root)) {
    return;
  }
  this->cpam_aug_map_ =
      CpamAugMap::multi_delete(std::move(this->cpam_aug_map_), A);
  return;
}

}  // namespace pspt

#endif  // PSPT_P_TREE_IMPL_P_BATCH_DELETE_HPP
