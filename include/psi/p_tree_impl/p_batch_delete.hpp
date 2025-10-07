#ifndef PSI_P_TREE_IMPL_P_BATCH_DELETE_HPP
#define PSI_P_TREE_IMPL_P_BATCH_DELETE_HPP

#define PTREE_TEMPLATE                                                  \
  template <typename Point, typename SplitRule, uint_fast8_t kSkHeight, \
            uint_fast8_t kImbaRatio>
#define PTREE_CLASS PTree<Point, SplitRule, kSkHeight, kImbaRatio>

#include "../p_tree.h"

namespace psi {

// NOTE: default batch delete
PTREE_TEMPLATE
template <typename Range>
void PTREE_CLASS::BatchDelete(Range&& In) {
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
PTREE_TEMPLATE
void PTREE_CLASS::BatchDelete_(Slice A) {
  if (!this->cpam_aug_map_.root ||
      !CpamAugMap::size(this->cpam_aug_map_.root)) {
    return;
  }
  this->cpam_aug_map_ =
      CpamAugMap::multi_delete(std::move(this->cpam_aug_map_), A);
  return;
}

}  // namespace psi

#undef PTREE_TEMPLATE
#undef PTREE_CLASS

#endif
