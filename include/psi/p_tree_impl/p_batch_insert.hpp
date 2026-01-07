#ifndef PSI_P_TREE_IMPL_P_BATCH_INSERT_HPP_
#define PSI_P_TREE_IMPL_P_BATCH_INSERT_HPP_

#include "../p_tree.h"
#include "parlay/slice.h"

namespace psi {
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchInsert(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  auto aux = Points::uninitialized(parlay::size(In));
  parlay::copy(In, parlay::make_slice(aux));
  Slice A = parlay::make_slice(aux);
  BatchInsert_(A);
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchInsert_(Slice A) {
  if (!this->cpam_aug_map_.root ||
      !CpamAugMap::size(this->cpam_aug_map_.root)) {
    return Build(std::forward<Slice>(A));
  }
  this->cpam_aug_map_ =
      CpamAugMap::multi_insert(std::move(this->cpam_aug_map_), A);
  return;
}

}  // namespace psi

#endif  // PSI_P_TREE_IMPL_P_BATCH_INSERT_HPP_
