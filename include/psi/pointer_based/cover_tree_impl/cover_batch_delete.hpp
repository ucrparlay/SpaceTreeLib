#ifndef PSI_POINTER_BASED_COVER_TREE_IMPL_COVER_BATCH_DELETE_HPP_
#define PSI_POINTER_BASED_COVER_TREE_IMPL_COVER_BATCH_DELETE_HPP_

#include "../cover_tree.h"

namespace psi {
namespace pointer_based {

// NOTE: default batch delete
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDelete(
    Range&& In) {
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
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDelete_(Slice A) {
  Points B = Points::uninitialized(A.size());
  this->root_ = BatchDeleteRecursive(this->root_, A, parlay::make_slice(B),
                                     this->tree_box_, 1);
  return;
}

}  // namespace pointer_based
}  // namespace psi

#endif  // PSI_POINTER_BASED_COVER_TREE_IMPL_COVER_BATCH_DELETE_HPP_