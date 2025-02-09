#ifndef PSTP_COVER_TREE_IMPL_COVER_BUILD_TREE_HPP_
#define PSTP_COVER_TREE_IMPL_COVER_BUILD_TREE_HPP_

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include <algorithm>

#include "pstp/cover_tree.h"
#include "pstp/dependence/tree_node.h"

namespace pstp {

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::Build_(Slice A) {
  Points B = Points::uninitialized(A.size());
  this->tree_box_ = BT::GetBox(A);
  assert(this->root_ != nullptr);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::Build_(
    Slice A, Box const& box) {
  assert(BT::WithinBox(BT::GetBox(A), box));

  Points B = Points::uninitialized(A.size());
  this->tree_box_ = box;
  // this->root_ = BuildRecursive(A, B.cut(0, A.size()), this->tree_box_);
  assert(this->root_ != nullptr);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range, typename... Args>
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::Build(Range&& In,
                                                               Args&&... args) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  BatchInsert(std::forward<Range>(In), std::forward<Args>(args)...);
}
}  // namespace pstp

#endif  // PSTP_COVER_TREE_IMPL_COVER_BUILD_TREE_HPP_
