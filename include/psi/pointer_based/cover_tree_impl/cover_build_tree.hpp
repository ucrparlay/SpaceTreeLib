#ifndef PSI_COVER_TREE_IMPL_COVER_BUILD_TREE_HPP_
#define PSI_COVER_TREE_IMPL_COVER_BUILD_TREE_HPP_

#define COVERTREE_TEMPLATE                                              \
  template <typename Point, typename SplitRule, uint_fast8_t kSkHeight, \
            uint_fast8_t kImbaRatio>
#define COVERTREE_CLASS CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include <algorithm>

#include "psi/cover_tree.h"
#include "../../dependence/tree_node.h"

namespace psi {

COVERTREE_TEMPLATE
void COVERTREE_CLASS::Build_(Slice A) {
  Points B = Points::uninitialized(A.size());
  this->tree_box_ = BT::GetBox(A);
  assert(this->root_ != nullptr);
  return;
}

COVERTREE_TEMPLATE
void COVERTREE_CLASS::Build_(Slice A, Box const& box) {
  assert(BT::WithinBox(BT::GetBox(A), box));

  Points B = Points::uninitialized(A.size());
  this->tree_box_ = box;
  // this->root_ = BuildRecursive(A, B.cut(0, A.size()), this->tree_box_);
  assert(this->root_ != nullptr);
  return;
}

COVERTREE_TEMPLATE
template <typename Range, typename... Args>
void COVERTREE_CLASS::Build(Range&& In, Args&&... args) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  BatchInsert(std::forward<Range>(In), std::forward<Args>(args)...);
}
}  // namespace psi

#undef COVERTREE_TEMPLATE
#undef COVERTREE_CLASS

#endif
