#ifndef PSTP_COVER_BATCH_INSERT_HPP_
#define PSTP_COVER_BATCH_INSERT_HPP_

#include <utility>

#include "../cover_tree.h"

namespace pstp {

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchInsert(
    Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);
  Slice A = parlay::make_slice(In);
  BatchInsert_(A);
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchInsert_(Slice A) {
  // Points B = Points::uninitialized(A.size());
  // this->tree_box_ = BT::GetBox(this->tree_box_, BT::GetBox(A));
  for (auto const& p : A) {
    // BUG: need to check the depth
    this->root_ = PointInsertRecursive(this->root_, p, 0);
  }
  assert(this->root_ != NULL);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
Node* CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::PointInsertRecursive(
    Node* T, Point const& p, DepthType dep) {
  if (T->is_leaf) {
    auto TL = static_cast<Leaf*>(T);
    if (BT::WithinCircle(p, Circle(TL->pts[0], (1 << dep)))) {
      Point tmp_p(p);
      std::ranges::swap(TL->pts[0], tmp_p);
      auto new_inter = AllocInteriorNode<Interior>(typename Interior::NodeArr(),
                                                   typename Interior::ST(),
                                                   typename Interior::AT());
      auto TI = static_cast<Interior*>(new_inter);
      TI->tree_nodes.push_back(T);
      TI->split.push_back(Circle(p, (1 << dep)));
      return new_inter;
    }
  }
  return T;
}

}  // namespace pstp

#endif  // PSTP_COVER_BATCH_INSERT_HPP_
