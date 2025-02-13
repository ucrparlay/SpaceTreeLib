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
  if (this->root_ == nullptr) {
    this->root_ = AllocEmptyLeafNode<Slice, Leaf>();
  }
  // TODO: we don't need the radius of the point in each level
  for (auto const& p : A) {
    // BUG: need to check the depth
    // BUG: need to compute the center
    // NOTE: the center is the center for all circles in one level
    bool flag = false;
    std::tie(this->root_, flag) =
        PointInsertRecursive(this->root_, this->center_, p, 0);
    if (flag == false) {
      assert(!BT::WithinCircle(p, this->root_cover_circle_));
      BuildUpwards(this->root_, p);
      std::tie(this->root_, flag) =
          PointInsertRecursive(this->root_, this->center_, p, 0);
      assert(flag == true);
    }
  }
  assert(this->root_ != NULL);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::BuildUpwards(
    CoverCircle& root_cc, Point const& p) {
  assert(!BT::WithinCircle(p, root_cc));
  Coord dis = BT::P2PDistance(p, root_cc.first);
  while (Num::Lt(static_cast<Coord>(1 << root_cc.second), dis)) {
    this->root_ = AllocInteriorNode<Interior>(
        CoverNodeArr(this->root_, 1), Splitter(root_cc.first, 1), AugType());
    root_cc.second++;
  }
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
typename CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::NodeBoolean
CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::PointInsertRecursive(
    Node* T, Point const& center, Point const& p, DepthType dep) {
  return {T, false};
}

}  // namespace pstp

#endif  // PSTP_COVER_BATCH_INSERT_HPP_
