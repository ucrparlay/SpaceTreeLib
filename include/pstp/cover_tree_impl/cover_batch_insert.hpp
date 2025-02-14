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

  for (auto const& p : A) {
    if (!BT::WithinCircle(p, this->root_cover_circle_)) {
      ExtendCoverRangeUpwards(this->root_cover_circle_, p);
    }

    assert(BT::WithinCircle(p, this->root_cover_circle_));
    bool flag = false;
    std::tie(this->root_, flag) =
        PointInsertRecursive(this->root_, p, this->root_cover_circle_.second);
    assert(flag == true);
  }

  assert(this->root_ != nullptr);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void CoverTree<Point, SplitRule, kSkHeight,
               kImbaRatio>::ExtendCoverRangeUpwards(CoverCircle& root_cc,
                                                    Point const& p) {
  assert(!BT::WithinCircle(p, root_cc));
  Coord dis = BT::P2PDistance(p, root_cc.first);
  while (Num::Lt(static_cast<Coord>(1 << root_cc.second), dis)) {
    this->root_ = AllocInteriorNode<Interior>(
        typename Interior::CoverNodeArr(1, this->root_),
        typename Interior::ST(1, root_cc.first), typename Interior::AT());
    root_cc.second++;
  }
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
typename CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::NodeBoolean
CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::PointInsertRecursive(
    Node* T, Point const& p, [[maybe_unused]] DepthType level) {
  if (T->is_leaf) {
    auto TL = static_cast<Leaf*>(T);
    if (TL->CapacityFull()) {
      auto cover_circle = BT::GetCircle(TL->pts.cut(0, TL->size));
      if (!BT::WithinCircle(p, cover_circle)) {
        return {T, false};  // cannot insert within
      }
      auto leaf_points = TL->pts;
      auto split_node = ShrinkCoverRangeDownwards(T, cover_circle, level);
      for (auto const& point : leaf_points) {
        bool flag = false;
        std::tie(split_node, flag) =
            PointInsertRecursive(split_node, point, level);
        assert(flag == true);
      }
      return NodeBoolean(split_node, true);
    } else {
      return {BT::template InsertPoints2Leaf<Leaf>(
                  T, parlay::make_slice(Points{p})),
              true};
    }
  }

  return {T, false};
}

}  // namespace pstp

#endif  // PSTP_COVER_BATCH_INSERT_HPP_
