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
  }
  assert(this->root_ != NULL);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
typename CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::NodeBoolean
CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::PointInsertRecursive(
    Node* T, Point const& center, Point const& p, DepthType dep) {
  if (T->is_leaf) {
    auto TL = static_cast<Leaf*>(T);
    if (TL->size == 0) {  // TODO: check usage again
      TL->pts[0] = p;
      TL->size = 1;
      return NodeBoolean(T, true);
    }
    if (BT::WithinCircle(p, Circle(TL->pts[0], (1 << dep)))) {
      Point tmp_p(p);
      std::ranges::swap(TL->pts[0], tmp_p);
      auto new_inter = AllocInteriorNode<Interior>(typename Interior::NodeArr(),
                                                   typename Interior::ST(),
                                                   typename Interior::AT());
      auto TI = static_cast<Interior*>(new_inter);
      TI->tree_nodes.push_back(T);
      TI->split.push_back(Circle(p, (1 << dep)));
      return NodeBoolean(new_inter, true);
    } else {  // the current node connot cover the point, needs separation
      assert(Num::Gt(BT::P2PDistance(TL->pts[0], p), (1 << dep)));
      return NodeBoolean(nullptr, false);
    }
  }

  auto TI = static_cast<Interior*>(T);
  BallsType cover_node_cnt = 0;
  parlay::sequence<BallsType> cover_nodes(TI->split.size());
  for (BallsType i = 0; i < TI->split.size(); i++) {
    if (BT::WithinCircle(p, TI->split[i])) {  // BUG: check the circle
      cover_nodes[cover_node_cnt++] = i;
    }
  }
  if (cover_node_cnt == 0) {
    return NodeBoolean(nullptr, false);
  }
  // try insert to next level
  for (BallsType i = 0; i < cover_node_cnt; i++) {
    auto res =
        PointInsertRecursive(TI->tree_nodes[cover_nodes[i]],
                             TI->split[cover_nodes[i]].first, p, dep - 1);
    if (res.second) {
      TI->tree_nodes[cover_nodes[i]] = res.first;
      return NodeBoolean(T, true);
    }
  }
  // try insert to current level
  if (BT::WithinCircle(p, Circle(center, (1 << dep)))) {
    // BUG: it could be a leaf node
    auto target = static_cast<Interior*>(TI->tree_nodes[cover_nodes[0]]);
    // target->tree_nodes.emplace_back(AllocEmptyLeafNode<Slice, Leaf>());
    target->tree_nodes.emplace_back(
        AllocNormalLeafNode<Slice, Leaf>(parlay::make_slice(Points(1, p))));
    target->split.push_back(Circle(p, (1 << (dep - 1))));
    return NodeBoolean(T, true);
  } else {  // fail, needs to insert previous level
    return NodeBoolean(T, false);
  }
}

}  // namespace pstp

#endif  // PSTP_COVER_BATCH_INSERT_HPP_
