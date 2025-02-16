#ifndef PSPT_COVER_BATCH_INSERT_HPP_
#define PSPT_COVER_BATCH_INSERT_HPP_

#include <utility>

#include "../cover_tree.h"
#include "dependence/tree_node.h"

namespace pspt {

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
        PointInsertRecursive(this->root_, p, this->root_cover_circle_);
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
  Coord dis = BT::P2PDistance(p, root_cc.center);
  while (Num::Lt(static_cast<Coord>(1 << root_cc.level), dis)) {
    root_cc.level++;
    this->root_ = AllocInteriorNode<Interior>(
        typename Interior::CoverNodeArr(1, this->root_),
        typename Interior::ST(1, root_cc.center),
        typename Interior::AT{
            .cover_circle = root_cc,
            .parallel_flag = decltype(Interior::AT::parallel_flag)()});
  }
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
Node* CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::
    ShrinkCoverRangeDownwards(Node* T, CoverCircle const& level_cover_circle) {
  assert(T->is_leaf);
  auto TL = static_cast<Leaf*>(T);
  Coord max_dis = std::numeric_limits<Coord>::lowest();
  // TODO: can replace it with sorting to use the fix thin leave wrap
  for (size_t i = 0; i < TL->size; ++i) {
    max_dis = std::max(max_dis,
                       BT::P2PDistance(TL->pts[i], level_cover_circle.center));
  }
  assert(Num::Leq(max_dis, static_cast<Coord>(1 << level_cover_circle.level)));

  auto cover_circle = level_cover_circle;
  while (Num::Leq(max_dis, static_cast<Coord>(1 << cover_circle.level))) {
    cover_circle.level--;
  }
  Node* node = AllocEmptyLeafNode<Slice, Leaf>();
  while (cover_circle.level < level_cover_circle.level) {
    node = AllocInteriorNode<Interior>(
        typename Interior::CoverNodeArr(1, node),
        typename Interior::ST(1, level_cover_circle.center),
        typename Interior::AT{
            .cover_circle = cover_circle,
            .parallel_flag = decltype(Interior::AT::parallel_flag)()});
    cover_circle.level++;
  }
  assert(cover_circle.level == level_cover_circle.level);
  return node;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
typename CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::NodeBoolean
CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::PointInsertRecursive(
    Node* T, Point const& p, CoverCircle const& level_cover_circle) {
  if (T->is_leaf) {
    auto TL = static_cast<Leaf*>(T);
    if (TL->CapacityFull()) {  // leaf is full
      bool flag = false;
      auto split_node = ShrinkCoverRangeDownwards(T, level_cover_circle);
      for (auto const& pt : TL->pts) {
        std::tie(split_node, flag) =
            PointInsertRecursive(split_node, pt, level_cover_circle);
        assert(flag == true);
      }
      FreeNode<Leaf>(T);
      std::tie(split_node, flag) =
          PointInsertRecursive(split_node, p, level_cover_circle);
      return NodeBoolean(split_node, flag);
    } else {
      return {BT::template InsertPoints2Leaf<Leaf>(
                  T, parlay::make_slice(Points{p})),
              true};
    }
  }

  auto TI = static_cast<Interior*>(T);
  parlay::sequence<size_t> near_node_idx_seq(TI->tree_nodes.size());
  size_t near_node_cnt = 0;
  for (size_t i = 0; i < TI->tree_nodes.size(); i++) {
    if (BT::P2PDistance(p, TI->split[i]) <=
        static_cast<Coord>(1 << TI->GetCoverCircle().level)) {
      near_node_idx_seq[near_node_cnt++] = i;
    }
  }
  if (near_node_cnt == 0) {
    return {T, false};
  }

  bool flag = false;
  for (size_t i = 0; i < near_node_cnt; i++) {
    Node* new_node;
    std::tie(new_node, flag) = PointInsertRecursive(
        TI->tree_nodes[near_node_idx_seq[i]], p, TI->GetCoverCircle());
    if (flag) {
      TI->tree_nodes[near_node_idx_seq[i]] = new_node;
      return {T, true};
    }
  }

  assert(flag == false);
  assert(BT::WithinCircle(p, TI->GetCoverCircle()));
  TI->tree_nodes.emplace_back(
      AllocNormalLeafNode<Slice, Leaf>(parlay::make_slice(Points(1, p))));
  TI->split.emplace_back(p);

  return {T, true};
}

}  // namespace pspt

#endif  // PSPT_COVER_BATCH_INSERT_HPP_
