#ifndef PSI_COVER_BATCH_INSERT_HPP_
#define PSI_COVER_BATCH_INSERT_HPP_

#define COVERTREE_TEMPLATE                                              \
  template <typename Point, typename SplitRule, uint_fast8_t kSkHeight, \
            uint_fast8_t kImbaRatio>
#define COVERTREE_CLASS CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>

#include <algorithm>
#include <numeric>
#include <utility>

#include "../cover_tree.h"
#include "../../dependence/tree_node.h"
#include "parlay/slice.h"

namespace psi {

COVERTREE_TEMPLATE
template <typename Range>
void COVERTREE_CLASS::BatchInsert(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);
  Slice A = parlay::make_slice(In);
  Point p;
  p.pnt[0] = 70818;
  p.pnt[1] = 11051;
  std::cout << std::ranges::count(A, p);
  BatchInsert_(A);
}

COVERTREE_TEMPLATE
void COVERTREE_CLASS::BatchInsert_(Slice A) {
  if (this->root_ == nullptr) {
    this->root_ = AllocNormalLeafNode<Slice, Leaf>(A.cut(0, 1));
    this->root_cover_circle_ = CoverCircle{A[0], 0};
    std::cout << "leaf " << A[0] << '\n';
  }

  // int idx = 0;
  for (auto const& p : A.cut(this->root_->size, A.size())) {
    // std::cout << "idx: " << idx++ << p << '\n';
    if (!BT::WithinCircle(p, this->root_cover_circle_)) {
      ExtendCoverRangeUpwards(
          this->root_cover_circle_,
          p);  // PERF: this will change the level of root_cover_circle as well
      // std::cout << "finish extend cover range upwards" << '\n';
    }

    assert(BT::WithinCircle(p, this->root_cover_circle_));
    bool flag = false;
    std::tie(this->root_, flag) =
        PointInsertRecursive(this->root_, p, this->root_cover_circle_);
    assert(flag == true);
  }

  std::cout << "bb: " << BT::GetBox(A).first << " " << BT::GetBox(A).second
            << '\n';
  auto input_circle = BT::template GetCircle<typename BT::NormalCircle>(A);
  std::cout << "reduced: " << input_circle << " "
            << BT::CircleWithinCircle(input_circle, this->root_cover_circle_)
            << '\n';
  std::cout << "root_cover_circle_: " << this->root_cover_circle_ << '\n';

  assert(this->root_ != nullptr);
  return;
}

COVERTREE_TEMPLATE
void CoverTree<Point, SplitRule, kSkHeight,
               kImbaRatio>::ExtendCoverRangeUpwards(CoverCircle& root_cc,
                                                    Point const& p) {
  assert(!BT::WithinCircle(p, root_cc));
  Coord dis = BT::P2PDistanceSquare(p, root_cc.center);
  while (Num::Lt(root_cc.GetRadiusSquare(), dis)) {
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

COVERTREE_TEMPLATE
Node* COVERTREE_CLASS::ShrinkCoverRangeDownwards(
    Node* T, CoverCircle const& level_cover_circle) {
  assert(T->is_leaf);
  assert(T->size == BT::kLeaveWrap);

  auto TL = static_cast<Leaf*>(T);
  auto max_dis = std::numeric_limits<Coord>::lowest();
  // TODO: can replace it with sorting to use the fix thin leave wrap
  parlay::sort_inplace(
      TL->pts.cut(0, TL->size), [&](auto const& p1, auto const& p2) {
        return BT::P2PDistanceSquare(p1, level_cover_circle.center) <
               BT::P2PDistanceSquare(p2, level_cover_circle.center);
      });
  // for (size_t i = 0; i < TL->size; ++i) {
  //   std::cout << "idx: " << i << TL->pts[i] << ' '
  //             << BT::P2PDistanceSquare(TL->pts[i], level_cover_circle.center)
  //             << '\n';
  // }
  // std::cout << '\n';
  // for (size_t i = 0; i < TL->size; ++i) {
  //   max_dis = std::max(
  //       max_dis, BT::P2PDistanceSquare(TL->pts[i],
  //       level_cover_circle.center));
  // }
  max_dis = BT::P2PDistanceSquare(TL->pts[BT::kSlimLeaveWrap - 1],
                                  level_cover_circle.center);
  assert(Num::Leq(max_dis, level_cover_circle.GetRadiusSquare()));
  assert(max_dis);  // TODO: handling the duplicates here
  // std::cout << max_dis << " " << level_cover_circle.level << '\n';

  // NOTE: first find the lowest level that separates the points
  // then build levels in a bottom-up fashion
  auto sep_cover_circle = level_cover_circle;
  while (Num::Leq(max_dis, sep_cover_circle.GetRadiusSquare())) {
    sep_cover_circle.level--;
  }
  assert(sep_cover_circle.center == level_cover_circle.center);
  assert(sep_cover_circle.level < level_cover_circle.level);
  // std::cout << sep_cover_circle.GetRadiusSquare() << " "
  //           << sep_cover_circle.level << '\n';

  // TODO: change it to alloc use thin_leave_wrap
  // Node* node = AllocNormalLeafNode<Slice, Leaf>(
  //     parlay::make_slice(Points{level_cover_circle.GetCenter()}));
  Node* node = AllocNormalLeafNode<Slice, Leaf>(parlay::make_slice(Points()));
  while (sep_cover_circle.level <= level_cover_circle.level) {
    // std::cout << "alloc interior node" << sep_cover_circle << '\n';
    node = AllocInteriorNode<Interior>(
        typename Interior::CoverNodeArr(1, node),
        typename Interior::ST(
            1, level_cover_circle.center),  // NOTE: ensure nesting
        typename Interior::AT{
            .cover_circle = sep_cover_circle,
            .parallel_flag = decltype(Interior::AT::parallel_flag)()});
    sep_cover_circle.level++;  // prepare the new level
  }
  assert(static_cast<Interior*>(node)->GetCoverCircle().level ==
         level_cover_circle.level);
  return node;
}

// BUG: the return value is always true
COVERTREE_TEMPLATE
typename COVERTREE_CLASS::NodeBoolean COVERTREE_CLASS::PointInsertRecursive(
    Node* T, Point const& p, CoverCircle const& level_cover_circle) {
  assert(BT::WithinCircle(p, level_cover_circle));

  if (T->is_leaf) {
    auto TL = static_cast<Leaf*>(T);
    assert(std::find(TL->pts.begin(), TL->pts.begin() + TL->size,
                     level_cover_circle.center) != TL->pts.end());
    if (TL->CapacityFull()) {  // leaf is full
      bool flag = false;
      auto split_node = ShrinkCoverRangeDownwards(T, level_cover_circle);
      assert(static_cast<Interior*>(split_node)->GetCoverCircle() ==
             level_cover_circle);
      // std::cout << "finish shrink cover range downwards" << '\n';
      for (size_t i = 0; i < TL->size; ++i) {
        // std::cout << "leaf_idx: " << i << TL->pts[i] << ' ';
        std::tie(split_node, flag) =
            PointInsertRecursive(split_node, TL->pts[i], level_cover_circle);
        // std::cout << "finish reinsert" << '\n';
        assert(flag == true);
      }
      assert(split_node->size == TL->size);
      FreeNode<Leaf>(T);
      std::tie(split_node, flag) =
          PointInsertRecursive(split_node, p, level_cover_circle);
      // PERF: should not increment the size of the node here, since the point
      // would be either inserted into an interior node or a leaf, in both
      // cases, the size would be updated

      auto TI = static_cast<Interior*>(split_node);
      auto child_size = std::accumulate(
          TI->tree_nodes.begin(), TI->tree_nodes.end(), 0,
          [](auto const& acc, auto const& node) { return acc + node->size; });
      assert(std::cmp_equal(child_size, TI->size));
      assert(TI->size == BT::kLeaveWrap + 1);

      return NodeBoolean(split_node, flag);
    } else {
      return {BT::template InsertPoints2Leaf<Leaf>(
                  T, parlay::make_slice(Points{p})),
              true};
    }
  }

  auto TI = static_cast<Interior*>(T);
  assert(std::ranges::count(TI->split, level_cover_circle.center) == 1);
  assert(TI->GetCoverCircle() == level_cover_circle);

  // find circles that covers the point
  parlay::sequence<size_t> near_node_idx_seq(
      TI->tree_nodes.size());  // TODO: maybe replace using std::views::filter
  size_t near_node_cnt = 0;
  for (size_t i = 0; i < TI->tree_nodes.size(); i++) {
    auto tmp_circle = CoverCircle{TI->split[i], TI->GetCoverCircle().level - 1};
    if (Num::Leq(BT::P2PDistanceSquare(p, TI->split[i]),
                 tmp_circle.GetRadiusSquare())) {
      near_node_idx_seq[near_node_cnt++] = i;
    }
  }

  // try to insert
  bool flag = false;
  for (size_t i = 0; i < near_node_cnt; i++) {
    Node* new_node;
    std::tie(new_node, flag) =
        PointInsertRecursive(TI->tree_nodes[near_node_idx_seq[i]], p,
                             CoverCircle{TI->split[near_node_idx_seq[i]],
                                         TI->GetCoverCircle().level - 1});
    assert(flag == true);
    if (flag) {
      TI->tree_nodes[near_node_idx_seq[i]] = new_node;
      TI->size++;

      auto child_size = std::accumulate(
          TI->tree_nodes.begin(), TI->tree_nodes.end(), 0,
          [](auto const& acc, auto const& node) { return acc + node->size; });
      assert(std::cmp_equal(child_size, TI->size));

      return {T, true};
    }
  }

  // if fails, then insert into the node itself
  assert(flag == false);
  assert(std::ranges::all_of(TI->split, [&](auto const& cc) {
    auto tmp_circle = CoverCircle{cc, TI->GetCoverCircle().level - 1};
    return Num::Gt(BT::P2PDistanceSquare(p, cc), tmp_circle.GetRadiusSquare());
  }));
  TI->tree_nodes.emplace_back(
      AllocNormalLeafNode<Slice, Leaf>(parlay::make_slice(Points{p})));
  TI->split.emplace_back(p);
  TI->size++;

  assert(std::all_of(
      TI->split.begin(), TI->split.end(), [&](auto const& center_i) {
        auto i = &center_i - &TI->split[0];
        auto tmp_circle = CoverCircle{center_i, level_cover_circle.level - 1};
        return std::all_of(TI->split.begin() + i + 1, TI->split.end(),
                           [&](auto const& center_j) {
                             // PERF: should be Gt as a point in a boundary
                             // should falls within that circle
                             return Num::Gt(
                                 BT::P2PDistanceSquare(center_i, center_j),
                                 tmp_circle.GetRadiusSquare());
                           });
      }));

  auto child_size = std::accumulate(
      TI->tree_nodes.begin(), TI->tree_nodes.end(), 0,
      [](auto const& acc, auto const& node) { return acc + node->size; });
  assert(std::cmp_equal(child_size, TI->size));
  return {T, true};
}

}  // namespace psi

#undef COVERTREE_TEMPLATE
#undef COVERTREE_CLASS

#endif
