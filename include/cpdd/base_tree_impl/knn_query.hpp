#pragma once
#include <algorithm>
#include <tuple>
#include <utility>

#include "../base_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {

// NOTE: distance between two Points
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Coord
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::P2PDistance(
    Point const& p, Point const& q) {
  Coord r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
    r += (p.pnt[i] - q.pnt[i]) * (p.pnt[i] - q.pnt[i]);
  }
  return r;
}

// NOTE: Distance between a Point and a Box
// return 0 when p is inside the box a
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Coord
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::P2BMinDistance(
    Point const& p,
    typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box const&
        a) {
  Coord r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(p.pnt[i], a.first.pnt[i])) {
      r += (a.first.pnt[i] - p.pnt[i]) * (a.first.pnt[i] - p.pnt[i]);
    } else if (Num::Gt(p.pnt[i], a.second.pnt[i])) {
      r += (p.pnt[i] - a.second.pnt[i]) * (p.pnt[i] - a.second.pnt[i]);
    }
  }
  return r;
}

// NOTE: Max distance between a Point and a Box
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Coord
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::P2BMaxDistance(
    Point const& p,
    typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box const&
        a) {
  Coord r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(p.pnt[i], (a.second.pnt[i] + a.first.pnt[i]) / 2)) {
      r += (a.second.pnt[i] - p.pnt[i]) * (a.second.pnt[i] - p.pnt[i]);
    } else {
      r += (p.pnt[i] - a.first.pnt[i]) * (p.pnt[i] - a.first.pnt[i]);
    }
  }
  return r;
}

// NOTE: early return the partial distance between p and q if it is larger than
// r else return the distance between p and q
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Coord
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::InterruptibleDistance(
    Point const& p, Point const& q, Coord up) {
  Coord r = 0;
  DimsType i = 0;
  if (kDim >= 6) {
    while (1) {
      r += (p.pnt[i] - q.pnt[i]) * (p.pnt[i] - q.pnt[i]);
      ++i;
      r += (p.pnt[i] - q.pnt[i]) * (p.pnt[i] - q.pnt[i]);
      ++i;
      r += (p.pnt[i] - q.pnt[i]) * (p.pnt[i] - q.pnt[i]);
      ++i;
      r += (p.pnt[i] - q.pnt[i]) * (p.pnt[i] - q.pnt[i]);
      ++i;

      if (Num::Gt(r, up)) {
        return r;
      }
      if (i + 4 > kDim) {
        break;
      }
    }
  }
  while (i < kDim) {
    r += (p.pnt[i] - q.pnt[i]) * (p.pnt[i] - q.pnt[i]);
    ++i;
  }
  return r;
}

// NOTE: KNN search for Point q
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Range>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::KNNLeaf(
    Node* T, Point const& q, kBoundedQueue<Point, Range>& bq) {
  assert(T->is_leaf);

  Leaf* TL = static_cast<Leaf*>(T);
  size_t i = 0;
  while (!bq.full() && i < TL->size) {
    bq.insert(std::make_pair(std::ref(TL->pts[(!TL->is_dummy) * i]),
                             P2PDistance(q, TL->pts[(!TL->is_dummy) * i])));
    i++;
  }
  while (i < TL->size) {
    Coord r =
        InterruptibleDistance(q, TL->pts[(!TL->is_dummy) * i], bq.top_value());
    if (Num::Lt(r, bq.top_value())) {
      bq.insert(std::make_pair(std::ref(TL->pts[(!TL->is_dummy) * i]), r));
    } else if (TL->is_dummy) {
      break;
    }
    i++;
  }
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsBinaryNode Interior, typename Range>
  requires(std::same_as<typename Interior::ST,
                        typename BaseTree<Point, DerivedTree, kSkHeight,
                                          kImbaRatio>::HyperPlane>)
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::KNNBinary(
    Node* T, Point const& q, kBoundedQueue<Point, Range>& bq,
    Box const& node_box, KNNLogger& logger) {
  logger.vis_node_num++;

  if (T->is_leaf) {
    KNNLeaf<Leaf>(T, q, bq);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);
  bool go_left = Num::Gt(TI->split.first - q.pnt[TI->split.second], 0);
  BoxCut box_cut(node_box, TI->split, go_left);
  logger.generate_box_num += 1;

  KNNBinary<Leaf, Interior>(go_left ? TI->left : TI->right, q, bq,
                            box_cut.GetFirstBoxCut(), logger);

  logger.check_box_num++;
  if (Num::Gt(P2BMinDistance(q, box_cut.GetSecondBoxCut()), bq.top_value()) &&
      bq.full()) {
    logger.skip_box_num++;
    return;
  }
  KNNBinary<Leaf, Interior>(go_left ? TI->right : TI->left, q, bq,
                            box_cut.GetBox(), logger);
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsBinaryNode Interior, typename Range>
  requires(std::same_as<
           typename Interior::ST,
           typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box>)
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::KNNBinary(
    Node* T, Point const& q, kBoundedQueue<Point, Range>& bq,
    KNNLogger& logger) {
  logger.vis_node_num++;

  if (T->is_leaf) {
    KNNLeaf<Leaf>(T, q, bq);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);
  Coord dist_left =
      TI->left->is_leaf
          ? static_cast<Coord>(0)
          : P2BMinDistance(q, static_cast<Interior*>(TI->left)->split);
  Coord dist_right =
      TI->right->is_leaf
          ? static_cast<Coord>(0)
          : P2BMinDistance(q, static_cast<Interior*>(TI->right)->split);
  bool go_left = Num::Leq(dist_left, dist_right);
  logger.generate_box_num += 1;

  KNNBinary<Leaf, Interior>(go_left ? TI->left : TI->right, q, bq, logger);

  logger.check_box_num++;
  if (Num::Gt(go_left ? dist_right : dist_left, bq.top_value()) && bq.full()) {
    logger.skip_box_num++;
    return;
  }
  KNNBinary<Leaf, Interior>(go_left ? TI->right : TI->left, q, bq, logger);
  return;
}

// NOTE: compute knn for multinode as if a binary node
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsMultiNode Interior, typename Range>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::KNNMultiExpand(
    Node* T, Point const& q, DimsType dim, BucketType idx,
    kBoundedQueue<Point, Range>& bq, Box const& node_box, KNNLogger& logger) {
  logger.vis_node_num++;

  if (T->size == 0) {
    return;
  }

  if (T->is_leaf) {
    KNNLeaf<Leaf>(T, q, bq);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);
  bool go_left =
      Num::Gt(TI->split[dim].first - q.pnt[TI->split[dim].second], 0);

  BucketType first_idx = (idx << 1) + static_cast<BucketType>(!go_left);
  BucketType second_idx = (idx << 1) + static_cast<BucketType>(go_left);
  bool reach_leaf = first_idx >= Interior::kRegions;
  Node* first_node =
      reach_leaf ? TI->tree_nodes[first_idx - Interior::kRegions] : T;
  Node* second_node =
      reach_leaf ? TI->tree_nodes[second_idx - Interior::kRegions] : T;
  if (reach_leaf) {
    first_idx = second_idx = 1;
  }

  BoxCut box_cut(node_box, TI->split[dim], go_left);
  logger.generate_box_num += 1;
  assert((dim + 1) % kDim != 0 || (first_idx == 1 && second_idx == 1));

  KNNMultiExpand<Leaf, Interior>(first_node, q, (dim + 1) % kDim, first_idx, bq,
                                 box_cut.GetFirstBoxCut(), logger);

  // NOTE: compute the other bounding box
  logger.check_box_num++;
  if (Num::Gt(P2BMinDistance(q, box_cut.GetSecondBoxCut()), bq.top_value()) &&
      bq.full()) {
    logger.skip_box_num++;
    return;
  }
  KNNMultiExpand<Leaf, Interior>(second_node, q, (dim + 1) % kDim, second_idx,
                                 bq, box_cut.GetBox(), logger);
  return;
}

// NOTE: compute KNN for multi-node by computing bounding boxes
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsMultiNode Interior, typename Range>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::KNNMulti(
    Node* T, Point const& q, kBoundedQueue<Point, Range>& bq,
    Box const& node_box, KNNLogger& logger) {
  logger.vis_node_num++;

  if (T->size == 0) {
    return;
  }

  if (T->is_leaf) {
    KNNLeaf<Leaf>(T, q, bq);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);

  BoxSeq regions(Interior::kRegions);
  std::array<std::pair<Coord, BucketType>, Interior::kRegions> dists;
  TI->ComputeSubregions(
      regions, node_box, 1,
      0);  // PERF: find better way to compute the bounding boxes
  logger.generate_box_num += Interior::kRegions * 2 - 1;

  std::ranges::generate(dists, [i = 0, &q, &regions]() mutable {
    auto r = std::make_pair(P2BMinDistance(q, regions[i]), i);
    i++;
    return r;
  });
  std::ranges::sort(dists, std::less<>(),
                    [&](auto const& box_pair) { return box_pair.first; });

  KNNMulti<Leaf, Interior>(TI->tree_nodes[dists[0].second], q, bq,
                           regions[dists[0].second], logger);
  for (BucketType i = 1; i < Interior::kRegions; ++i) {
    logger.check_box_num++;
    if (Num::Gt(dists[i].first, bq.top_value()) && bq.full()) {
      logger.skip_box_num++;
      continue;
    }
    KNNMulti<Leaf, Interior>(TI->tree_nodes[dists[i].second], q, bq,
                             regions[dists[i].second], logger);
  }

  return;
}

}  // namespace cpdd
