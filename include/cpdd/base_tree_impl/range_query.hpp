#pragma once

#include <utility>

#include "../base_tree.h"

namespace cpdd {
template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
template <typename Leaf>
size_t BaseTree<Point, DerivedTree, kBDO>::RangeCountRectangleLeaf(
    Node* T, Box const& query_box) {
  assert(T->is_leaf);

  Leaf* TL = static_cast<Leaf*>(T);
  size_t cnt = 0;
  if (TL->is_dummy) {
    if (WithinBox(TL->pts[0], query_box)) {
      cnt = TL->size;
    }
  } else {
    for (size_t i = 0; i < TL->size; i++) {
      if (WithinBox(TL->pts[i], query_box)) {
        cnt++;
      }
    }
  }
  return cnt;
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
template <typename Leaf, IsBinaryNode Interior>
size_t BaseTree<Point, DerivedTree, kBDO>::RangeCountRectangle(
    Node* T, Box const& query_box, Box const& node_box,
    RangeQueryLogger& logger) {
  logger.vis_node_num++;
  if (T->is_leaf) {
    return RangeCountRectangleLeaf<Leaf>(T, query_box);
  }

  Interior* TI = static_cast<Interior*>(T);
  logger.generate_box_num++;

  size_t l, r;
  auto recurse = [&](Node* Ts, Box const& box, size_t& counter) -> void {
    if (!BoxIntersectBox(box, query_box)) {
      logger.skip_box_num++;
      counter = 0;
    } else if (WithinBox(box, query_box)) {
      logger.full_box_num++;
      counter = Ts->size;
    } else {
      counter = RangeCountRectangle<Leaf, Interior>(Ts, query_box, box, logger);
    }
  };

  BoxCut box_cut(node_box, TI->split, true);
  recurse(TI->left, box_cut.GetFirstBoxCut(), l);
  recurse(TI->right, box_cut.GetSecondBoxCut(), r);

  return l + r;
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
template <typename Leaf, IsMultiNode Interior>
size_t BaseTree<Point, DerivedTree, kBDO>::RangeCountRectangle(
    Node* T, Box const& query_box, Box const& node_box, DimsType dim,
    BucketType idx, RangeQueryLogger& logger) {
  logger.vis_node_num++;
  if (T->is_leaf) {
    return RangeCountRectangleLeaf<Leaf>(T, query_box);
  }

  Interior* TI = static_cast<Interior*>(T);
  logger.generate_box_num++;

  auto recurse = [&query_box, &logger](Node* Ts, Box const& box,
                                       size_t& counter, DimsType next_dim,
                                       DimsType next_idx) -> void {
    if (!BoxIntersectBox(box, query_box)) {
      logger.skip_box_num++;
      counter = 0;
    } else if (WithinBox(box, query_box)) {
      logger.full_box_num++;
      // NOTE: when reach leaf, Ts is the children
      counter = static_cast<Interior*>(Ts)->ReduceSums(next_idx);
    } else {
      counter = RangeCountRectangle<Leaf, Interior>(Ts, query_box, box,
                                                    next_dim, next_idx, logger);
    }
  };

  size_t l, r;
  idx <<= 1;
  bool reach_leaf = idx >= Interior::kRegions;
  BoxCut box_cut(node_box, TI->split[dim], true);

  // NOTE: visit left half
  assert(TI->split[dim].second == dim);

  recurse(reach_leaf ? TI->tree_nodes[idx - Interior::kRegions] : T,
          box_cut.GetFirstBoxCut(), l, (dim + 1) % kDim, reach_leaf ? 1 : idx);

  idx |= 1;
  recurse(reach_leaf ? TI->tree_nodes[idx - Interior::kRegions] : T,
          box_cut.GetSecondBoxCut(), r, (dim + 1) % kDim, reach_leaf ? 1 : idx);

  return l + r;
}

// TODO: as range_count_rectangle
template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
template <typename Leaf, IsBinaryNode Interior>
size_t BaseTree<Point, DerivedTree, kBDO>::RangeCountRadius(
    Node* T, Circle const& cl, Box const& node_box) {
  if (!circle_intersect_box(cl, node_box)) return 0;
  if (within_circle(node_box, cl)) return T->size;

  if (T->is_leaf) {
    size_t cnt = 0;
    Leaf* TL = static_cast<Leaf*>(T);
    for (int i = 0; i < TL->size; i++) {
      if (within_circle(TL->pts[(!TL->is_dummy) * i], cl)) {
        cnt++;
      }
    }
    return cnt;
  }

  Interior* TI = static_cast<Interior*>(T);
  Box lbox(node_box), rbox(node_box);
  lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
  rbox.first.pnt[TI->split.second] = TI->split.first;

  size_t l, r;
  parlay::par_do_if(
      TI->size >= kSerialBuildCutoff,
      [&] { l = RangeCountRadius(TI->left, cl, lbox); },
      [&] { r = RangeCountRadius(TI->right, cl, rbox); });

  return l + r;
};

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
template <typename Leaf, typename Range>
void BaseTree<Point, DerivedTree, kBDO>::RangeQueryLeaf(Node* T, Range Out,
                                                        size_t& s,
                                                        Box const& query_box) {
  assert(T->is_leaf);

  Leaf* TL = static_cast<Leaf*>(T);
  if (TL->is_dummy) {
    if (WithinBox(TL->pts[0], query_box)) {
      for (size_t i = 0; i < TL->size; i++) {
        Out[s++] = TL->pts[0];
      }
    }
  } else {
    for (size_t i = 0; i < TL->size; i++)
      if (WithinBox(TL->pts[i], query_box)) {
        Out[s++] = TL->pts[i];
      }
  }
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
template <typename Leaf, IsBinaryNode Interior, typename Range>
void BaseTree<Point, DerivedTree, kBDO>::RangeQuerySerialRecursive(
    Node* T, Range Out, size_t& s, Box const& query_box, Box const& node_box,
    RangeQueryLogger& logger) {
  logger.vis_node_num++;
  if (T->is_leaf) {
    RangeQueryLeaf<Leaf>(T, Out, s, query_box);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);
  logger.generate_box_num++;

  auto recurse = [&](Node* Ts, Box const& box) -> void {
    if (!BoxIntersectBox(box, query_box)) {
      logger.skip_box_num++;
      return;
    } else if (WithinBox(box, query_box)) {
      logger.full_box_num++;
      FlattenRec<Leaf, Interior>(Ts, Out.cut(s, s + Ts->size));
      s += Ts->size;
      return;
    } else {
      RangeQuerySerialRecursive<Leaf, Interior>(Ts, Out, s, query_box, box,
                                                logger);
      return;
    }
  };

  BoxCut box_cut(node_box, TI->split, true);
  recurse(TI->left, box_cut.GetFirstBoxCut());
  recurse(TI->right, box_cut.GetSecondBoxCut());

  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
template <typename Leaf, IsMultiNode Interior, typename Range>
void BaseTree<Point, DerivedTree, kBDO>::RangeQuerySerialRecursive(
    Node* T, Range Out, size_t& s, Box const& query_box, Box const& node_box,
    DimsType dim, BucketType idx, RangeQueryLogger& logger) {
  logger.vis_node_num++;
  if (T->is_leaf) {
    RangeQueryLeaf<Leaf>(T, Out, s, query_box);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);

  auto recurse = [&query_box, &s, &Out, &logger](Node* Ts, Box const& box,
                                                 DimsType next_dim,
                                                 BucketType next_idx) -> void {
    if (!BoxIntersectBox(box, query_box)) {
      logger.skip_box_num++;
      return;
    } else if (WithinBox(box, query_box)) {
      logger.full_box_num++;
      size_t candidate_size = static_cast<Interior*>(Ts)->ReduceSums(next_idx);
      PartialFlatten<Leaf, Interior>(Ts, Out.cut(s, s + candidate_size),
                                     next_idx);
      s += candidate_size;
      return;
    } else {
      RangeQuerySerialRecursive<Leaf, Interior>(Ts, Out, s, query_box, box,
                                                next_dim, next_idx, logger);
      return;
    }
  };

  idx <<= 1;
  bool reach_leaf = idx >= Interior::kRegions;
  logger.generate_box_num++;
  BoxCut box_cut(node_box, TI->split[dim], true);

  recurse(reach_leaf ? TI->tree_nodes[idx - Interior::kRegions] : T,
          box_cut.GetFirstBoxCut(), (dim + 1) % kDim, reach_leaf ? 1 : idx);

  // NOTE: visit right
  idx |= 1;
  recurse(reach_leaf ? TI->tree_nodes[idx - Interior::kRegions] : T,
          box_cut.GetSecondBoxCut(), (dim + 1) % kDim, reach_leaf ? 1 : idx);

  return;
}
}  // namespace cpdd
