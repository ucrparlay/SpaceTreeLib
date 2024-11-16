#pragma once

#include <algorithm>
#include <utility>

#include "../base_tree.h"

namespace cpdd {
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf>
size_t BaseTree<Point, DerivedTree, kSkHeight,
                kImbaRatio>::RangeCountRectangleLeaf(Node* T,
                                                     Box const& query_box) {
  assert(T->is_leaf);

  Leaf* TL = static_cast<Leaf*>(T);
  if (TL->is_dummy) {
    return WithinBox(TL->pts[0], query_box) ? TL->size : 0;
  } else {
    return std::ranges::count_if(
        TL->pts.begin(), TL->pts.begin() + TL->size,
        [&](auto const& p) { return WithinBox(p, query_box); });
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsBinaryNode Interior>
  requires(std::same_as<typename Interior::ST,
                        typename BaseTree<Point, DerivedTree, kSkHeight,
                                          kImbaRatio>::HyperPlane>)
size_t BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RangeCountRectangle(
    Node* T, Box const& query_box, Box const& node_box,
    RangeQueryLogger& logger) {
  logger.vis_node_num++;
  if (T->is_leaf) {
    return RangeCountRectangleLeaf<Leaf>(T, query_box);
  }

  Interior* TI = static_cast<Interior*>(T);

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
  logger.generate_box_num++;
  recurse(TI->left, box_cut.GetFirstBoxCut(), l);
  recurse(TI->right, box_cut.GetSecondBoxCut(), r);

  return l + r;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsBinaryNode Interior>
  requires(std::same_as<
           typename Interior::ST,
           typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box>)
size_t BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RangeCountRectangle(
    Node* T, Box const& query_box, RangeQueryLogger& logger) {
  logger.vis_node_num++;
  if (T->is_leaf) {
    return RangeCountRectangleLeaf<Leaf>(T, query_box);
  }

  Interior* TI = static_cast<Interior*>(T);
  auto const& node_box = TI->split;

  if (!BoxIntersectBox(node_box, query_box)) {
    logger.skip_box_num++;
    return 0;
  } else if (WithinBox(node_box, query_box)) {
    logger.full_box_num++;
    return TI->size;
  }

  size_t l, r;
  l = RangeCountRectangle<Leaf, Interior>(TI->left, query_box, logger);
  r = RangeCountRectangle<Leaf, Interior>(TI->right, query_box, logger);

  return l + r;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsMultiNode Interior>
size_t BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RangeCountRectangle(
    Node* T, Box const& query_box, Box const& node_box, DimsType dim,
    BucketType idx, RangeQueryLogger& logger) {
  logger.vis_node_num++;
  if (T->is_leaf) {
    return RangeCountRectangleLeaf<Leaf>(T, query_box);
  }

  Interior* TI = static_cast<Interior*>(T);

  auto recurse = [&query_box, &logger](Node* Ts, Box const& box,
                                       size_t& counter, DimsType next_dim,
                                       DimsType next_idx) -> void {
    if (!BoxIntersectBox(box, query_box)) {
      logger.skip_box_num++;
      counter = 0;
    } else if (WithinBox(box, query_box)) {
      logger.full_box_num++;
      // NOTE: when reach leaf, Ts is the children
      counter = static_cast<Interior*>(Ts)->MergeSize(next_idx);
    } else {
      counter = RangeCountRectangle<Leaf, Interior>(Ts, query_box, box,
                                                    next_dim, next_idx, logger);
    }
  };

  size_t l, r;
  idx <<= 1;
  bool reach_leaf = idx >= Interior::kRegions;
  BoxCut box_cut(node_box, TI->split[dim], true);
  logger.generate_box_num++;

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
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsBinaryNode Interior>
size_t BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RangeCountRadius(
    Node* T, Circle const& cl, Box const& node_box) {
  if (!circle_intersect_box(cl, node_box)) return 0;
  if (within_circle(node_box, cl)) return T->size;

  if (T->is_leaf) {
    size_t cnt = 0;
    Leaf* TL = static_cast<Leaf*>(T);
    if (TL->is_dummy) {
      if (within_circle(TL->pts[0], cl)) {
        cnt += TL->size;
      }
    } else {
      std::ranges::for_each(TL->pts.begin(), TL->pts.begin() + TL->size,
                            [&](auto const& p) {
                              if (within_circle(p, cl)) {
                                cnt++;
                              }
                            });
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

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Range>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RangeQueryLeaf(
    Node* T, Range Out, size_t& s, Box const& query_box) {
  assert(T->is_leaf);

  Leaf* TL = static_cast<Leaf*>(T);
  if (TL->is_dummy) {
    if (WithinBox(TL->pts[0], query_box)) {
      std::fill_n(Out.begin() + s, TL->size, TL->pts[0]);
      s += TL->size;
    }
  } else {
    auto result = std::ranges::copy_if(
        TL->pts.begin(), TL->pts.begin() + TL->size, Out.begin() + s,
        [&](auto const& p) { return WithinBox(p, query_box); });
    s += std::ranges::distance(Out.begin() + s, result.out);
  }
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsBinaryNode Interior, typename Range>
  requires(std::same_as<typename Interior::ST,
                        typename BaseTree<Point, DerivedTree, kSkHeight,
                                          kImbaRatio>::HyperPlane>)
void BaseTree<Point, DerivedTree, kSkHeight,
              kImbaRatio>::RangeQuerySerialRecursive(Node* T, Range Out,
                                                     size_t& s,
                                                     Box const& query_box,
                                                     Box const& node_box,
                                                     RangeQueryLogger& logger) {
  logger.vis_node_num++;
  if (T->is_leaf) {
    RangeQueryLeaf<Leaf>(T, Out, s, query_box);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);

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
  logger.generate_box_num++;
  recurse(TI->left, box_cut.GetFirstBoxCut());
  recurse(TI->right, box_cut.GetSecondBoxCut());

  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsBinaryNode Interior, typename Range>
  requires(std::same_as<
           typename Interior::ST,
           typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box>)
void BaseTree<Point, DerivedTree, kSkHeight,
              kImbaRatio>::RangeQuerySerialRecursive(Node* T, Range Out,
                                                     size_t& s,
                                                     Box const& query_box,
                                                     RangeQueryLogger& logger) {
  logger.vis_node_num++;
  if (T->is_leaf) {
    RangeQueryLeaf<Leaf>(T, Out, s, query_box);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);
  auto const& node_box = TI->split;

  if (!BoxIntersectBox(node_box, query_box)) {
    logger.skip_box_num++;
    return;
  } else if (WithinBox(node_box, query_box)) {
    logger.full_box_num++;
    FlattenRec<Leaf, Interior>(T, Out.cut(s, s + T->size));
    s += T->size;
    return;
  }

  RangeQuerySerialRecursive<Leaf, Interior>(TI->left, Out, s, query_box,
                                            logger);
  RangeQuerySerialRecursive<Leaf, Interior>(TI->right, Out, s, query_box,
                                            logger);

  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsMultiNode Interior, typename Range>
void BaseTree<Point, DerivedTree, kSkHeight,
              kImbaRatio>::RangeQuerySerialRecursive(Node* T, Range Out,
                                                     size_t& s,
                                                     Box const& query_box,
                                                     Box const& node_box,
                                                     DimsType dim,
                                                     BucketType idx,
                                                     RangeQueryLogger& logger) {
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
      size_t candidate_size = static_cast<Interior*>(Ts)->MergeSize(next_idx);
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
