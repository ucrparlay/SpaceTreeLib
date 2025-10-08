#ifndef PSI_BASE_TREE_IMPL_KNN_QUERY_HPP
#define PSI_BASE_TREE_IMPL_KNN_QUERY_HPP

#include <algorithm>
#include <utility>

#include "../base_tree.h"
#include "parlay/primitives.h"
#include "psi/dependence/tree_node.h"

#define BASETREE_TEMPLATE                                                 \
  template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight, \
            uint_fast8_t kImbaRatio>
#define BASETREE_CLASS BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>

namespace psi {

// NOTE: distance between two Points
// TODO: change the name to P2PDistanceSquare to avoid ambiguous
BASETREE_TEMPLATE
inline typename BASETREE_CLASS::DisType BASETREE_CLASS::P2PDistanceSquare(
    Point const& p, Point const& q) {
  constexpr uint_fast8_t kDim = Point::GetDim();
  DisType r = 0;

  if constexpr (kDim == 2) {
    DisType x = static_cast<DisType>(p.pnt[0]) - static_cast<DisType>(q.pnt[0]);
    DisType y = static_cast<DisType>(p.pnt[1]) - static_cast<DisType>(q.pnt[1]);
    r = x * x + y * y;
  } else if constexpr (kDim == 3) {
    DisType x = static_cast<DisType>(p.pnt[0]) - static_cast<DisType>(q.pnt[0]);
    DisType y = static_cast<DisType>(p.pnt[1]) - static_cast<DisType>(q.pnt[1]);
    DisType z = static_cast<DisType>(p.pnt[2]) - static_cast<DisType>(q.pnt[2]);
    r = x * x + y * y + z * z;
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
    }
  }
  return r;
}

// NOTE: Distance between a Point and a Box
// return 0 when p is inside the box a
BASETREE_TEMPLATE
inline typename BASETREE_CLASS::DisType BASETREE_CLASS::P2BMinDistanceSquare(
    Point const& p, typename BASETREE_CLASS::Box const& a) {
  DisType r = 0;
  // NOTE: the distance is 0 when p is inside the box
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(p.pnt[i], a.first.pnt[i])) {
      r += (static_cast<DisType>(a.first.pnt[i]) -
            static_cast<DisType>(p.pnt[i])) *
           (static_cast<DisType>(a.first.pnt[i]) -
            static_cast<DisType>(p.pnt[i]));
    } else if (Num::Gt(p.pnt[i], a.second.pnt[i])) {
      r += (static_cast<DisType>(p.pnt[i]) -
            static_cast<DisType>(a.second.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) -
            static_cast<DisType>(a.second.pnt[i]));
    } else {  // will not count the dis if p is inside the box in dimension i
      ;
    }
  }
  return r;
}

// NOTE: Max distance between a Point and a Box
BASETREE_TEMPLATE
inline typename BASETREE_CLASS::DisType BASETREE_CLASS::P2BMaxDistanceSquare(
    Point const& p, typename BASETREE_CLASS::Box const& a) {
  DisType r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(p.pnt[i], (a.second.pnt[i] + a.first.pnt[i]) / 2)) {
      r += (static_cast<DisType>(a.second.pnt[i]) -
            static_cast<DisType>(p.pnt[i])) *
           (static_cast<DisType>(a.second.pnt[i]) -
            static_cast<DisType>(p.pnt[i]));
    } else {
      r += (static_cast<DisType>(p.pnt[i]) -
            static_cast<DisType>(a.first.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) -
            static_cast<DisType>(a.first.pnt[i]));
    }
  }
  return r;
}

BASETREE_TEMPLATE
inline double BaseTree<Point, DerivedTree, kSkHeight,
                       kImbaRatio>::P2CMinDistance(Point const& p,
                                                   Point const& center,
                                                   DisType const r) {
  // return Num_Comparator<double>::Max(
  //     0.0, std::sqrt(P2PDistanceSquare(p, center)) - static_cast<double>(r));
  return std::sqrt(P2PDistanceSquare(p, center)) - static_cast<double>(r);
}

BASETREE_TEMPLATE
template <typename CircleType>
inline double BaseTree<Point, DerivedTree, kSkHeight,
                       kImbaRatio>::P2CMinDistance(Point const& p,
                                                   CircleType const& cl) {
  return P2CMinDistance(p, cl.GetCenter(), cl.GetRadius());
}

// NOTE: early return the partial distance between p and q if it is larger than
// r else return the distance between p and q
BASETREE_TEMPLATE
inline typename BASETREE_CLASS::DisType BASETREE_CLASS::InterruptibleDistance(
    Point const& p, Point const& q, DisType up) {
  DisType r = 0;
  DimsType i = 0;
  if (kDim >= 6) {
    while (1) {
      r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
      ++i;
      r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
      ++i;
      r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
      ++i;
      r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
      ++i;

      if (Num_Comparator<DisType>::Gt(r, up)) {
        return r;
      }
      if (i + 4 > kDim) {
        break;
      }
    }
  }
  while (i < kDim) {
    r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
         (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
    ++i;
  }
  return r;
}

// NOTE: KNN search for Point q
BASETREE_TEMPLATE
template <typename Leaf, typename Range>
void BASETREE_CLASS::KNNLeaf(Node* T, Point const& q,
                             kBoundedQueue<Point, Range>& bq) {
  assert(T->is_leaf);

  Leaf* TL = static_cast<Leaf*>(T);
  size_t i = 0;
  while (!bq.full() && i < TL->size) {
    bq.insert(
        std::make_pair(std::ref(TL->pts[(!TL->is_dummy) * i]),
                       P2PDistanceSquare(q, TL->pts[(!TL->is_dummy) * i])));
    i++;
  }
  while (i < TL->size) {
    auto r =
        InterruptibleDistance(q, TL->pts[(!TL->is_dummy) * i], bq.top_value());
    if (Num::Lt(r, bq.top_value())) {  // PERF: the queue is full, no need to
                                       // insert points with equal distances
      bq.insert(std::make_pair(std::ref(TL->pts[(!TL->is_dummy) * i]), r));
    } else if (TL->is_dummy) {
      break;
    }
    i++;
  }
  return;
}

BASETREE_TEMPLATE
template <typename Leaf, IsBinaryNode Interior, typename Range>
void BASETREE_CLASS::KNNBinary(Node* T, Point const& q,
                               kBoundedQueue<Point, Range>& bq,
                               Box const& node_box, KNNLogger& logger)
  requires std::same_as<typename Interior::ST,
                        typename BASETREE_CLASS::HyperPlane>
{
  if (T->is_leaf) {
    logger.vis_leaf_num++;
    KNNLeaf<Leaf>(T, q, bq);
    return;
  }

  logger.vis_interior_num++;
  Interior* TI = static_cast<Interior*>(T);
  bool go_left = Num::Gt(TI->split.first - q.pnt[TI->split.second], 0);
  BoxCut box_cut(node_box, TI->split, go_left);
  logger.generate_box_num += 1;

  KNNBinary<Leaf, Interior>(go_left ? TI->left : TI->right, q, bq,
                            box_cut.GetFirstBoxCut(), logger);

  logger.check_box_num++;
  if (Num::Gt(P2BMinDistanceSquare(q, box_cut.GetSecondBoxCut()),
              bq.top_value()) &&
      bq.full()) {
    logger.skip_box_num++;
    return;
  }
  KNNBinary<Leaf, Interior>(go_left ? TI->right : TI->left, q, bq,
                            box_cut.GetBox(), logger);
  return;
}

BASETREE_TEMPLATE
template <typename Leaf, IsBinaryNode Interior, typename Range>
void BASETREE_CLASS::KNNBinaryBox(Node* T, Point const& q,
                                  kBoundedQueue<Point, Range>& bq,
                                  KNNLogger& logger) {
  if (bq.size() &&
      Num::Gt(P2BMinDistanceSquare(q, RetriveBox<Leaf, Interior>(T)),
              bq.top_value()) &&
      bq.full()) {
    logger.skip_box_num++;
    return;
  }

  if (T->is_leaf) {
    logger.vis_leaf_num++;
    KNNLeaf<Leaf>(T, q, bq);
    return;
  }

  logger.vis_interior_num++;
  Interior* TI = static_cast<Interior*>(T);
  Coord dist_left =
      P2BMinDistanceSquare(q, RetriveBox<Leaf, Interior>(TI->left));
  Coord dist_right =
      P2BMinDistanceSquare(q, RetriveBox<Leaf, Interior>(TI->right));
  bool go_left = Num::Leq(dist_left, dist_right);

  KNNBinaryBox<Leaf, Interior>(go_left ? TI->left : TI->right, q, bq, logger);

  logger.check_box_num++;
  if (Num::Gt(go_left ? dist_right : dist_left, bq.top_value()) && bq.full()) {
    logger.skip_box_num++;
    return;
  }
  KNNBinaryBox<Leaf, Interior>(go_left ? TI->right : TI->left, q, bq, logger);
  return;
}

// NOTE: compute knn for multinode as if a binary node
BASETREE_TEMPLATE
template <typename Leaf, IsMultiNode Interior, typename Range>
void BASETREE_CLASS::KNNMultiExpand(Node* T, Point const& q, DimsType dim,
                                    BucketType idx,
                                    kBoundedQueue<Point, Range>& bq,
                                    Box const& node_box, KNNLogger& logger) {
  if (T->size == 0) {
    return;
  }

  if (T->is_leaf) {
    logger.vis_leaf_num++;
    KNNLeaf<Leaf>(T, q, bq);
    return;
  }

  logger.vis_interior_num++;

  Interior* TI = static_cast<Interior*>(T);
  auto const& split = Interior::EqualSplit() ? TI->split[dim] : TI->split[idx];
  bool go_left = Num::Gt(split.first - q.pnt[split.second], 0);

  BucketType first_idx = (idx << 1) + static_cast<BucketType>(!go_left);
  BucketType second_idx = (idx << 1) + static_cast<BucketType>(go_left);
  bool reach_leaf =
      first_idx >= Interior::GetRegions();  // whether reach the skeleton leaf
  Node* first_node =
      reach_leaf ? TI->tree_nodes[first_idx - Interior::GetRegions()] : T;
  Node* second_node =
      reach_leaf ? TI->tree_nodes[second_idx - Interior::GetRegions()] : T;
  if (reach_leaf) {
    first_idx = second_idx = 1;
  }

  BoxCut box_cut(node_box, split, go_left);
  logger.generate_box_num += 1;
  assert((dim + 1) % kDim != 0 || (first_idx == 1 && second_idx == 1));

  KNNMultiExpand<Leaf, Interior>(first_node, q, (dim + 1) % kDim, first_idx, bq,
                                 box_cut.GetFirstBoxCut(), logger);

  // NOTE: compute the other bounding box
  logger.check_box_num++;
  if (Num::Gt(P2BMinDistanceSquare(q, box_cut.GetSecondBoxCut()),
              bq.top_value()) &&
      bq.full()) {
    logger.skip_box_num++;
    return;
  }
  KNNMultiExpand<Leaf, Interior>(second_node, q, (dim + 1) % kDim, second_idx,
                                 bq, box_cut.GetBox(), logger);
  return;
}

// NOTE: compute knn for multinode as if a binary node
BASETREE_TEMPLATE
template <typename Leaf, IsMultiNode Interior, typename Range>
void BASETREE_CLASS::KNNMultiExpandBox(Node* T, Point const& q, DimsType dim,
                                       BucketType idx,
                                       kBoundedQueue<Point, Range>& bq,
                                       KNNLogger& logger) {
  if (T->size == 0) {
    return;
  }

  if (idx == 1 && bq.size() &&
      Num::Gt(P2BMinDistanceSquare(q, RetriveBox<Leaf, Interior>(T)),
              bq.top_value()) &&
      bq.full()) {
    logger.skip_box_num++;
    return;
  }
  if (T->is_leaf) {
    logger.vis_leaf_num++;
    KNNLeaf<Leaf>(T, q, bq);
    return;
  }

  logger.vis_interior_num++;

  Interior* TI = static_cast<Interior*>(T);
  auto const& split = Interior::EqualSplit() ? TI->split[dim] : TI->split[idx];
  bool go_left = Num::Gt(split.first - q.pnt[split.second], 0);

  BucketType first_idx = (idx << 1) + static_cast<BucketType>(!go_left);
  BucketType second_idx = (idx << 1) + static_cast<BucketType>(go_left);
  bool reach_leaf =
      first_idx >= Interior::GetRegions();  // whether reach the skeleton leaf
  Node* first_node =
      reach_leaf ? TI->tree_nodes[first_idx - Interior::GetRegions()] : T;
  Node* second_node =
      reach_leaf ? TI->tree_nodes[second_idx - Interior::GetRegions()] : T;
  if (reach_leaf) {
    first_idx = second_idx = 1;
  }
  auto second_box =
      second_node->is_leaf
          ? RetriveBox<Leaf, Interior>(second_node)
          : static_cast<Interior*>(second_node)->GetBoxById(second_idx);

  assert((dim + 1) % kDim != 0 || (first_idx == 1 && second_idx == 1));

  KNNMultiExpandBox<Leaf, Interior>(first_node, q, (dim + 1) % kDim, first_idx,
                                    bq, logger);

  // NOTE: compute the other bounding box
  logger.check_box_num++;
  if (Num::Gt(P2BMinDistanceSquare(q, second_box), bq.top_value()) &&
      bq.full()) {
    logger.skip_box_num++;
    return;
  }
  KNNMultiExpandBox<Leaf, Interior>(second_node, q, (dim + 1) % kDim,
                                    second_idx, bq, logger);
  return;
}

// NOTE: compute KNN for multi-node by computing bounding boxes
BASETREE_TEMPLATE
template <typename Leaf, IsMultiNode Interior, typename Range>
void BASETREE_CLASS::KNNMulti(Node* T, Point const& q,
                              kBoundedQueue<Point, Range>& bq,
                              KNNLogger& logger) {
  if (T->size == 0) {
    return;
  }

  if (bq.size() &&
      Num::Gt(P2BMinDistanceSquare(q, RetriveBox<Leaf, Interior>(T)),
              bq.top_value()) &&
      bq.full()) {
    logger.skip_box_num++;
    return;
  }
  if (T->is_leaf) {
    logger.vis_leaf_num++;
    KNNLeaf<Leaf>(T, q, bq);
    return;
  }

  logger.vis_interior_num++;
  Interior* TI = static_cast<Interior*>(T);

  // TODO: may change to use parlay::sequence
  std::array<std::pair<Coord, BucketType>, Interior::GetRegions()> dists;

  std::ranges::generate(dists, [i = 0, &q, TI]() mutable {
    auto r = std::make_pair(
        P2BMinDistanceSquare(q, RetriveBox<Leaf, Interior>(TI->tree_nodes[i])),
        i);
    i++;
    return r;
  });
  std::ranges::sort(dists, std::less<>(),
                    [&](auto const& box_pair) { return box_pair.first; });

  KNNMulti<Leaf, Interior>(TI->tree_nodes[dists[0].second], q, bq, logger);
  for (BucketType i = 1; i < Interior::GetRegions(); ++i) {
    logger.check_box_num++;
    if (Num::Gt(dists[i].first, bq.top_value()) && bq.full()) {
      logger.skip_box_num++;
      continue;
    }
    KNNMulti<Leaf, Interior>(TI->tree_nodes[dists[i].second], q, bq, logger);
  }

  return;
}

BASETREE_TEMPLATE
template <typename Leaf, IsBinaryNode BN, IsMultiNode MN, typename Range>
void BASETREE_CLASS::KNNMix(Node* T, Point const& q, DimsType dim,
                            BucketType idx, kBoundedQueue<Point, Range>& bq,
                            Box const& node_box, KNNLogger& logger) {
  if (T->size == 0) {
    return;
  }

  if (T->is_leaf) {
    logger.vis_leaf_num++;
    KNNLeaf<Leaf>(T, q, bq);
    return;
  }

  logger.vis_interior_num++;
  bool go_left(false);
  HyperPlane* split(nullptr);
  Node *first_node(nullptr), *second_node(nullptr);
  BucketType first_idx(0), second_idx(0);
  if (dynamic_cast<BN*>(T)) {  // NOTE: Binary node
    auto TI = static_cast<BN*>(T);
    assert(TI->GetParallelFlagIniStatus() && !TI->aug.value());
    split = &(TI->split);
    go_left = Num::Gt(split->first - q.pnt[split->second], 0);
    first_node = go_left ? TI->left : TI->right;
    second_node = go_left ? TI->right : TI->left;
    dim = 0;
    first_idx = second_idx = 1;
  } else if (dynamic_cast<MN*>(T)) {  // NOTE: Multi node
    auto TI = static_cast<MN*>(T);
    assert(TI->GetParallelFlagIniStatus() && TI->aug.value());
    assert(!MN::EqualSplit());
    split = MN::EqualSplit() ? &(TI->split[dim]) : &(TI->split[idx]);
    go_left = Num::Gt(split->first - q.pnt[split->second], 0);

    first_idx = (idx << 1) + static_cast<BucketType>(!go_left);
    second_idx = (idx << 1) + static_cast<BucketType>(go_left);
    bool reach_leaf =
        first_idx >= MN::GetRegions();  // whether reach the skeleton leaf
    first_node = reach_leaf ? TI->tree_nodes[first_idx - MN::GetRegions()] : T;
    second_node =
        reach_leaf ? TI->tree_nodes[second_idx - MN::GetRegions()] : T;
    if (reach_leaf) {
      first_idx = second_idx = 1;
    }
    dim = (dim + 1) % kDim;
  } else {
    assert(false);
  }

  BoxCut box_cut(node_box, *split, go_left);
  logger.generate_box_num += 1;
  KNNMix<Leaf, BN, MN>(first_node, q, dim, first_idx, bq,
                       box_cut.GetFirstBoxCut(), logger);
  logger.check_box_num++;
  if (Num::Gt(P2BMinDistanceSquare(q, box_cut.GetSecondBoxCut()),
              bq.top_value()) &&
      bq.full()) {
    logger.skip_box_num++;
    return;
  }
  KNNMix<Leaf, BN, MN>(second_node, q, dim, second_idx, bq, box_cut.GetBox(),
                       logger);
  return;
}

BASETREE_TEMPLATE
template <typename Leaf, IsDynamicNode Interior, typename Range>
void BASETREE_CLASS::KNNCover(Node* T, Point const& q,
                              kBoundedQueue<Point, Range>& bq,
                              KNNLogger& logger) {
  using CoverCircle = typename Interior::CircleType;

  if (T->size == 0) {
    return;
  }

  if (T->is_leaf) {
    logger.vis_leaf_num++;
    KNNLeaf<Leaf>(T, q, bq);
    return;
  }

  logger.vis_interior_num++;
  Interior* TI = static_cast<Interior*>(T);
  auto sorted_idx = parlay::sort(
      parlay::tabulate(TI->tree_nodes.size(), [&](auto i) { return i; }),
      [&](auto const i1, auto const i2) {
        return Num::Lt(P2PDistanceSquare(q, TI->split[i1]),
                       P2PDistanceSquare(q, TI->split[i2]));
      });

  KNNCover<Leaf, Interior>(TI->tree_nodes[sorted_idx[0]], q, bq, logger);
  for (size_t i = 1; i < TI->tree_nodes.size(); ++i) {
    logger.check_box_num++;
    // auto p2cdis = P2CMinDistance(
    //     q,
    //     CoverCircle{TI->split[sorted_idx[i]], TI->GetCoverCircle().level -
    //     1});
    // auto sqrtbqtop = std::sqrt(bq.top_value());
    if (Num_Comparator<double>::Gt(
            P2CMinDistance(q, CoverCircle{TI->split[sorted_idx[i]],
                                          TI->GetCoverCircle().level - 1}),
            std::sqrt(bq.top_value())) &&
        bq.full()) {
      // if (Num_Comparator<double>::Gt(p2cdis, sqrtbqtop) && bq.full()) {
      logger.skip_box_num++;
      // auto next_node = TI->tree_nodes[sorted_idx[i]];
      // Points p_seq(next_node->size);
      // FlattenRec<Leaf, Interior>(TI->tree_nodes[sorted_idx[i]],
      //                            parlay::make_slice(p_seq));
      // for (auto const& p : p_seq) {
      //   if (P2PDistanceSquare(p, q) <= bq.top_value()) {
      //     std::cout << p << q
      //               << CoverCircle{TI->split[sorted_idx[i]],
      //                              TI->GetCoverCircle().level - 1}
      //               << " ";
      //     std::cout << "p2cdis: " << p2cdis << " ";
      //     std::cout << "sqrtbqtop: " << sqrtbqtop << " ";
      //     std::cout << "p2pdis: " << std::sqrt(P2PDistanceSquare(p, q)) << "
      //     "
      //               << std::flush;
      //     assert(WithinCircle(p, CoverCircle{TI->split[sorted_idx[i]],
      //                                        TI->GetCoverCircle().level -
      //                                        1}));
      //     assert(false);
      //   }
      // }
      continue;
    }
    KNNCover<Leaf, Interior>(TI->tree_nodes[sorted_idx[i]], q, bq, logger);
  }
  return;
}

}  // namespace psi

#undef BASETREE_TEMPLATE
#undef BASETREE_CLASS

#endif  // PSI_BASE_TREE_IMPL_KNN_QUERY_HPP
