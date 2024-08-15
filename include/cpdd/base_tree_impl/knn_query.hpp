#pragma once
#include <algorithm>
#include <tuple>
#include <utility>
#include "../base_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {

// NOTE: distance between two Points
template<typename Point, uint8_t kBDO>
inline typename BaseTree<Point, kBDO>::Coord BaseTree<Point, kBDO>::P2PDistance(
    const Point& p, const Point& q) {
    Coord r = 0;
    for (DimsType i = 0; i < kDim; ++i) {
        r += (p.pnt[i] - q.pnt[i]) * (p.pnt[i] - q.pnt[i]);
    }
    return r;
}

// NOTE: Distance between a Point and a Box
// return 0 when p is inside the box a
template<typename Point, uint8_t kBDO>
inline typename BaseTree<Point, kBDO>::Coord
BaseTree<Point, kBDO>::P2BMinDistance(
    const Point& p, const typename BaseTree<Point, kBDO>::Box& a) {
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
template<typename Point, uint8_t kBDO>
inline typename BaseTree<Point, kBDO>::Coord
BaseTree<Point, kBDO>::P2BMaxDistance(
    const Point& p, const typename BaseTree<Point, kBDO>::Box& a) {
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
template<typename Point, uint8_t kBDO>
inline typename BaseTree<Point, kBDO>::Coord
BaseTree<Point, kBDO>::InterruptibleDistance(const Point& p, const Point& q,
                                             Coord up) {
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
template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Range>
void BaseTree<Point, kBDO>::KNNLeaf(Node* T, const Point& q,
                                    kBoundedQueue<Point, Range>& bq,
                                    const Box& node_box) {
    assert(T->is_leaf);

    Leaf* TL = static_cast<Leaf*>(T);
    int i = 0;
    while (!bq.full() && i < TL->size) {
        bq.insert(std::make_pair(std::ref(TL->pts[(!TL->is_dummy) * i]),
                                 P2PDistance(q, TL->pts[(!TL->is_dummy) * i])));
        i++;
    }
    while (i < TL->size) {
        Coord r = InterruptibleDistance(q, TL->pts[(!TL->is_dummy) * i],
                                        bq.top_value());
        if (Num::Lt(r, bq.top_value())) {
            bq.insert(
                std::make_pair(std::ref(TL->pts[(!TL->is_dummy) * i]), r));
        } else if (TL->is_dummy) {
            break;
        }
        i++;
    }
    return;
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, IsBinaryNode Interior, typename Range>
void BaseTree<Point, kBDO>::KNNBinary(Node* T, const Point& q,
                                      kBoundedQueue<Point, Range>& bq,
                                      const Box& node_box, KNNLogger& logger) {
    logger.vis_node_num++;

    if (T->is_leaf) {
        KNNLeaf<Leaf>(T, q, bq, node_box);
        return;
    }

    Interior* TI = static_cast<Interior*>(T);
    bool go_left = Num::Gt(TI->split.first - q.pnt[TI->split.second], 0);

    // Box first_box(node_box), second_box(node_box);
    Box next_box(node_box);
    Coord* mod_dim = go_left ? &next_box.second.pnt[TI->split.second]
                             : &next_box.first.pnt[TI->split.second];
    auto split = TI->split.first;
    std::ranges::swap(split, *mod_dim);
    logger.generate_box_num += 1;

    KNNBinary<Leaf, Interior>(go_left ? TI->left : TI->right, q, bq, next_box,
                              logger);

    logger.check_box_num++;
    std::ranges::swap(split, *mod_dim);
    mod_dim = go_left ? &next_box.first.pnt[TI->split.second]
                      : &next_box.second.pnt[TI->split.second];
    *mod_dim = split;
    if (Num::Gt(P2BMinDistance(q, next_box), bq.top_value()) && bq.full()) {
        logger.skip_box_num++;
        return;
    }
    KNNBinary<Leaf, Interior>(go_left ? TI->right : TI->left, q, bq, next_box,
                              logger);
    return;
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, IsMultiNode Interior, typename Range>
void BaseTree<Point, kBDO>::KNNMultiExpand(Node* T, const Point& q,
                                           DimsType dim, BucketType idx,
                                           kBoundedQueue<Point, Range>& bq,
                                           const Box& node_box,
                                           KNNLogger& logger) {
    logger.vis_node_num++;

    if (T->size == 0) {
        return;
    }

    if (T->is_leaf) {
        KNNLeaf<Leaf>(T, q, bq, node_box);
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

    // Box first_box(node_box), second_box(node_box);
    Box next_box(node_box);
    logger.generate_box_num += 1;
    Coord* mod_dim = go_left ? &next_box.second.pnt[TI->split[dim].second]
                             : &next_box.first.pnt[TI->split[dim].second];
    auto split = TI->split[dim].first;
    std::ranges::swap(split, *mod_dim);
    // first_box.second.pnt[TI->split[dim].second] = TI->split[dim].first;
    // second_box.first.pnt[TI->split[dim].second] = TI->split[dim].first;
    // if (!go_left) {
    //     std::ranges::swap(first_box, second_box);
    // }

    assert(dim != 0 || (first_idx == 1 && second_idx == 1));

    KNNMultiExpand<Leaf, Interior>(first_node, q, (dim + 1) % kDim, first_idx,
                                   bq, next_box, logger);

    // NOTE: compute the other bounding box
    logger.check_box_num++;
    std::ranges::swap(split, *mod_dim);
    mod_dim = go_left ? &next_box.first.pnt[TI->split[dim].second]
                      : &next_box.second.pnt[TI->split[dim].second];
    *mod_dim = split;
    if (Num::Gt(P2BMinDistance(q, next_box), bq.top_value()) && bq.full()) {
        logger.skip_box_num++;
        return;
    }
    KNNMultiExpand<Leaf, Interior>(second_node, q, (dim + 1) % kDim, second_idx,
                                   bq, next_box, logger);
    return;
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, IsMultiNode Interior, typename Range>
void BaseTree<Point, kBDO>::KNNMulti(Node* T, const Point& q,
                                     kBoundedQueue<Point, Range>& bq,
                                     const Box& node_box, KNNLogger& logger) {
    logger.vis_node_num++;

    if (T->size == 0) {
        return;
    }

    if (T->is_leaf) {
        KNNLeaf<Leaf>(T, q, bq, node_box);
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
                      [&](const auto& box_pair) { return box_pair.first; });

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
