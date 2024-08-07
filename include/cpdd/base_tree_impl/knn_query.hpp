#pragma once
#include <algorithm>
#include <tuple>
#include "../base_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {

// NOTE: distance between two Points
template<typename Point, uint8_t kBDO>
inline typename BaseTree<Point, kBDO>::Coord BaseTree<Point, kBDO>::P2PDistance(
    const Point& p, const Point& q, const DimsType DIM) {
    Coord r = 0;
    for (DimsType i = 0; i < DIM; ++i) {
        r += (p.pnt[i] - q.pnt[i]) * (p.pnt[i] - q.pnt[i]);
    }
    return r;
}

// NOTE: Distance between a Point and a Box
// return 0 when p is inside the box a
template<typename Point, uint8_t kBDO>
inline typename BaseTree<Point, kBDO>::Coord
BaseTree<Point, kBDO>::P2BMinDistance(
    const Point& p, const typename BaseTree<Point, kBDO>::Box& a,
    const DimsType DIM) {
    Coord r = 0;
    for (DimsType i = 0; i < DIM; ++i) {
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
    const Point& p, const typename BaseTree<Point, kBDO>::Box& a,
    const DimsType DIM) {
    Coord r = 0;
    for (DimsType i = 0; i < DIM; ++i) {
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
                                             Coord up, DimsType DIM) {
    Coord r = 0;
    DimsType i = 0;
    if (DIM >= 6) {
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
            if (i + 4 > DIM) {
                break;
            }
        }
    }
    while (i < DIM) {
        r += (p.pnt[i] - q.pnt[i]) * (p.pnt[i] - q.pnt[i]);
        ++i;
    }
    return r;
}

// NOTE: KNN search for Point q
template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Range>
void BaseTree<Point, kBDO>::KNNLeaf(Node* T, const Point& q, const DimsType DIM,
                                    kBoundedQueue<Point, Range>& bq,
                                    const Box& node_box, size_t& vis_node_num) {
    assert(T->is_leaf);
    vis_node_num++;

    Leaf* TL = static_cast<Leaf*>(T);
    int i = 0;
    while (!bq.full() && i < TL->size) {
        bq.insert(
            std::make_pair(std::ref(TL->pts[(!TL->is_dummy) * i]),
                           P2PDistance(q, TL->pts[(!TL->is_dummy) * i], DIM)));
        i++;
    }
    while (i < TL->size) {
        Coord r = InterruptibleDistance(q, TL->pts[(!TL->is_dummy) * i],
                                        bq.top_value(), DIM);
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
                                      const DimsType DIM,
                                      kBoundedQueue<Point, Range>& bq,
                                      const Box& node_box,
                                      size_t& vis_node_num) {
    vis_node_num++;

    if (T->is_leaf) {
        KNNLeaf<Leaf>(T, q, DIM, bq, node_box, vis_node_num);
        return;
    }

    Interior* TI = static_cast<Interior*>(T);
    bool go_left = Num::Gt(TI->split.first - q.pnt[TI->split.second], 0);

    Box first_box(node_box), second_box(node_box);

    if (go_left) {  // NOTE: go left child
        first_box.second.pnt[TI->split.second] = TI->split.first;
        second_box.first.pnt[TI->split.second] = TI->split.first;
    } else {  // NOTE: go right child
        first_box.first.pnt[TI->split.second] = TI->split.first;
        second_box.second.pnt[TI->split.second] = TI->split.first;
    }

    KNNBinary<Leaf, Interior>(go_left ? TI->left : TI->right, q, DIM, bq,
                              first_box, vis_node_num);
    if (Num::Gt(P2BMinDistance(q, second_box, DIM), bq.top_value()) &&
        bq.full()) {
        return;
    }
    KNNBinary<Leaf, Interior>(go_left ? TI->right : TI->left, q, DIM, bq,
                              second_box, vis_node_num);
    return;
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, IsMultiNode Interior, typename Range>
void BaseTree<Point, kBDO>::KNNMulti(Node* T, const Point& q,
                                     const DimsType DIM,
                                     kBoundedQueue<Point, Range>& bq,
                                     const Box& node_box,
                                     size_t& vis_node_num) {
    // TODO: using MN::kNodeRegions
    constexpr size_t const kNodeRegions =
        std::tuple_size_v<typename Interior::Nodes>;

    vis_node_num++;

    if (T->size == 0) {
        return;
    }

    if (T->is_leaf) {
        KNNLeaf<Leaf>(T, q, DIM, bq, node_box, vis_node_num);
        return;
    }

    Interior* TI = static_cast<Interior*>(T);

    BoxSeq regions(kNodeRegions);
    std::array<std::pair<Coord, BucketType>, kNodeRegions> dists;
    TI->ComputeSubregions(regions, node_box, 1, 0);

    std::ranges::generate(dists, [i = 0, &q, &regions, &DIM]() mutable {
        auto r = std::make_pair(P2BMinDistance(q, regions[i], DIM), i);
        i++;
        return r;
    });
    std::ranges::sort(dists, std::less<>(),
                      [&](const auto& box_pair) { return box_pair.first; });

    KNNMulti<Leaf, Interior>(TI->tree_nodes[dists[0].second], q, DIM, bq,
                             regions[dists[0].second], vis_node_num);
    for (BucketType i = 1; i < kNodeRegions; ++i) {
        if (Num::Gt(dists[i].first, bq.top_value()) && bq.full()) {
            continue;
        }
        KNNMulti<Leaf, Interior>(TI->tree_nodes[dists[i].second], q, DIM, bq,
                                 regions[dists[i].second], vis_node_num);
    }

    return;
}

}  // namespace cpdd
