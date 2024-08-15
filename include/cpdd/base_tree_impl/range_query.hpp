#pragma once

#include "../base_tree.h"
#include <utility>

namespace cpdd {
template<typename Point, uint8_t kBDO>
template<typename Leaf>
size_t BaseTree<Point, kBDO>::RangeCountRectangleLeaf(Node* T,
                                                      const Box& query_box,
                                                      const Box& node_box) {
    assert(T->is_leaf);

    Leaf* TL = static_cast<Leaf*>(T);
    size_t cnt = 0;
    if (TL->is_dummy) {
        if (WithinBox(TL->pts[0], query_box)) {
            cnt = TL->size;
        }
    } else {
        for (int i = 0; i < TL->size; i++) {
            if (WithinBox(TL->pts[i], query_box)) {
                cnt++;
            }
        }
    }
    return cnt;
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, IsBinaryNode Interior>
size_t BaseTree<Point, kBDO>::RangeCountRectangle(Node* T, const Box& query_box,
                                                  const Box& node_box) {
    if (T->is_leaf) {
        return RangeCountRectangleLeaf<Leaf>(T, query_box, node_box);
    }

    Interior* TI = static_cast<Interior*>(T);
    Box abox(node_box);

    size_t l, r;
    auto recurse = [&](Node* Ts, const Box& box, size_t& counter) -> void {
        if (!BoxIntersectBox(box, query_box)) {
            counter = 0;
        } else if (WithinBox(box, query_box)) {
            counter = Ts->size;
        } else {
            counter = RangeCountRectangle<Leaf, Interior>(Ts, query_box, box);
        }
    };

    auto& mod_dim = abox.second.pnt[TI->split.second];
    auto split = TI->split.first;
    std::ranges::swap(mod_dim, split);
    recurse(TI->left, abox, l);

    std::ranges::swap(mod_dim, split);
    abox.first.pnt[TI->split.second] = split;
    recurse(TI->right, abox, r);

    return l + r;
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, IsMultiNode Interior>
size_t BaseTree<Point, kBDO>::RangeCountRectangle(Node* T, const Box& query_box,
                                                  const Box& node_box,
                                                  DimsType dim,
                                                  BucketType idx) {
    if (T->is_leaf) {
        return RangeCountRectangleLeaf<Leaf>(T, query_box, node_box);
    }

    Interior* TI = static_cast<Interior*>(T);
    Box abox(node_box);

    auto recurse = [&query_box](Node* Ts, const Box& box, size_t& counter,
                                DimsType next_dim, DimsType next_idx) -> void {
        if (!BoxIntersectBox(box, query_box)) {
            counter = 0;
        } else if (WithinBox(box, query_box)) {
            // NOTE: when reach leaf, Ts is the children
            counter = static_cast<Interior*>(Ts)->ReduceSums(next_idx);
        } else {
            counter = RangeCountRectangle<Leaf, Interior>(Ts, query_box, box,
                                                          next_dim, next_idx);
        }
    };

    size_t l, r;
    idx <<= 1;
    bool reach_leaf = idx >= Interior::kRegions;

    // NOTE: visit left half
    assert(TI->split[dim].second == dim);
    auto& mod_dim = abox.second.pnt[TI->split[dim].second];
    auto split = TI->split[dim].first;
    std::ranges::swap(mod_dim, split);
    recurse(reach_leaf ? TI->tree_nodes[idx - Interior::kRegions] : T, abox, l,
            (dim + 1) % kDim, reach_leaf ? 1 : idx);

    idx |= 1;
    std::ranges::swap(mod_dim, split);
    assert(abox == node_box);
    abox.first.pnt[TI->split[dim].second] = split;
    recurse(reach_leaf ? TI->tree_nodes[idx - Interior::kRegions] : T, abox, r,
            (dim + 1) % kDim, reach_leaf ? 1 : idx);

    return l + r;
}

// TODO: as range_count_rectangle
template<typename Point, uint8_t kBDO>
template<typename Leaf, IsBinaryNode Interior>
size_t BaseTree<Point, kBDO>::RangeCountRadius(Node* T, const Circle& cl,
                                               const Box& node_box) {
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

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Range>
void BaseTree<Point, kBDO>::RangeQueryLeaf(Node* T, Range Out, size_t& s,
                                           const Box& query_box,
                                           const Box& node_box) {
    assert(T->is_leaf);

    Leaf* TL = static_cast<Leaf*>(T);
    if (TL->is_dummy) {
        if (WithinBox(TL->pts[0], query_box)) {
            for (int i = 0; i < TL->size; i++)
                Out[s++] = TL->pts[0];
        }
    } else {
        for (int i = 0; i < TL->size; i++)
            if (WithinBox(TL->pts[i], query_box)) {
                Out[s++] = TL->pts[i];
            }
    }
    return;
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, IsBinaryNode Interior, typename Range>
void BaseTree<Point, kBDO>::RangeQuerySerialRecursive(Node* T, Range Out,
                                                      size_t& s,
                                                      const Box& query_box,
                                                      const Box& node_box) {
    if (T->is_leaf) {
        RangeQueryLeaf<Leaf>(T, Out, s, query_box, node_box);
        return;
    }

    Interior* TI = static_cast<Interior*>(T);
    Box abox(node_box);

    auto recurse = [&](Node* Ts, const Box& box) -> void {
        if (!BoxIntersectBox(box, query_box)) {
            return;
        } else if (WithinBox(box, query_box)) {
            FlattenRec<Leaf, Interior>(Ts, Out.cut(s, s + Ts->size));
            s += Ts->size;
            return;
        } else {
            RangeQuerySerialRecursive<Leaf, Interior>(Ts, Out, s, query_box,
                                                      box);
            return;
        }
    };

    auto& mod_dim = abox.second.pnt[TI->split.second];
    auto split = TI->split.first;
    std::ranges::swap(mod_dim, split);
    recurse(TI->left, abox);

    std::ranges::swap(mod_dim, split);
    abox.first.pnt[TI->split.second] = split;
    recurse(TI->right, abox);

    return;
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, IsMultiNode Interior, typename Range>
void BaseTree<Point, kBDO>::RangeQuerySerialRecursive(
    Node* T, Range Out, size_t& s, const Box& query_box, const Box& node_box,
    DimsType dim, BucketType idx) {
    if (T->is_leaf) {
        RangeQueryLeaf<Leaf>(T, Out, s, query_box, node_box);
        return;
    }

    Interior* TI = static_cast<Interior*>(T);

    auto recurse = [&query_box, &s, &Out](Node* Ts, const Box& box,
                                          DimsType next_dim,
                                          BucketType next_idx) -> void {
        if (!BoxIntersectBox(box, query_box)) {
            return;
        } else if (WithinBox(box, query_box)) {
            size_t candidate_size =
                static_cast<Interior*>(Ts)->ReduceSums(next_idx);
            PartialFlatten<Leaf, Interior>(Ts, Out.cut(s, s + candidate_size),
                                           next_idx);
            s += candidate_size;
            return;
        } else {
            RangeQuerySerialRecursive<Leaf, Interior>(Ts, Out, s, query_box,
                                                      box, next_dim, next_idx);
            return;
        }
    };

    idx <<= 1;
    bool reach_leaf = idx >= Interior::kRegions;
    Box abox(node_box);

    // NOTE: visit left half
    auto& mod_dim = abox.second.pnt[TI->split[dim].second];
    auto split = TI->split[dim].first;
    std::ranges::swap(mod_dim, split);
    recurse(reach_leaf ? TI->tree_nodes[idx - Interior::kRegions] : T, abox,
            (dim + 1) % kDim, reach_leaf ? 1 : idx);

    // NOTE: visit right
    idx |= 1;
    std::ranges::swap(mod_dim, split);
    abox.first.pnt[TI->split[dim].second] = split;
    recurse(reach_leaf ? TI->tree_nodes[idx - Interior::kRegions] : T, abox,
            (dim + 1) % kDim, reach_leaf ? 1 : idx);

    return;
}
}  // namespace cpdd
