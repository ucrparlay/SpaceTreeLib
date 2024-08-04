#pragma once

#include "../base_tree.h"
#include <utility>

namespace cpdd {

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior>
size_t BaseTree<Point, kBDO>::RangeCountRectangle(Node* T, const Box& query_box,
                                                  const Box& node_box,
                                                  size_t& vis_leaf_num,
                                                  size_t& vis_inter_num) {
    if (T->is_leaf) {
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

    Interior* TI = static_cast<Interior*>(T);
    // Box lbox( node_box ), rbox( node_box );
    Box abox(node_box);
    // lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
    // rbox.first.pnt[TI->split.second] = TI->split.first;

    size_t l, r;
    auto recurse = [&](Node* Ts, const Box& bx, size_t& counter) -> void {
        if (!BoxIntersectBox(bx, query_box)) {
            counter = 0;
        } else if (WithinBox(bx, query_box)) {
            counter = Ts->size;
        } else {
            counter = RangeCountRectangle<Leaf, Interior>(
                Ts, query_box, bx, vis_leaf_num, vis_inter_num);
        }
    };

    auto& mod_dim = abox.second.pnt[TI->split.second];
    auto split = TI->split.first;
    std::swap(mod_dim, split);
    recurse(TI->left, abox, l);

    std::swap(mod_dim, split);
    abox.first.pnt[TI->split.second] = split;
    recurse(TI->right, abox, r);

    return l + r;
}

// TODO as range_count_rectangle
template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior>
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
        return std::move(cnt);
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

    return std::move(l + r);
};

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior, typename Range>
void BaseTree<Point, kBDO>::RangeQuerySerialRecursive(Node* T, Range Out,
                                                      size_t& s,
                                                      const Box& query_box,
                                                      const Box& node_box) {
    if (T->is_leaf) {
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

    Interior* TI = static_cast<Interior*>(T);
    // Box lbox( node_box ), rbox( node_box );
    Box abox(node_box);
    // lbox.second.pnt[TI->split.second] = TI->split.first;  //* loose
    // rbox.first.pnt[TI->split.second] = TI->split.first;

    auto recurse = [&](Node* Ts, const Box& bx) -> void {
        if (!BoxIntersectBox(bx, query_box)) {
            return;
        } else if (WithinBox(bx, query_box)) {
            FlattenRec<Leaf, Interior>(Ts, Out.cut(s, s + Ts->size));
            s += Ts->size;
            return;
        } else {
            RangeQuerySerialRecursive<Leaf, Interior>(Ts, Out, s, query_box,
                                                      bx);
            return;
        }
    };

    auto& mod_dim = abox.second.pnt[TI->split.second];
    auto split = TI->split.first;
    std::swap(mod_dim, split);
    recurse(TI->left, abox);

    std::swap(mod_dim, split);
    abox.first.pnt[TI->split.second] = split;
    recurse(TI->right, abox);

    return;
}
}  // namespace cpdd
