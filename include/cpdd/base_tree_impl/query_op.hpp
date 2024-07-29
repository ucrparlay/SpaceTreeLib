#pragma once

#include "../base_tree.h"
#include <utility>

namespace cpdd {

// NOTE: NN search
template<typename Point, uint8_t kBDO>
inline typename BaseTree<Point, kBDO>::Coord BaseTree<Point, kBDO>::P2PDistance(
    const Point& p, const Point& q, const DimsType DIM) {
    Coord r = 0;
    for (DimsType i = 0; i < DIM; ++i) {
        r += (p.pnt[i] - q.pnt[i]) * (p.pnt[i] - q.pnt[i]);
    }
    return std::move(r);
}

// NOTE: distance between a Point and a Box
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
template<typename Leaf, typename Interior, typename StoreType>
void BaseTree<Point, kBDO>::KNNRec(Node* T, const Point& q, const DimsType DIM,
                                   kBoundedQueue<Point, StoreType>& bq,
                                   const Box& node_box, size_t& vis_node_num) {
    vis_node_num++;

    if (T->is_leaf) {
        Leaf* TL = static_cast<Leaf*>(T);
        int i = 0;
        while (!bq.full() && i < TL->size) {
            bq.insert(std::make_pair(
                std::ref(TL->pts[(!TL->is_dummy) * i]),
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

    Interior* TI = static_cast<Interior*>(T);
    auto GoLeftTester = [&]() -> bool {
        return Num::Gt(TI->split.first - q.pnt[TI->split.second], 0);
    };

    Box first_box(node_box), second_box(node_box);

    if (GoLeftTester()) {  //* go left child
        first_box.second.pnt[TI->split.second] = TI->split.first;
        second_box.first.pnt[TI->split.second] = TI->split.first;
    } else {  //* go right child
        first_box.first.pnt[TI->split.second] = TI->split.first;
        second_box.second.pnt[TI->split.second] = TI->split.first;
    }

    KNNRec<Leaf, Interior>(GoLeftTester() ? TI->left : TI->right, q, DIM, bq,
                           first_box, vis_node_num);
    if (Num::Gt(P2BMinDistance(q, second_box, DIM), bq.top_value()) &&
        bq.full()) {
        return;
    }
    KNNRec<Leaf, Interior>(GoLeftTester() ? TI->right : TI->left, q, DIM, bq,
                           second_box, vis_node_num);
    return;
}

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
template<typename Leaf, typename Interior, typename StoreType>
void BaseTree<Point, kBDO>::RangeQuerySerialRecursive(Node* T, StoreType Out,
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
