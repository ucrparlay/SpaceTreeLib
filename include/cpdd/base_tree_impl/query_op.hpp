#pragma once

#include "../base_tree.h"
#include <utility>

namespace cpdd {

//* NN search
template<typename Point>
inline typename BaseTree<Point>::Coord BaseTree<Point>::P2PDistance(
    const Point& p, const Point& q, const DimsType DIM) {
    Coord r = 0;
    for (DimsType i = 0; i < DIM; ++i) {
        r += (p.pnt[i] - q.pnt[i]) * (p.pnt[i] - q.pnt[i]);
    }
    return std::move(r);
}

//* distance between a Point and a Box
template<typename Point>
inline typename BaseTree<Point>::Coord BaseTree<Point>::P2BMinDistance(
    const Point& p, const typename BaseTree<Point>::Box& a,
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

template<typename Point>
inline typename BaseTree<Point>::Coord BaseTree<Point>::P2BMaxDistance(
    const Point& p, const typename BaseTree<Point>::Box& a,
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

//* early return the partial distance between p and q if it is larger than r
//* else return the distance between p and q
template<typename Point>
inline typename BaseTree<Point>::Coord BaseTree<Point>::InterruptibleDistance(
    const Point& p, const Point& q, Coord up, DimsType DIM) {
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

template<typename Point>
template<typename Leaf, typename Interior, typename StoreType>
void BaseTree<Point>::KNNRec(Node* T, const Point& q, const DimsType DIM,
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
}  // namespace cpdd
