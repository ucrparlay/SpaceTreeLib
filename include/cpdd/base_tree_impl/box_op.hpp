#pragma once
#include "../base_tree.h"

namespace cpdd {

template<typename Point, uint_fast8_t kBDO>
inline bool BaseTree<Point, kBDO>::LegalBox(const Box& bx) {
    if (bx == GetEmptyBox()) return true;
    for (DimsType i = 0; i < bx.first.get_dim(); ++i) {
        if (Num::Gt(bx.first.pnt[i], bx.second.pnt[i])) {
            return false;
        }
    }
    return true;
}

template<typename Point, uint_fast8_t kBDO>
inline bool BaseTree<Point, kBDO>::WithinBox(const Box& a, const Box& b) {
    assert(LegalBox(a));
    assert(LegalBox(b));
    for (DimsType i = 0; i < a.first.get_dim(); ++i) {
        if (Num::Lt(a.first.pnt[i], b.first.pnt[i]) ||
            Num::Gt(a.second.pnt[i], b.second.pnt[i])) {
            return false;
        }
    }
    return true;
}

template<typename Point, uint_fast8_t kBDO>
inline bool BaseTree<Point, kBDO>::WithinBox(const Point& p, const Box& bx) {
    assert(LegalBox(bx));
    for (DimsType i = 0; i < p.get_dim(); ++i) {
        if (Num::Lt(p.pnt[i], bx.first.pnt[i]) ||
            Num::Gt(p.pnt[i], bx.second.pnt[i])) {
            return false;
        }
    }
    return true;
}

template<typename Point, uint_fast8_t kBDO>
inline bool BaseTree<Point, kBDO>::BoxIntersectBox(const Box& a, const Box& b) {
    assert(LegalBox(a) && LegalBox(b));
    for (DimsType i = 0; i < a.first.get_dim(); ++i) {
        if (Num::Lt(a.second.pnt[i], b.first.pnt[i]) ||
            Num::Gt(a.first.pnt[i], b.second.pnt[i])) {
            return false;
        }
    }
    return true;
}

template<typename Point, uint_fast8_t kBDO>
inline typename BaseTree<Point, kBDO>::Box
BaseTree<Point, kBDO>::GetEmptyBox() {
    return Box(Point(std::numeric_limits<Coord>::max()),
               Point(std::numeric_limits<Coord>::min()));
}

template<typename Point, uint_fast8_t kBDO>
typename BaseTree<Point, kBDO>::Box BaseTree<Point, kBDO>::GetBox(
    const Box& x, const Box& y) {
    return Box(x.first.minCoords(y.first), x.second.maxCoords(y.second));
}

template<typename Point, uint_fast8_t kBDO>
typename BaseTree<Point, kBDO>::Box BaseTree<Point, kBDO>::GetBox(Slice V) {
    if (V.size() == 0) {
        return GetEmptyBox();
    } else {
        auto minmax = [&](const Box& x, const Box& y) {
            return Box(x.first.minCoords(y.first),
                       x.second.maxCoords(y.second));
        };
        auto boxes = parlay::delayed_seq<Box>(
            V.size(), [&](size_t i) { return Box(V[i].pnt, V[i].pnt); });
        return parlay::reduce(boxes, parlay::make_monoid(minmax, boxes[0]));
    }
}

template<typename Point, uint_fast8_t kBDO>
typename BaseTree<Point, kBDO>::Box BaseTree<Point, kBDO>::GetBox(Node* T) {
    Points wx = Points::uninitialized(T->size);
    Flatten(T, parlay::make_slice(wx));
    return GetBox(parlay::make_slice(wx));
}

template<typename Point, uint_fast8_t kBDO>
inline bool BaseTree<Point, kBDO>::WithinCircle(const Box& bx,
                                                const Circle& cl) {
    //* the logical is same as p2b_max_distance <= radius
    Coord r = 0;
    for (DimsType i = 0; i < cl.first.get_dim(); ++i) {
        if (Num::Lt(cl.first.pnt[i],
                    (bx.first.pnt[i] + bx.second.pnt[i]) / 2)) {
            r += (bx.second.pnt[i] - cl.first.pnt[i]) *
                 (bx.second.pnt[i] - cl.first.pnt[i]);
        } else {
            r += (cl.first.pnt[i] - bx.first.pnt[i]) *
                 (cl.first.pnt[i] - bx.first.pnt[i]);
        }
        if (Num::Gt(r, cl.second * cl.second)) return false;
    }
    assert(Num::Leq(r, cl.second * cl.second));
    return true;
}

template<typename Point, uint_fast8_t kBDO>
inline bool BaseTree<Point, kBDO>::WithinCircle(const Point& p,
                                                const Circle& cl) {
    Coord r = 0;
    for (DimsType i = 0; i < cl.first.get_dim(); ++i) {
        r += (p.pnt[i] - cl.first.pnt[i]) * (p.pnt[i] - cl.first.pnt[i]);
        if (Num::Gt(r, cl.second * cl.second)) return false;
    }
    assert(Num::Leq(r, cl.second * cl.second));
    return true;
}

template<typename Point, uint_fast8_t kBDO>
inline bool BaseTree<Point, kBDO>::CircleIntersectBox(const Circle& cl,
                                                      const Box& bx) {
    //* the logical is same as p2b_min_distance > radius
    Coord r = 0;
    for (DimsType i = 0; i < cl.first.get_dim(); ++i) {
        if (Num::Lt(cl.first.pnt[i], bx.first.pnt[i])) {
            r += (bx.first.pnt[i] - cl.first.pnt[i]) *
                 (bx.first.pnt[i] - cl.first.pnt[i]);
        } else if (Num::Gt(cl.first.pnt[i], bx.second.pnt[i])) {
            r += (cl.first.pnt[i] - bx.second.pnt[i]) *
                 (cl.first.pnt[i] - bx.second.pnt[i]);
        }
        if (Num::Gt(r, cl.second * cl.second)) return false;  //* not intersect
    }
    assert(Num::Leq(r, cl.second * cl.second));
    return true;
}

}  // namespace cpdd
