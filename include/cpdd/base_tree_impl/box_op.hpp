#pragma once
#include "../base_tree.h"

namespace cpdd {

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
inline bool BaseTree<Point, DerivedTree, kBDO>::LegalBox(Box const& bx) {
  if (bx == GetEmptyBox()) return true;
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Gt(bx.first.pnt[i], bx.second.pnt[i])) {
      return false;
    }
  }
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
inline bool BaseTree<Point, DerivedTree, kBDO>::WithinBox(Box const& a,
                                                          Box const& b) {
  assert(LegalBox(a));
  assert(LegalBox(b));
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(a.first.pnt[i], b.first.pnt[i]) ||
        Num::Gt(a.second.pnt[i], b.second.pnt[i])) {
      return false;
    }
  }
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
inline bool BaseTree<Point, DerivedTree, kBDO>::WithinBox(Point const& p,
                                                          Box const& bx) {
  assert(LegalBox(bx));
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(p.pnt[i], bx.first.pnt[i]) ||
        Num::Gt(p.pnt[i], bx.second.pnt[i])) {
      return false;
    }
  }
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
inline bool BaseTree<Point, DerivedTree, kBDO>::BoxIntersectBox(Box const& a,
                                                                Box const& b) {
  assert(LegalBox(a) && LegalBox(b));
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(a.second.pnt[i], b.first.pnt[i]) ||
        Num::Gt(a.first.pnt[i], b.second.pnt[i])) {
      return false;
    }
  }
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
inline typename BaseTree<Point, DerivedTree, kBDO>::Box
BaseTree<Point, DerivedTree, kBDO>::GetEmptyBox() {
  return Box(Point(std::numeric_limits<Coord>::max()),
             Point(std::numeric_limits<Coord>::lowest()));
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
typename BaseTree<Point, DerivedTree, kBDO>::Box
BaseTree<Point, DerivedTree, kBDO>::GetBox(Box const& x, Box const& y) {
  return Box(x.first.MinCoords(y.first), x.second.MaxCoords(y.second));
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
typename BaseTree<Point, DerivedTree, kBDO>::Box
BaseTree<Point, DerivedTree, kBDO>::GetBox(Slice V) {
  if (V.size() == 0) {
    return GetEmptyBox();
  } else {
    auto minmax = [&](Box const& x, Box const& y) {
      return Box(x.first.MinCoords(y.first), x.second.MaxCoords(y.second));
    };
    auto boxes = parlay::delayed_seq<Box>(
        V.size(), [&](size_t i) { return Box(V[i].pnt, V[i].pnt); });
    return parlay::reduce(boxes, parlay::make_monoid(minmax, boxes[0]));
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
template <typename Leaf, typename Interior>
typename BaseTree<Point, DerivedTree, kBDO>::Box
BaseTree<Point, DerivedTree, kBDO>::GetBox(Node* T) {
  // TODO: we can just traverse the tree to get the box
  Points wx = Points::uninitialized(T->size);
  FlattenRec<Leaf, Interior>(T, parlay::make_slice(wx));
  return GetBox(parlay::make_slice(wx));
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
typename BaseTree<Point, DerivedTree, kBDO>::Box
BaseTree<Point, DerivedTree, kBDO>::GetBox(BoxSeq const& box_seq) {
  Box box = GetEmptyBox();
  for (auto const& b : box_seq) {
    box = GetBox(box, b);
  }
  return std::move(box);
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
inline bool BaseTree<Point, DerivedTree, kBDO>::WithinCircle(Box const& bx,
                                                             Circle const& cl) {
  //* the logical is same as p2b_max_distance <= radius
  Coord r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(cl.first.pnt[i], (bx.first.pnt[i] + bx.second.pnt[i]) / 2)) {
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

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
inline bool BaseTree<Point, DerivedTree, kBDO>::WithinCircle(Point const& p,
                                                             Circle const& cl) {
  Coord r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
    r += (p.pnt[i] - cl.first.pnt[i]) * (p.pnt[i] - cl.first.pnt[i]);
    if (Num::Gt(r, cl.second * cl.second)) return false;
  }
  assert(Num::Leq(r, cl.second * cl.second));
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
inline bool BaseTree<Point, DerivedTree, kBDO>::CircleIntersectBox(
    Circle const& cl, Box const& bx) {
  //* the logical is same as p2b_min_distance > radius
  Coord r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
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
