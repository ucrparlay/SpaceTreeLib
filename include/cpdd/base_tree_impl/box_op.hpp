#pragma once
#include "../base_tree.h"
#include "cpdd/dependence/concepts.h"

namespace cpdd {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Coord
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetBoxMid(DimsType const d,
                                                               Box const& box) {
  if constexpr (std::is_integral_v<Coord>) {
    // NOTE: since the points on the box line, will be put in right, therefore
    // the box should always be rounded up as well. Consideing the example (1,
    // 1) and (1, 2), in order to split the points, the new box should have
    // split at y=2
    return static_cast<Coord>(
        std::ceil(static_cast<double>(box.first.pnt[d] + box.second.pnt[d]) /
                  static_cast<double>(2)));
  } else if constexpr (std::is_floating_point_v<Coord>) {
    return (box.first.pnt[d] + box.second.pnt[d]) / 2;
  } else {
    return (box.first.pnt[d] + box.second.pnt[d]) / 2;
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::LegalBox(
    Box const& bx) {
  if (bx == GetEmptyBox()) return true;
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Gt(bx.first.pnt[i], bx.second.pnt[i])) {
      return false;
    }
  }
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::WithinBox(
    Box const& a, Box const& b) {
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

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::SameBox(
    Box const& a, Box const& b) {
  assert(LegalBox(a));
  assert(LegalBox(b));
  for (DimsType i = 0; i < kDim; ++i) {
    if (!Num::Eq(a.first.pnt[i], b.first.pnt[i]) ||
        !Num::Eq(a.second.pnt[i], b.second.pnt[i])) {
      return false;
    }
  }
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::WithinBox(
    Point const& p, Box const& bx) {
  assert(LegalBox(bx));
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(p.pnt[i], bx.first.pnt[i]) ||
        Num::Gt(p.pnt[i], bx.second.pnt[i])) {
      return false;
    }
  }
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::BoxIntersectBox(Box const& a, Box const& b) {
  assert(LegalBox(a) && LegalBox(b));
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(a.second.pnt[i], b.first.pnt[i]) ||
        Num::Gt(a.first.pnt[i], b.second.pnt[i])) {
      return false;
    }
  }
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetEmptyBox() {
  return Box(Point(std::numeric_limits<Coord>::max()),
             Point(std::numeric_limits<Coord>::lowest()));
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetBox(Box const& x,
                                                            Box const& y) {
  return Box(x.first.MinCoords(y.first), x.second.MaxCoords(y.second));
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetBox(Slice V) {
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

// NOTE: this function omit the possibility that T contains the bounding box --
// it will always try to reduce a bounding box in the tree T
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetBox(Node* T) {
  if (T->is_leaf) {
    Leaf* TL = static_cast<Leaf*>(T);
    if (TL->is_dummy) {
      return Box(TL->pts[0], TL->pts[0]);
    }
    return GetBox(TL->pts.cut(0, TL->size));
  }
  Interior* TI = static_cast<Interior*>(T);
  if constexpr (IsBinaryNode<Interior>) {
    Box const& left_box = GetBox<Leaf, Interior>(TI->left);
    Box const& right_box = GetBox<Leaf, Interior>(TI->right);
    return GetBox(left_box, right_box);
  } else if constexpr (IsMultiNode<Interior>) {
    BoxSeq return_box_seq(Interior::GetRegions());
    for (size_t i = 0; i < Interior::GetRegions(); i++) {
      return_box_seq[i] = GetBox<Leaf, Interior>(TI->tree_nodes[i]);
    }
    return GetBox(return_box_seq);
  } else {
    assert(false);
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetBox(
    BoxSeq const& box_seq) {
  Box box = GetEmptyBox();
  for (auto const& b : box_seq) {
    box = GetBox(box, b);
  }
  return std::move(box);
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::WithinCircle(
    Box const& bx, Circle const& cl) {
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

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::WithinCircle(
    Point const& p, Circle const& cl) {
  Coord r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
    r += (p.pnt[i] - cl.first.pnt[i]) * (p.pnt[i] - cl.first.pnt[i]);
    if (Num::Gt(r, cl.second * cl.second)) return false;
  }
  assert(Num::Leq(r, cl.second * cl.second));
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::CircleIntersectBox(Circle const& cl,
                                                     Box const& bx) {
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
