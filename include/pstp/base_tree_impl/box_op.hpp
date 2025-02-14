#ifndef PSTP_BASE_TREE_IMPL_BOX_OP_HPP_
#define PSTP_BASE_TREE_IMPL_BOX_OP_HPP_

#include "../base_tree.h"
#include "pstp/dependence/concepts.h"

namespace pstp {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Coord
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetBoxMid(DimsType const d,
                                                               Box const& box) {
  if constexpr (std::is_integral_v<Coord>) {
    // NOTE: since the points on the box line will be put in right, therefore
    // the box should always be rounded up as well.
    // Consideing the example (1, 1) and (1, 2), in order to split the points,
    // the new box should have split at y=2
    return static_cast<Coord>(
        std::ceil(static_cast<double>(box.first.pnt[d] + box.second.pnt[d]) /
                  static_cast<double>(2)));
  } else {
    return (box.first.pnt[d] + box.second.pnt[d]) / 2;
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::LegalBox(
    Box const& bx) {
  // TODO: remove it
  if (bx == GetEmptyBox()) return true;

  if constexpr (kDim == 2) {
    return !Num::Gt(bx.first.pnt[0], bx.second.pnt[0]) &&
           !Num::Gt(bx.first.pnt[1], bx.second.pnt[1]);
  } else if constexpr (kDim == 3) {
    return !Num::Gt(bx.first.pnt[0], bx.second.pnt[0]) &&
           !Num::Gt(bx.first.pnt[1], bx.second.pnt[1]) &&
           !Num::Gt(bx.first.pnt[2], bx.second.pnt[2]);
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      if (Num::Gt(bx.first.pnt[i], bx.second.pnt[i])) {
        return false;
      }
    }
    return true;
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::WithinBox(
    Box const& a, Box const& b) {
  constexpr uint_fast8_t kDim = Point::GetDim();
  assert(LegalBox(a));
  assert(LegalBox(b));

  if constexpr (kDim == 2) {
    return !Num::Lt(a.first.pnt[0], b.first.pnt[0]) &&
           !Num::Lt(a.first.pnt[1], b.first.pnt[1]) &&
           !Num::Gt(a.second.pnt[0], b.second.pnt[0]) &&
           !Num::Gt(a.second.pnt[1], b.second.pnt[1]);
  } else if constexpr (kDim == 3) {
    return !Num::Lt(a.first.pnt[0], b.first.pnt[0]) &&
           !Num::Lt(a.first.pnt[1], b.first.pnt[1]) &&
           !Num::Lt(a.first.pnt[2], b.first.pnt[2]) &&
           !Num::Gt(a.second.pnt[0], b.second.pnt[0]) &&
           !Num::Gt(a.second.pnt[1], b.second.pnt[1]) &&
           !Num::Gt(a.second.pnt[2], b.second.pnt[2]);
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      if (Num::Lt(a.first.pnt[i], b.first.pnt[i]) ||
          Num::Gt(a.second.pnt[i], b.second.pnt[i])) {
        return false;
      }
    }
    return true;
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::SameBox(
    Box const& a, Box const& b) {
  assert(LegalBox(a));
  assert(LegalBox(b));

  if constexpr (kDim == 2) {
    return Num::Eq(a.first.pnt[0], b.first.pnt[0]) &&
           Num::Eq(a.first.pnt[1], b.first.pnt[1]) &&
           Num::Eq(a.second.pnt[0], b.second.pnt[0]) &&
           Num::Eq(a.second.pnt[1], b.second.pnt[1]);
  } else if constexpr (kDim == 3) {
    return Num::Eq(a.first.pnt[0], b.first.pnt[0]) &&
           Num::Eq(a.first.pnt[1], b.first.pnt[1]) &&
           Num::Eq(a.first.pnt[2], b.first.pnt[2]) &&
           Num::Eq(a.second.pnt[0], b.second.pnt[0]) &&
           Num::Eq(a.second.pnt[1], b.second.pnt[1]) &&
           Num::Eq(a.second.pnt[2], b.second.pnt[2]);
    for (DimsType i = 0; i < kDim; ++i) {
      if (!Num::Eq(a.first.pnt[i], b.first.pnt[i]) ||
          !Num::Eq(a.second.pnt[i], b.second.pnt[i])) {
        return false;
      }
    }
    return true;
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::WithinBox(
    Point const& p, Box const& bx) {
  assert(LegalBox(bx));

  if constexpr (kDim == 2) {
    return !Num::Lt(p.pnt[0], bx.first.pnt[0]) &&
           !Num::Lt(p.pnt[1], bx.first.pnt[1]) &&
           !Num::Gt(p.pnt[0], bx.second.pnt[0]) &&
           !Num::Gt(p.pnt[1], bx.second.pnt[1]);
  } else if constexpr (kDim == 3) {
    return !Num::Lt(p.pnt[0], bx.first.pnt[0]) &&
           !Num::Lt(p.pnt[1], bx.first.pnt[1]) &&
           !Num::Lt(p.pnt[2], bx.first.pnt[2]) &&
           !Num::Gt(p.pnt[0], bx.second.pnt[0]) &&
           !Num::Gt(p.pnt[1], bx.second.pnt[1]) &&
           !Num::Gt(p.pnt[2], bx.second.pnt[2]);
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      if (Num::Lt(p.pnt[i], bx.first.pnt[i]) ||
          Num::Gt(p.pnt[i], bx.second.pnt[i])) {
        return false;
      }
    }
    return true;
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::BoxIntersectBox(Box const& a, Box const& b) {
  assert(LegalBox(a) && LegalBox(b));

  if constexpr (kDim == 2) {
    return !(Num::Lt(a.second.pnt[0], b.first.pnt[0]) ||
             Num::Gt(a.first.pnt[0], b.second.pnt[0]) ||
             Num::Lt(a.second.pnt[1], b.first.pnt[1]) ||
             Num::Gt(a.first.pnt[1], b.second.pnt[1]));
  } else if constexpr (kDim == 3) {
    return !(Num::Lt(a.second.pnt[0], b.first.pnt[0]) ||
             Num::Gt(a.first.pnt[0], b.second.pnt[0]) ||
             Num::Lt(a.second.pnt[1], b.first.pnt[1]) ||
             Num::Gt(a.first.pnt[1], b.second.pnt[1]) ||
             Num::Lt(a.second.pnt[2], b.first.pnt[2]) ||
             Num::Gt(a.first.pnt[2], b.second.pnt[2]));
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      if (Num::Lt(a.second.pnt[i], b.first.pnt[i]) ||
          Num::Gt(a.first.pnt[i], b.second.pnt[i])) {
        return false;
      }
    }
    return true;
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::IsBoxLineInDimension(Box const& box,
                                                       DimsType d) {
  assert(LegalBox(box));
  return Num::Eq(box.first[d], box.second[d]);
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::VerticalLineSplitBox(Coord const& line,
                                                       Box const& box,
                                                       DimsType d) {
  assert(LegalBox(box));
  return VerticalLineIntersectBoxExclude(line, box, d) ||
         (VerticalLineOnBoxRightEdge(line, box, d) &&
          !IsBoxLineInDimension(box, d));
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::VerticalLineOnBoxLeftEdge(Coord const& line,
                                                            Box const& box,
                                                            DimsType d) {
  assert(LegalBox(box));
  return Num::Eq(line, box.first.pnt[d]);
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::VerticalLineOnBoxRightEdge(Coord const& line,
                                                             Box const& box,
                                                             DimsType d) {
  assert(LegalBox(box));
  return Num::Eq(line, box.second.pnt[d]);
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::VerticalLineOnBoxEdge(Coord const& line,
                                                        Box const& box,
                                                        DimsType d) {
  assert(LegalBox(box));
  return VerticalLineOnBoxLeftEdge(line, box, d) ||
         VerticalLineOnBoxRightEdge(line, box, d);
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::VerticalLineIntersectBox(Coord const& line,
                                                           Box const& box,
                                                           DimsType d) {
  assert(LegalBox(box));
  return Num::Geq(line, box.first.pnt[d]) && Num::Leq(line, box.second.pnt[d]);
}

// NOTE: if the line @line is one the boundary of the box, then it will be
// considered as not intersect
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::
    VerticalLineIntersectBoxExclude(Coord const& line, Box const& box,
                                    DimsType d) {
  assert(LegalBox(box));
  return Num::Gt(line, box.first.pnt[d]) && Num::Lt(line, box.second.pnt[d]);
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
Point BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetBoxCenter(
    Box const& box) {
  Point center;
  if constexpr (kDim == 2) {
    center.pnt[0] = (box.first.pnt[0] + box.second.pnt[0]) / 2;
    center.pnt[1] = (box.first.pnt[1] + box.second.pnt[1]) / 2;
  } else if constexpr (kDim == 3) {
    center.pnt[0] = (box.first.pnt[0] + box.second.pnt[0]) / 2;
    center.pnt[1] = (box.first.pnt[1] + box.second.pnt[1]) / 2;
    center.pnt[2] = (box.first.pnt[2] + box.second.pnt[2]) / 2;
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      center.pnt[i] = GetBoxMid(i, box);
    }
  }
  return center;
}
}  // namespace pstp

#endif  // PSTP_BASE_TREE_IMPL_BOX_OP_HPP_
