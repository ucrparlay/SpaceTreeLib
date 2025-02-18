#ifndef PSPT_BASE_TREE_IMPL_CIRCLE_OP_HPP_
#define PSPT_BASE_TREE_IMPL_CIRCLE_OP_HPP_

#include "../base_tree.h"
#include "pspt/dependence/concepts.h"

namespace pspt {
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename CircleType>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::LegalCircle(
    CircleType const& cl) {
  return Num::Geq(cl.GetRadius(), 0);
}

// NOTE: point within the circle
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename CircleType>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::WithinCircle(
    Point const& p, CircleType const& cl) {
  Coord r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
    r +=
        (p.pnt[i] - cl.GetCenter().pnt[i]) * (p.pnt[i] - cl.GetCenter().pnt[i]);
    if (Num::Gt(r, cl.GetRadius() * cl.GetRadius())) return false;
  }
  assert(Num::Leq(r, cl.GetRadius() * cl.GetRadius()));
  return true;
}

// NOTE: box within the circle
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename CircleType>
inline bool BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::WithinCircle(
    Box const& box, CircleType const& cl) {
  assert(LegalBox(box));
  assert(LegalCircle(cl));

  // NOTE: the logical is same as p2b_max_distance <= radius
  Coord r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(cl.GetCenter().pnt[i],
                (box.first.pnt[i] + box.second.pnt[i]) / 2)) {
      r += (box.second.pnt[i] - cl.GetCenter().pnt[i]) *
           (box.second.pnt[i] - cl.GetCenter().pnt[i]);
    } else {
      r += (cl.GetCenter().pnt[i] - box.first.pnt[i]) *
           (cl.GetCenter().pnt[i] - box.first.pnt[i]);
    }
    if (Num::Gt(r, cl.GetRadius() * cl.GetRadius())) return false;
  }
  assert(Num::Leq(r, cl.GetRadius() * cl.GetRadius()));
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename CircleType>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::CircleIntersectBox(CircleType const& cl,
                                                     Box const& box) {
  assert(LegalBox(box));
  assert(LegalCircle(cl));
  // NOTE: the logical is same as p2b_min_distance > radius
  Coord r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(cl.GetCenter().pnt[i], box.first.pnt[i])) {
      r += (box.first.pnt[i] - cl.GetCenter().pnt[i]) *
           (box.first.pnt[i] - cl.GetCenter().pnt[i]);
    } else if (Num::Gt(cl.GetCenter().pnt[i], box.second.pnt[i])) {
      r += (cl.GetCenter().pnt[i] - box.second.pnt[i]) *
           (cl.GetCenter().pnt[i] - box.second.pnt[i]);
    }
    if (Num::Gt(r, cl.GetRadius() * cl.GetRadius()))
      return false;  //* not intersect
  }
  assert(Num::Leq(r, cl.GetRadius() * cl.GetRadius()));
  return true;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename CircleType>
inline CircleType
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetCircle(Box const& box) {
  assert(LegalBox(box) && box != GetEmptyBox());

  Coord r = 0;
  if constexpr (kDim == 2) {
    r = (box.second.pnt[0] - box.first.pnt[0]) *
            (box.second.pnt[0] - box.first.pnt[0]) +
        (box.second.pnt[1] - box.first.pnt[1]) *
            (box.second.pnt[1] - box.first.pnt[1]);
  } else if constexpr (kDim == 3) {
    r = (box.second.pnt[0] - box.first.pnt[0]) *
            (box.second.pnt[0] - box.first.pnt[0]) +
        (box.second.pnt[1] - box.first.pnt[1]) *
            (box.second.pnt[1] - box.first.pnt[1]) +
        (box.second.pnt[2] - box.first.pnt[2]) *
            (box.second.pnt[2] - box.first.pnt[2]);
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      r += (box.second.pnt[i] - box.first.pnt[i]) *
           (box.second.pnt[i] - box.first.pnt[i]);
    }
  }

  return CircleType{GetBoxCenter(box),
                    CircleType::ComputeRadius(std::sqrt(r) / 2.0)};
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename CircleType>
inline CircleType
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetCircle(Slice V) {
  return GetCircle<CircleType>(GetBox(V));
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename CircleType, typename CircleSeq>
inline CircleType BaseTree<Point, DerivedTree, kSkHeight,
                           kImbaRatio>::GetCircle(CircleSeq const& circle_seq) {
  CircleType circle = circle_seq[0];
  for (size_t i = 1; i < circle_seq.size(); i++) {
    circle = GetCircle(circle, circle_seq[i]);
  }
  return circle;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename CircleType>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::CircleWithinCircle(CircleType const& a,
                                                     CircleType const& b) {
  return Num::Leq(a.GetRadius(), b.GetRadius()) &&
         Num::Leq(P2PDistance(a.GetCenter(), b.GetCenter()),
                  b.GetRadius() - a.GetRadius());
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename CircleType>
inline CircleType BaseTree<Point, DerivedTree, kSkHeight,
                           kImbaRatio>::GetCircle(CircleType const& a,
                                                  CircleType const& b) {
  if (CircleWithinCircle(a, b)) {
    return b;
  } else if (CircleWithinCircle(b, a)) {
    return a;
  }

  Point center;
  for (DimsType i = 0; i < kDim; ++i) {
    center.pnt[i] = (a.GetCenter().pnt[i] + b.GetCenter().pnt[i]) / 2;
  }
  Coord radius = P2PDistance(a.GetCenter(), b.GetCenter()) +
                 std::max(a.GetRadius(), b.GetRadius());

  return CircleType{center, CircleType::ComputeRadius(radius)};
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename CircleType>
inline CircleType BaseTree<Point, DerivedTree, kSkHeight,
                           kImbaRatio>::GetCircle(Point const& p,
                                                  CircleType const& cl) {
  if (PointWithinCircle(p, cl)) {
    return cl;
  }
  Coord radius = P2PDistance(p, cl.GetCenter());
  return CircleType{cl.GetCenter(), CircleType::ComputeRadius(radius)};
}

}  // namespace pspt

#endif  // PSPT_BASE_TREE_IMPL_CIRCLE_OP_HPP_
