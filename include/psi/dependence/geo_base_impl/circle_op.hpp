#ifndef PSI_DEPENDENCE_GEO_BASE_IMPL_CIRCLE_OP_HPP_
#define PSI_DEPENDENCE_GEO_BASE_IMPL_CIRCLE_OP_HPP_

#include <cmath>

#include "dependence/concepts.h"
#include "dependence/geo_base.h"

namespace psi {
template <class TypeTrait>
template <typename CircleType>
inline bool GeoBase<TypeTrait>::LegalCircle(CircleType const& cl) {
  return Num::Geq(cl.GetRadius(), 0);
}

// NOTE: point within the circle
template <class TypeTrait>
template <typename CircleType>
inline bool GeoBase<TypeTrait>::WithinCircle(Point const& p,
                                             CircleType const& cl) {
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
template <class TypeTrait>
template <typename CircleType>
inline bool GeoBase<TypeTrait>::WithinCircle(Box const& box,
                                             CircleType const& cl) {
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
    if (Num::Gt(r, cl.GetRadiusSquare())) return false;
  }
  assert(Num::Leq(r, cl.GetRadiusSquare()));
  return true;
}

template <class TypeTrait>
template <typename CircleType>
inline bool GeoBase<TypeTrait>::CircleIntersectBox(CircleType const& cl,
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

template <class TypeTrait>
template <typename CircleType>
inline CircleType GeoBase<TypeTrait>::GetCircle(Box const& box) {
  assert(LegalBox(box) && box != GetEmptyBox());

  Coord r = 0;
  if constexpr (kDim == 2) {
    r = Num::DivideTwoCeil(box.second.pnt[0] - box.first.pnt[0]) *
            Num::DivideTwoCeil(box.second.pnt[0] - box.first.pnt[0]) +
        Num::DivideTwoCeil(box.second.pnt[1] - box.first.pnt[1]) *
            Num::DivideTwoCeil(box.second.pnt[1] - box.first.pnt[1]);
  } else if constexpr (kDim == 3) {
    r = Num::DivideTwoCeil(box.second.pnt[0] - box.first.pnt[0]) *
            Num::DivideTwoCeil(box.second.pnt[0] - box.first.pnt[0]) +
        Num::DivideTwoCeil(box.second.pnt[1] - box.first.pnt[1]) *
            Num::DivideTwoCeil(box.second.pnt[1] - box.first.pnt[1]) +
        Num::DivideTwoCeil(box.second.pnt[2] - box.first.pnt[2]) *
            Num::DivideTwoCeil(box.second.pnt[2] - box.first.pnt[2]);
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      r += Num::DivideTwoCeil(box.second.pnt[i] - box.first.pnt[i]) *
           Num::DivideTwoCeil(box.second.pnt[i] - box.first.pnt[i]);
    }
  }

  // std::cout << "r: " << std::sqrt(r) << '\n';
  // assert(WithinCircle(
  //     box,
  //     CircleType{GetBoxCenter(box),
  //     CircleType::ComputeRadius(std::sqrt(r))}));

  return CircleType{GetBoxCenter(box), CircleType::ComputeRadius(std::sqrt(r))};
}

template <class TypeTrait>
template <typename CircleType>
inline CircleType GeoBase<TypeTrait>::GetCircle(Slice V) {
  return GetCircle<CircleType>(GetBox(V));
}

template <class TypeTrait>
template <typename CircleType>
inline CircleType GeoBase<TypeTrait>::GetCircle(
    parlay::sequence<CircleType> const& circle_seq) {
  CircleType circle = circle_seq[0];
  for (size_t i = 1; i < circle_seq.size(); i++) {
    circle = GetCircle(circle, circle_seq[i]);
  }
  return circle;
}

template <class TypeTrait>
template <typename CircleType1, typename CircleType2>
inline bool GeoBase<TypeTrait>::CircleWithinCircle(CircleType1 const& a,
                                                   CircleType2 const& b) {
  assert(LegalCircle(a) && LegalCircle(b));
  return Num::Leq(a.GetRadius(), b.GetRadius()) &&
         Num::Leq(
             P2PDistanceSquare(a.GetCenter(), b.GetCenter()),
             (b.GetRadius() - a.GetRadius()) * (b.GetRadius() - a.GetRadius()));
}

template <class TypeTrait>
template <typename CircleType>
inline CircleType GeoBase<TypeTrait>::GetCircle(CircleType const& a,
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
  double radius =
      std::sqrt(P2PDistanceSquare(a.GetCenter(), b.GetCenter())) / 2.0 +
      std::max(a.GetRadius(), b.GetRadius());

  return CircleType{center, CircleType::ComputeRadius(radius)};
}

template <class TypeTrait>
template <typename CircleType>
inline CircleType GeoBase<TypeTrait>::GetCircle(Point const& p,
                                                CircleType const& cl) {
  if (PointWithinCircle(p, cl)) {
    return cl;
  }
  Coord radius = P2PDistanceSquare(p, cl.GetCenter());
  return CircleType{cl.GetCenter(),
                    CircleType::ComputeRadius(std::sqrt(radius))};
}

}  // namespace psi

#endif  // PSI_DEPENDENCE_GEO_BASE_IMPL_CIRCLE_OP_HPP_