#ifndef PSI_DEPENDENCE_GEO_BASE_IMPL_DISTANCE_OP_HPP_
#define PSI_DEPENDENCE_GEO_BASE_IMPL_DISTANCE_OP_HPP_

#include <algorithm>
#include <utility>

#include "dependence/geo_base.h"

namespace psi {

// NOTE: distance between two Points
// TODO: change the name to P2PDistanceSquare to avoid ambiguous
template <class TypeTrait>
inline typename GeoBase<TypeTrait>::DisType
GeoBase<TypeTrait>::P2PDistanceSquare(Point const& p,
                                                    Point const& q) {
  constexpr uint_fast8_t kDim = Point::GetDim();
  DisType r = 0;

  if constexpr (kDim == 2) {
    DisType x = static_cast<DisType>(p.pnt[0]) - static_cast<DisType>(q.pnt[0]);
    DisType y = static_cast<DisType>(p.pnt[1]) - static_cast<DisType>(q.pnt[1]);
    r = x * x + y * y;
  } else if constexpr (kDim == 3) {
    DisType x = static_cast<DisType>(p.pnt[0]) - static_cast<DisType>(q.pnt[0]);
    DisType y = static_cast<DisType>(p.pnt[1]) - static_cast<DisType>(q.pnt[1]);
    DisType z = static_cast<DisType>(p.pnt[2]) - static_cast<DisType>(q.pnt[2]);
    r = x * x + y * y + z * z;
  } else {
    for (DimsType i = 0; i < kDim; ++i) {
      r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
    }
  }
  return r;
}

// NOTE: Distance between a Point and a Box
// return 0 when p is inside the box a
template <class TypeTrait>
inline typename GeoBase<TypeTrait>::DisType
GeoBase<TypeTrait>::P2BMinDistanceSquare(
    Point const& p, typename GeoBase<TypeTrait>::Box const& a) {
  DisType r = 0;
  // NOTE: the distance is 0 when p is inside the box
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(p.pnt[i], a.first.pnt[i])) {
      r += (static_cast<DisType>(a.first.pnt[i]) -
            static_cast<DisType>(p.pnt[i])) *
           (static_cast<DisType>(a.first.pnt[i]) -
            static_cast<DisType>(p.pnt[i]));
    } else if (Num::Gt(p.pnt[i], a.second.pnt[i])) {
      r += (static_cast<DisType>(p.pnt[i]) -
            static_cast<DisType>(a.second.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) -
            static_cast<DisType>(a.second.pnt[i]));
    } else {  // will not count the dis if p is inside the box in dimension i
      ;
    }
  }
  return r;
}

// NOTE: Max distance between a Point and a Box
template <class TypeTrait>
inline typename GeoBase<TypeTrait>::DisType
GeoBase<TypeTrait>::P2BMaxDistanceSquare(
    Point const& p, typename GeoBase<TypeTrait>::Box const& a) {
  DisType r = 0;
  for (DimsType i = 0; i < kDim; ++i) {
    if (Num::Lt(p.pnt[i], (a.second.pnt[i] + a.first.pnt[i]) / 2)) {
      r += (static_cast<DisType>(a.second.pnt[i]) -
            static_cast<DisType>(p.pnt[i])) *
           (static_cast<DisType>(a.second.pnt[i]) -
            static_cast<DisType>(p.pnt[i]));
    } else {
      r += (static_cast<DisType>(p.pnt[i]) -
            static_cast<DisType>(a.first.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) -
            static_cast<DisType>(a.first.pnt[i]));
    }
  }
  return r;
}

template <class TypeTrait>
inline double GeoBase<TypeTrait>::P2CMinDistance(
    Point const& p, Point const& center, DisType const r) {
  // return Num_Comparator<double>::Max(
  //     0.0, std::sqrt(P2PDistanceSquare(p, center)) - static_cast<double>(r));
  return std::sqrt(P2PDistanceSquare(p, center)) - static_cast<double>(r);
}

template <class TypeTrait>
template <typename CircleType>
inline double GeoBase<TypeTrait>::P2CMinDistance(
    Point const& p, CircleType const& cl) {
  return P2CMinDistance(p, cl.GetCenter(), cl.GetRadius());
}

// NOTE: early return the partial distance between p and q if it is larger than
// r else return the distance between p and q
template <class TypeTrait>
inline typename GeoBase<TypeTrait>::DisType
GeoBase<TypeTrait>::InterruptibleDistance(Point const& p,
                                                        Point const& q,
                                                        DisType up) {
  DisType r = 0;
  DimsType i = 0;
  if (kDim >= 6) {
    while (1) {
      r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
      ++i;
      r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
      ++i;
      r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
      ++i;
      r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
           (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
      ++i;

      if (Num_Comparator<DisType>::Gt(r, up)) {
        return r;
      }
      if (i + 4 > kDim) {
        break;
      }
    }
  }
  while (i < kDim) {
    r += (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i])) *
         (static_cast<DisType>(p.pnt[i]) - static_cast<DisType>(q.pnt[i]));
    ++i;
  }
  return r;
}

}  // namespace psi

#endif  // PSI_DEPENDENCE_GEO_BASE_IMPL_DISTANCE_OP_HPP_