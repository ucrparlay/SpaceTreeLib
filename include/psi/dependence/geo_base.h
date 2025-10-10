#ifndef PSI_DEPENDENCE_GEO_BASE_H_
#define PSI_DEPENDENCE_GEO_BASE_H_

#include <sys/types.h>

#include <cstdint>
#include <cstdio>
#include <type_traits>

#include "comparator.h"
#include "concepts.h"
#include "tree_node.h"
#include "type_trait.h"

namespace psi {

template <class TypeTrait>
class GeoBase {
 public:
  using Box = typename TypeTrait::Box;
  using BoxSeq = typename TypeTrait::BoxSeq;
  using Coord = typename TypeTrait::Coord;
  using DimsType = typename TypeTrait::DimsType;
  using DisType = typename TypeTrait::DisType;
  using Node = typename TypeTrait::Node;
  using Point = typename TypeTrait::Point;
  using Slice = typename TypeTrait::Slice;
  using Num = typename TypeTrait::Num;
  using DepthType = int;

  static constexpr uint_fast8_t kDim = Point::GetDim();
  static constexpr uint_fast16_t kSerialBuildCutoff = 1024;

  struct NormalCircle {
    Point const& GetCenter() const { return center; }

    Coord GetRadius() const { return radius; }

    Coord GetRadiusSquare() const { return radius * radius; }

    static Coord ComputeRadius(double const& r) {
      return static_cast<Coord>(r);
    }

    friend std::ostream& operator<<(std::ostream& o, NormalCircle const& cl) {
      o << "{ " << cl.center << ", " << cl.radius << "}";
      return o;
    }

    Point center;
    Coord radius;
  };

  struct CoverCircle {
    Point const& GetCenter() const { return center; }

    Coord GetRadius() const { return static_cast<Coord>(1 << level); }

    Coord GetRadiusSquare() const { return GetRadius() * GetRadius(); }

    static DepthType ComputeRadius(auto const& r) {
      return Num_Comparator<decltype(r)>::IsZero(r)
                 ? -1
                 : static_cast<DepthType>(std::ceil(std::log2(r)));
    }

    bool operator==(CoverCircle const& cl) const {
      return center == cl.center && level == cl.level;
    }

    friend std::ostream& operator<<(std::ostream& o, CoverCircle const& cl) {
      o << "{ " << cl.center << ", " << cl.level << "}";
      return o;
    }

    Point center;
    DepthType level;
  };

  //============================================================================
  // SECTION 9: GEOMETRY UTILITIES - BOX OPERATIONS
  // Implementation: base_tree_impl/box_op.hpp
  //============================================================================

  static inline Coord GetBoxMid(DimsType const d, Box const& bx);
  static inline bool LegalBox(Box const& bx);
  static inline bool WithinBox(Box const& a, Box const& b);
  static inline bool SameBox(Box const& a, Box const& b);
  static inline bool WithinBox(Point const& p, Box const& bx);
  static inline bool BoxIntersectBox(Box const& a, Box const& b);
  static inline bool IsBoxLineInDimension(Box const& box, DimsType d);
  static inline bool VerticalLineSplitBox(Coord const& l, Box const& box,
                                          DimsType d);
  static inline bool VerticalLineOnBoxLeftEdge(Coord const& l, Box const& box,
                                               DimsType d);
  static inline bool VerticalLineOnBoxRightEdge(Coord const& l, Box const& box,
                                                DimsType d);
  static inline bool VerticalLineOnBoxEdge(Coord const& l, Box const& box,
                                           DimsType d);
  static inline bool VerticalLineIntersectBox(Coord const& l, Box const& box,
                                              DimsType d);
  static inline bool VerticalLineIntersectBoxExclude(Coord const& l,
                                                     Box const& box,
                                                     DimsType d);

  static inline Box GetEmptyBox();

  static inline Point GetBoxCenter(Box const& box);

  static Box GetBox(Box const& x, Box const& y);

  template <typename SliceType>
  static Box GetBoxFromSlice(SliceType const V);

  template <typename Range>
  static Box GetBox(Range&& range)
    requires(parlay::is_random_access_range_v<Range>);

  static Box GetBoxFromBoxSeq(BoxSeq const& box_seq);

  template <typename Leaf, typename Interior>
  static Box GetBox(Node* T);

  //============================================================================
  // SECTION 10: GEOMETRY UTILITIES - CIRCLE OPERATIONS
  // Implementation: base_tree_impl/circle_op.hpp
  //============================================================================

  template <typename CircleType>
  static bool LegalCircle(CircleType const& cl);

  template <typename CircleType>
  static inline bool WithinCircle(Point const& p, CircleType const& cl);

  template <typename CircleType>
  static inline bool WithinCircle(Box const& box, CircleType const& cl);

  template <typename CircleType1, typename CircleType2>
  static inline bool CircleWithinCircle(CircleType1 const& a,
                                        CircleType2 const& b);

  template <typename CircleType>
  static inline bool CircleIntersectBox(CircleType const& cl, Box const& box);

  template <typename CircleType>
  static inline bool CircleIntersectCircle(CircleType const& a,
                                           CircleType const& b);

  template <typename CircleType>
  static inline CircleType GetCircle(Box const& box);

  template <typename CircleType>
  static inline CircleType GetCircle(Slice V);

  template <typename CircleType>
  static inline CircleType GetCircle(CircleType const& a, CircleType const& b);

  template <typename CircleType>
  static inline CircleType GetCircle(Point const& p, CircleType const& cl);

  template <typename CircleType>
  static inline CircleType GetCircle(
      parlay::sequence<CircleType> const& circle_seq);

  //============================================================================
  // SECTION 11: GEOMETRY UTILITIES - DISTANCE METRICS
  // Implementation: base_tree_impl/knn_query.hpp
  //============================================================================

  static inline DisType P2PDistanceSquare(Point const& p, Point const& q);

  static inline DisType P2BMinDistanceSquare(Point const& p, Box const& a);

  static inline DisType P2BMaxDistanceSquare(Point const& p, Box const& a);

  static inline double P2CMinDistance(Point const& p, Point const& center,
                                      DisType const r);

  template <typename CircleType>
  static inline double P2CMinDistance(Point const& p, CircleType const& cl);

  static inline DisType InterruptibleDistance(Point const& p, Point const& q,
                                              DisType up);
};

}  // namespace psi

#include "geo_base_impl/box_op.hpp"       // Section 9
#include "geo_base_impl/circle_op.hpp"    // Section 10
#include "geo_base_impl/distance_op.hpp"  // Section 11

#endif  // PSI_DEPENDENCE_GEO_BASE_H_