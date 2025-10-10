#ifndef PSI_DEPENDENCE_TYPE_TRAIT_H_
#define PSI_DEPENDENCE_TYPE_TRAIT_H_

#include <cstdint>
#include <type_traits>
#include <utility>

#include "basic_point.h"
#include "comparator.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"

namespace psi {

template <typename PointType, typename NodeType = void,
          typename SplitRuleType = void, typename LeafAug = void,
          typename InteriorAug = void, uint_fast8_t kSkHeightVal = 6,
          uint_fast8_t kImbaRatioVal = 30>
class TypeTrait {
 public:
  static constexpr uint_fast8_t const kSkHeight = kSkHeightVal;
  static constexpr uint_fast8_t const kImbaRatio = kImbaRatioVal;

  using Point = PointType;
  using BasicPoint = typename Point::BP;
  using TemplatePoint = Point;
  using Coord = typename Point::Coord;
  using Coords = typename Point::Coords;
  using DisType = typename Point::DisType;

  using DimsType = uint_fast8_t;
  using DepthType = int;
  using BucketType =
      std::conditional_t<(kSkHeight > 7), uint_fast16_t, uint_fast8_t>;
  using BallsType = uint_fast32_t;
  using IDType = uint_fast32_t;

  using Num = Num_Comparator<Coord>;
  using Slice = parlay::slice<Point*, Point*>;
  using ConstSlice = parlay::slice<Point const*, Point const*>;
  using Points = parlay::sequence<Point>;
  using ConstPoints = parlay::sequence<Point> const;
  using PointsIter = typename parlay::sequence<Point>::iterator;
  using BucketSeq = parlay::sequence<BucketType>;
  using BallSeq = parlay::sequence<BallsType>;

  using HyperPlane = std::pair<Coord, DimsType>;
  using HyperPlaneSeq = parlay::sequence<HyperPlane>;
  using Box = std::pair<BasicPoint, BasicPoint>;
  using BoxSeq = parlay::sequence<Box>;

  using Node = NodeType;
  using NodeBoolean = std::pair<Node*, bool>;
  using NodeBox = std::pair<Node*, Box>;
  using NodeBoxSeq = parlay::sequence<NodeBox>;
  using NodeTag = std::pair<Node*, uint_fast8_t>;
  using NodeTagSeq = parlay::sequence<NodeTag>;

  using LeafAugType = LeafAug;
  using InteriorAugType = InteriorAug;
  using SplitRule = SplitRuleType;
};

}  // namespace psi

#endif  // PSI_DEPENDENCE_TYPE_TRAIT_H_