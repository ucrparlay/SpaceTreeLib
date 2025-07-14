#ifndef PSPT_P_TREE_IMPL_P_OVERRIDE_HPP_
#define PSPT_P_TREE_IMPL_P_OVERRIDE_HPP_

#include <type_traits>
#include <utility>

#include "../p_tree.h"
#include "pspt/dependence/loggers.h"
#include "pspt/dependence/tree_node.h"

namespace pspt {
#define FT long

// return the sqr distance between a point and a mbr
template <class Point, class MBR>
auto point_mbr_sqrdis(Point p, MBR& mbr) {
  FT dx = max(max(mbr.first[0] - p[0], (FT)0.0),
              max(p[0] - mbr.second[0], (FT)0.0));
  FT dy = max(max(mbr.first[1] - p[1], (FT)0.0),
              max(p[1] - mbr.second[1], (FT)0.0));
  return dx * dx + dy * dy;
}

// return the sqr distance between two points
template <class Point>
auto point_point_sqrdis(Point const& lhs, Point const& rhs) {
  return (lhs.pnt[0] - rhs.pnt[0]) * (lhs.pnt[0] - rhs.pnt[0]) +
         (lhs.pnt[1] - rhs.pnt[1]) * (lhs.pnt[1] - rhs.pnt[1]);
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
auto PTree<Point, SplitRule, kSkHeight, kImbaRatio>::KNN(
    Node* T, Point const& q, kBoundedQueue<Point, Range>& bq) {
  KNNLogger logger;
  auto f = [&](auto const cur_pt) { return BT::P2PDistanceSquare(cur_pt, q); };

  auto f2 = [&](Box const cur_mbr) {
    return BT::P2BMinDistanceSquare(q, cur_mbr);
  };

  // using nn_pair = std::pair<Point, FT>;
  // struct nn_pair_cmp {
  //   bool operator()(nn_pair& lhs, nn_pair& rhs) {
  //     return lhs.second < rhs.second ||
  //            (lhs.second == rhs.second && lhs.first < rhs.first);
  //   }
  // };
  //
  // std::priority_queue<nn_pair, std::vector<nn_pair>, nn_pair_cmp> nn_res;
  CpamAugMap::template knn_filter<BT>(this->cpam_aug_map_, q, bq, logger);
  // return nn_res.top().second;

  // CpamAugMap::template knn<BT>(this->cpam_aug_map_, q, bq, logger);
  return logger;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::Flatten(Range&& Out) {
  CpamAugMap::entries(parlay::make_slice(Out));
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
auto PTree<Point, SplitRule, kSkHeight, kImbaRatio>::RangeCount(
    Box const& query_box) {
  RangeQueryLogger logger;

  auto size = CpamAugMap::template range_count_filter2<BT>(this->cpam_aug_map_,
                                                           query_box, logger);

  return std::make_pair(size, logger);
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
auto PTree<Point, SplitRule, kSkHeight, kImbaRatio>::RangeQuery(
    Box const& query_box, Range&& Out) {
  RangeQueryLogger logger;
  size_t cnt = 0;
  CpamAugMap::template range_report_filter2<BT>(
      this->cpam_aug_map_, query_box, cnt, parlay::make_slice(Out), logger);
  return std::make_pair(cnt, logger);
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
constexpr void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::DeleteTree() {
  cpam_aug_map_.clear();
  return;
}

}  // namespace pspt

#endif  // PSPT_P_TREE_IMPL_P_OVERRIDE_HPP_
