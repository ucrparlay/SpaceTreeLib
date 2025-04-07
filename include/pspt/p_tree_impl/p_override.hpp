#ifndef PSPT_P_TREE_IMPL_P_OVERRIDE_HPP_
#define PSPT_P_TREE_IMPL_P_OVERRIDE_HPP_

#include <type_traits>
#include <utility>

#include "../p_tree.h"
#include "pspt/dependence/loggers.h"
#include "pspt/dependence/tree_node.h"

namespace pspt {

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
auto PTree<Point, SplitRule, kSkHeight, kImbaRatio>::KNN(
    Node* T, Point const& q, kBoundedQueue<Point, Range>& bq) {
  KNNLogger logger;
  // BT::template KNNMix<Leaf, KdInteriorNode, CompressInterior>(
  //     T, q, 0, 1, bq, this->tree_box_, logger);
  // BT::template KNNBinary<Leaf, Interior>(T, q, bq, this->tree_box_, logger);
  return logger;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::Flatten(Range&& Out) {
  // BT::template FlattenRec<Leaf, Interior>(this->root_,
  // parlay::make_slice(Out));
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
auto PTree<Point, SplitRule, kSkHeight, kImbaRatio>::RangeCount(
    Box const& query_box) {
  RangeQueryLogger logger;

  // auto f = [&](auto const& cur) { return mbr_mbr_relation(cur, query_mbr); };
  // auto f2 = [&](auto const& cur) { return point_in_mbr(cur, query_mbr); };
  // auto res = CpamAugMap::range_count_filter2<BaseTree>(zCPAM, f, f2);
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
  cpam_aug_map_.finish();
  return;
}

}  // namespace pspt

#endif  // PSPT_P_TREE_IMPL_P_OVERRIDE_HPP_
