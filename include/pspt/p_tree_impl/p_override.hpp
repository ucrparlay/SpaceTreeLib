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
  CpamAugMap::template knn<BT>(this->cpam_aug_map_, q, bq, logger);
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
  cpam_aug_map_.finish();
  return;
}

}  // namespace pspt

#endif  // PSPT_P_TREE_IMPL_P_OVERRIDE_HPP_
