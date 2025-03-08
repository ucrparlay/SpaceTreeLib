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
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::Compress2Multi() {
  this->root_ = BT::template Compress2Multi<KdInteriorNode, CompressInterior>(
      this->root_);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
auto PTree<Point, SplitRule, kSkHeight, kImbaRatio>::KNN(
    Node* T, Point const& q, kBoundedQueue<Point, Range>& bq) {
  KNNLogger logger;
  // BT::template KNNMix<Leaf, KdInteriorNode, CompressInterior>(
  //     T, q, 0, 1, bq, this->tree_box_, logger);
  BT::template KNNBinary<Leaf, Interior>(T, q, bq, this->tree_box_, logger);
  return logger;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::Flatten(Range&& Out) {
  BT::template FlattenRec<Leaf, Interior>(this->root_, parlay::make_slice(Out));
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
auto PTree<Point, SplitRule, kSkHeight, kImbaRatio>::RangeCount(Box const& bx) {
  RangeQueryLogger logger;
  size_t size = BT::template RangeCountRectangle<Leaf, Interior>(
      this->root_, bx, this->tree_box_, logger);
  return std::make_pair(size, logger);
}

// template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
//           uint_fast8_t kImbaRatio>
// auto PTree<Point, SplitRule, kSkHeight, kImbaRatio>::RangeCount(
//     Circle const& cl) {
//   return BT::template RangeCountRadius<Leaf, Interior>(this->root_, cl,
//                                                        this->tree_box_);
// }

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
auto PTree<Point, SplitRule, kSkHeight, kImbaRatio>::RangeQuery(
    Box const& query_box, Range&& Out) {
  RangeQueryLogger logger;
  size_t s = 0;
  BT::template RangeQuerySerialRecursive<Leaf, Interior>(
      this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_,
      logger);
  return std::make_pair(s, logger);
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
constexpr void PTree<Point, SplitRule, kSkHeight, kImbaRatio>::DeleteTree() {
  BT::template DeleteTreeWrapper<Leaf, Interior>();
}

}  // namespace pspt

#endif  // PSPT_P_TREE_IMPL_P_OVERRIDE_HPP_
