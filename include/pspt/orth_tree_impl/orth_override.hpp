#ifndef PSPT_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_
#define PSPT_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_

#include <tuple>
#include <type_traits>
#include <utility>

#include "../orth_tree.h"
#include "pspt/dependence/tree_node.h"

namespace pspt {
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
template <typename Range>
auto OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::KNN(
    Node* T, Point const& q, kBoundedQueue<Point, Range>& bq) {
  KNNLogger logger;
  // BT::template KNNMulti<Leaf, Interior>(T, q, DIM, bq, this->tree_box_,
  //                                       vis_node_num, generate_box_num,
  //                                       check_box_num);
  // BT::template KNNBinary<Leaf, KdInteriorNode>(T, q, DIM, bq,
  // this->tree_box_,
  //                                              logger);
  BT::template KNNMultiExpand<Leaf, Interior>(T, q, 0, 1, bq, this->tree_box_,
                                              logger);
  return logger;
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
template <typename Range>
void OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::Flatten(
    Range&& Out) {
  BT::template FlattenRec<Leaf, Interior>(this->root_, parlay::make_slice(Out));
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
auto OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::RangeCount(
    Box const& bx) {
  RangeQueryLogger logger;
  size_t s = BT::template RangeCountRectangle<Leaf, Interior>(
      this->root_, bx, this->tree_box_, 0, 1, logger);
  return std::make_pair(s, logger);
}

// template <typename Point, typename SplitRule, uint_fast8_t kMD,
//           uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
// auto OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::RangeCount(
//     Circle const& cl) {
//   return BT::template RangeCountRadius<Leaf, Interior>(this->root_, cl,
//                                                        this->tree_box_);
// }

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
template <typename Range>
auto OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::RangeQuery(
    Box const& query_box, Range&& Out) {
  RangeQueryLogger logger;
  size_t s = 0;
  BT::template RangeQuerySerialRecursive<Leaf, Interior>(
      this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_, 0, 1,
      logger);
  return std::make_pair(s, logger);
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
constexpr void
OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::DeleteTree() {
  BT::template DeleteTreeWrapper<Leaf, Interior>();
}

}  // namespace pspt

#endif  // PSPT_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_
