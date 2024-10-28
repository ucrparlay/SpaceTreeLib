#include <type_traits>
#include <utility>

#include "../kd_tree.h"
#include "cpdd/dependence/loggers.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template <typename Point, typename SplitRule, uint_fast8_t kBDO>
template <typename Range>
auto KdTree<Point, SplitRule, kBDO>::KNN(Node* T, Point const& q,
                                         kBoundedQueue<Point, Range>& bq) {
  KNNLogger logger;
  BT::template KNNBinary<Leaf, Interior>(T, q, bq, this->tree_box_, logger);
  return logger;
}

template <typename Point, typename SplitRule, uint_fast8_t kBDO>
template <typename Range>
void KdTree<Point, SplitRule, kBDO>::Flatten(Range&& Out) {
  BT::template FlattenRec<Leaf, Interior>(this->root_, parlay::make_slice(Out));
}

template <typename Point, typename SplitRule, uint_fast8_t kBDO>
auto KdTree<Point, SplitRule, kBDO>::RangeCount(Box const& bx) {
  RangeQueryLogger logger;
  size_t size = BT::template RangeCountRectangle<Leaf, Interior>(
      this->root_, bx, this->tree_box_, logger);
  return std::make_pair(size, logger);
}

template <typename Point, typename SplitRule, uint_fast8_t kBDO>
auto KdTree<Point, SplitRule, kBDO>::RangeCount(Circle const& cl) {
  return BT::template RangeCountRadius<Leaf, Interior>(this->root_, cl,
                                                       this->tree_box_);
}

template <typename Point, typename SplitRule, uint_fast8_t kBDO>
template <typename Range>
auto KdTree<Point, SplitRule, kBDO>::RangeQuery(Box const& query_box,
                                                Range&& Out) {
  RangeQueryLogger logger;
  size_t s = 0;
  BT::template RangeQuerySerialRecursive<Leaf, Interior>(
      this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_,
      logger);
  return std::make_pair(s, logger);
}

template <typename Point, typename SplitRule, uint_fast8_t kBDO>
constexpr void KdTree<Point, SplitRule, kBDO>::DeleteTree() {
  BT::template DeleteTreeWrapper<Leaf, Interior>();
}

}  // namespace cpdd
