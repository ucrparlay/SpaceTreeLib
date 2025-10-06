#ifndef PSI_R_TREE_IMPL_R_OVERRIDE_HPP_
#define PSI_R_TREE_IMPL_R_OVERRIDE_HPP_

#include <utility>

#include "../r_tree.h"

#define RTREE_TEMPLATE template <typename Point, typename SplitRule, \
    uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
#define RTREE_CLASS RTree<Point, SplitRule, kSkHeight, kImbaRatio>

#include "psi/dependence/loggers.h"
#include "psi/dependence/tree_node.h"

namespace psi {
RTREE_TEMPLATE

template <typename Range>
auto RTREE_CLASS::KNN(
    Node* T, Point const& q, kBoundedQueue<Point, Range>& bq) {
  KNNLogger logger;
  BT::template KNNBinary<Leaf, Interior>(T, q, bq, logger);
  return logger;
}

RTREE_TEMPLATE

template <typename Range>
void RTREE_CLASS::Flatten(Range&& Out) {
  BT::template FlattenRec<Leaf, Interior>(this->root_, parlay::make_slice(Out));
  return;
}

RTREE_TEMPLATE

auto RTREE_CLASS::RangeCount(
    Box const& query_box) {
  RangeQueryLogger logger;
  size_t size = BT::template RangeCountRectangle<Leaf, Interior>(
      this->root_, query_box, this->tree_box_, logger);
  return std::make_pair(size, logger);
}

RTREE_TEMPLATE

auto RTREE_CLASS::RangeCount(
    Circle const& cl) {
  return BT::template RangeCountRadius<Leaf, Interior>(this->root_, cl,
                                                       this->tree_box_);
}

RTREE_TEMPLATE

template <typename Range>
auto RTREE_CLASS::RangeQuery(
    Box const& query_box, Range&& Out) {
  RangeQueryLogger logger;
  size_t s = 0;
  BT::template RangeQuerySerialRecursive<Leaf, Interior>(
      this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_,
      logger);
  return std::make_pair(s, logger);
}

RTREE_TEMPLATE

constexpr void RTREE_CLASS::DeleteTree() {
  BT::template DeleteTreeWrapper<Leaf, Interior>();
}


#undef RTREE_TEMPLATE
#undef RTREE_CLASS
}  // namespace psi

#endif  // PSI_R_TREE_IMPL_R_OVERRIDE_HPP_
