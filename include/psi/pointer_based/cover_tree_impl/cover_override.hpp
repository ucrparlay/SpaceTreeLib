#ifndef PSI_COVER_TREE_IMPL_COVER_OVERRIDE_HPP_
#define PSI_COVER_TREE_IMPL_COVER_OVERRIDE_HPP_

#define COVERTREE_TEMPLATE                                              \
  template <typename Point, typename SplitRule, uint_fast8_t kSkHeight, \
            uint_fast8_t kImbaRatio>
#define COVERTREE_CLASS CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>

#include <tuple>
#include <type_traits>
#include <utility>

#include "../cover_tree.h"
#include "../../dependence/tree_node.h"

namespace psi {
COVERTREE_TEMPLATE
template <typename Range>
auto COVERTREE_CLASS::KNN(Node* T, Point const& q,
                          kBoundedQueue<Point, Range>& bq) {
  KNNLogger logger;
  BT::template KNNCover<Leaf, Interior>(T, q, bq, logger);
  return logger;
}

COVERTREE_TEMPLATE
template <typename Range>
void COVERTREE_CLASS::Flatten(Range&& Out) {
  BT::template FlattenRec<Leaf, Interior>(this->root_, parlay::make_slice(Out));
}

COVERTREE_TEMPLATE
auto COVERTREE_CLASS::RangeCount(Box const& bx) {
  RangeQueryLogger logger;
  size_t s = BT::template RangeCountRectangle<Leaf, Interior>(
      this->root_, bx, this->tree_box_, 0, 1, logger);
  return std::make_pair(s, logger);
}

// template <typename Point, typename SplitRule,
//           uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
// auto COVERTREE_CLASS::RangeCount(
//     Circle const& cl) {
//   return BT::template RangeCountRadius<Leaf, Interior>(this->root_, cl,
//                                                        this->tree_box_);
// }

COVERTREE_TEMPLATE
template <typename Range>
auto COVERTREE_CLASS::RangeQuery(Box const& query_box, Range&& Out) {
  RangeQueryLogger logger;
  size_t s = 0;
  BT::template RangeQuerySerialRecursive<Leaf, Interior>(
      this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_, 0, 1,
      logger);
  return std::make_pair(s, logger);
}

COVERTREE_TEMPLATE
constexpr void COVERTREE_CLASS::DeleteTree() {
  BT::template DeleteTreeWrapper<Leaf, Interior>();
}

}  // namespace psi

#undef COVERTREE_TEMPLATE
#undef COVERTREE_CLASS

#endif
