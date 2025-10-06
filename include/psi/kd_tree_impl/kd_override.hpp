#ifndef PSI_KD_TREE_IMPL_KD_OVERRIDE_HPP_
#define PSI_KD_TREE_IMPL_KD_OVERRIDE_HPP_

#include <type_traits>
#include <utility>

#include "../kd_tree.h"
#include "psi/dependence/loggers.h"
#include "psi/dependence/tree_node.h"

#define KDTREE_TEMPLATE                                               \
  template <typename Point, typename SplitRule, typename LeafAugType, \
            typename InteriorAugType, uint_fast8_t kSkHeight,         \
            uint_fast8_t kImbaRatio>
#define KDTREE_CLASS \
  KdTree<Point, SplitRule, LeafAugType, InteriorAugType, kSkHeight, kImbaRatio>

namespace psi {
KDTREE_TEMPLATE
void KDTREE_CLASS::Compress2Multi() {
  this->root_ = BT::template Compress2Multi<KdInteriorNode, CompressInterior>(
      this->root_);
  return;
}

KDTREE_TEMPLATE
template <typename Range>
auto KDTREE_CLASS::KNN(Node* T, Point const& q,
                       kBoundedQueue<Point, Range>& bq) {
  KNNLogger logger;
  // BT::template KNNMix<Leaf, KdInteriorNode, CompressInterior>(
  //     T, q, 0, 1, bq, this->tree_box_, logger);
  if constexpr (HasBox<typename Interior::AT>) {
    BT::template KNNBinaryBox<Leaf, Interior>(T, q, bq, logger);
  } else {
    BT::template KNNBinary<Leaf, Interior>(T, q, bq, this->tree_box_, logger);
  }
  return logger;
}

KDTREE_TEMPLATE
template <typename Range>
void KDTREE_CLASS::Flatten(Range&& Out) {
  BT::template FlattenRec<Leaf, Interior>(this->root_, parlay::make_slice(Out));
  return;
}

KDTREE_TEMPLATE
auto KDTREE_CLASS::RangeCount(Box const& bx) {
  RangeQueryLogger logger;
  size_t size = BT::template RangeCountRectangle<Leaf, Interior>(
      this->root_, bx, this->tree_box_, logger);
  return std::make_pair(size, logger);
}

// template <typename Point, typename SplitRule, typename LeafAugType, typename
// InteriorAugType,  uint_fast8_t kSkHeight,
//           uint_fast8_t kImbaRatio>
// auto KdTree<Point, SplitRule, LeafAugType, InteriorAugType, kSkHeight ,
// kImbaRatio>::RangeCount(
//     Circle const& cl) {
//   return BT::template RangeCountRadius<Leaf, Interior>(this->root_, cl,
//                                                        this->tree_box_);
// }

KDTREE_TEMPLATE
template <typename Range>
auto KDTREE_CLASS::RangeQuery(Box const& query_box, Range&& Out) {
  RangeQueryLogger logger;
  size_t s = 0;
  BT::template RangeQuerySerialRecursive<Leaf, Interior>(
      this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_,
      logger);
  return std::make_pair(s, logger);
}

KDTREE_TEMPLATE
constexpr void KDTREE_CLASS::DeleteTree() {
  BT::template DeleteTreeWrapper<Leaf, Interior>();
}

}  // namespace psi

#undef KDTREE_TEMPLATE
#undef KDTREE_CLASS

#endif  // PSI_KD_TREE_IMPL_KD_OVERRIDE_HPP_
