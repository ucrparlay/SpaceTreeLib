#ifndef PSI_POINTER_BASED_KD_TREE_IMPL_KD_OVERRIDE_HPP_
#define PSI_POINTER_BASED_KD_TREE_IMPL_KD_OVERRIDE_HPP_

#include <type_traits>
#include <utility>

#include "../../dependence/loggers.h"
#include "../../dependence/tree_node.h"
#include "../kd_tree.h"

namespace psi {
template <typename TypeTrait>
void KdTree<TypeTrait>::Compress2Multi() {
  this->root_ = BT::template Compress2Multi<KdInteriorNode, CompressInterior>(
      this->root_);
  return;
}

template <typename TypeTrait>
template <typename Range>
auto KdTree<TypeTrait>::KNN(Node* T, Point const& q,
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

template <typename TypeTrait>
template <typename Range>
void KdTree<TypeTrait>::Flatten(Range&& Out) {
  BT::template FlattenRec<Leaf, Interior>(this->root_, parlay::make_slice(Out));
  return;
}

template <typename TypeTrait>
auto KdTree<TypeTrait>::RangeCount(Box const& bx) {
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

template <typename TypeTrait>
template <typename Range>
auto KdTree<TypeTrait>::RangeQuery(Box const& query_box, Range&& Out) {
  RangeQueryLogger logger;
  size_t s = 0;
  BT::template RangeQuerySerialRecursive<Leaf, Interior>(
      this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_,
      logger);
  return std::make_pair(s, logger);
}

template <typename TypeTrait>
constexpr void KdTree<TypeTrait>::DeleteTree() {
  this->template DeleteTreeWrapper<Leaf, Interior>();
}

}  // namespace psi

#endif  // PSI_POINTER_BASED_KD_TREE_IMPL_KD_OVERRIDE_HPP_