#ifndef PSI_POINTER_VIEW_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_
#define PSI_POINTER_VIEW_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_

#include <tuple>
#include <type_traits>
#include <utility>

#include "../../dependence/concepts.h"
#include "../../dependence/tree_node.h"
#include "../orth_tree.h"

namespace psi {
namespace pointer_view {
template <typename TypeTrait>
template <typename Range>
auto OrthTree<TypeTrait>::KNN(Node* T, Point const& q,
                         kBoundedQueue<Point, Range>& bq) {
  KNNLogger logger;
  // BT::template KNNMulti<Leaf, Interior>(T, q, DIM, bq, this->tree_box_,
  //                                       vis_node_num, generate_box_num,
  //                                       check_box_num);
  // BT::template KNNBinary<Leaf, KdInteriorNode>(T, q, DIM, bq,
  // this->tree_box_,
  //                                              logger);
  if constexpr (HasBox<typename Interior::AT>) {
    BT::template KNNMulti<Leaf, Interior>(T, q, bq, logger);
    // BT::template KNNMultiExpandBox<Leaf, Interior>(T, q, 0, 1, bq, logger);
  } else {
    BT::template KNNMultiExpand<Leaf, Interior>(T, q, 0, 1, bq, this->tree_box_,
                                                logger);
  }
  return logger;
}

template <typename TypeTrait>
template <typename Range>
void OrthTree<TypeTrait>::Flatten(Range&& Out) {
  BT::template FlattenRec<Leaf, Interior>(this->root_, parlay::make_slice(Out));
}

template <typename TypeTrait>
auto OrthTree<TypeTrait>::RangeCount(Box const& bx) {
  RangeQueryLogger logger;
  size_t s = BT::template RangeCountRectangle<Leaf, Interior>(
      this->root_, bx, this->tree_box_, 0, 1, logger);
  return std::make_pair(s, logger);
}

// template <typename Point, typename SplitRule, typename LeafAugType, typename
// InteriorAugType, uint_fast8_t kMD,
//           uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
// auto OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::RangeCount(
//     Circle const& cl) {
//   return BT::template RangeCountRadius<Leaf, Interior>(this->root_, cl,
//                                                        this->tree_box_);
// }

template <typename TypeTrait>
template <typename Range>
auto OrthTree<TypeTrait>::RangeQuery(Box const& query_box, Range&& Out) {
  RangeQueryLogger logger;
  size_t s = 0;
  BT::template RangeQuerySerialRecursive<Leaf, Interior>(
      this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_, 0, 1,
      logger);
  return std::make_pair(s, logger);
}

template <typename TypeTrait>
constexpr void OrthTree<TypeTrait>::DeleteTree() {
  BT::template DeleteTreeWrapper<Leaf, Interior>();
  this->fixed_box = false;
}

}  // namespace pointer_view
}  // namespace psi

 
 

#endif  // PSI_POINTER_VIEW_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_
