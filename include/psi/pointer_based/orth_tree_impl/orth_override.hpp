#ifndef PSI_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_
#define PSI_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_

#include <tuple>
#include <type_traits>
#include <utility>

#include "../orth_tree.h"
#include "../../dependence/concepts.h"
#include "../../dependence/tree_node.h"

#define ORTHTREE_TEMPLATE                                             \
  template <typename Point, typename SplitRule, typename LeafAugType, \
            typename InteriorAugType, uint_fast8_t kMD,               \
            uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
#define ORTHTREE_CLASS                                                     \
  OrthTree<Point, SplitRule, LeafAugType, InteriorAugType, kMD, kSkHeight, \
           kImbaRatio>

namespace psi {
ORTHTREE_TEMPLATE
template <typename Range>
auto ORTHTREE_CLASS::KNN(Node* T, Point const& q,
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

ORTHTREE_TEMPLATE
template <typename Range>
void ORTHTREE_CLASS::Flatten(Range&& Out) {
  BT::template FlattenRec<Leaf, Interior>(this->root_, parlay::make_slice(Out));
}

ORTHTREE_TEMPLATE
auto ORTHTREE_CLASS::RangeCount(Box const& bx) {
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

ORTHTREE_TEMPLATE
template <typename Range>
auto ORTHTREE_CLASS::RangeQuery(Box const& query_box, Range&& Out) {
  RangeQueryLogger logger;
  size_t s = 0;
  BT::template RangeQuerySerialRecursive<Leaf, Interior>(
      this->root_, parlay::make_slice(Out), s, query_box, this->tree_box_, 0, 1,
      logger);
  return std::make_pair(s, logger);
}

ORTHTREE_TEMPLATE
constexpr void ORTHTREE_CLASS::DeleteTree() {
  BT::template DeleteTreeWrapper<Leaf, Interior>();
  this->fixed_box = false;
}

}  // namespace psi

#undef ORTHTREE_TEMPLATE
#undef ORTHTREE_CLASS

#endif  // PSI_ORTH_TREE_IMPL_ORTH_OVERRIDE_HPP_
