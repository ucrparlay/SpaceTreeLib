#ifndef PSTP_BASE_TREE_IMPL_TREE_OP_BUILD_INNER_TREE_HPP_
#define PSTP_BASE_TREE_IMPL_TREE_OP_BUILD_INNER_TREE_HPP_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#include "../../base_tree.h"

namespace pstp {
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsBinaryNode Interior, typename Base>
Base BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::BuildInnerTree(
    BucketType idx, HyperPlaneSeq const& pivots,
    parlay::sequence<Base> const& tree_nodes) {
  if (idx > kPivotNum) {
    assert(idx - kPivotNum - 1 < kBucketNum);
    return tree_nodes[idx - kPivotNum - 1];
  }

  Base const L = BuildInnerTree<Leaf, Interior>(idx << 1, pivots, tree_nodes);
  Base const R =
      BuildInnerTree<Leaf, Interior>(idx << 1 | 1, pivots, tree_nodes);

  if constexpr (std::same_as<HyperPlane, typename Interior::ST>) {
    return AllocInteriorNode<Interior>(L, R, pivots[idx],
                                       typename Interior::AT());
  } else if constexpr (std::same_as<Box, typename Interior::ST>) {
    return AllocInteriorNode<Interior>(
        L, R, GetBox(GetSplit<Leaf, Interior>(L), GetSplit<Leaf, Interior>(R)),
        typename Interior::AT());
  } else {
    static_assert(0);
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsMultiNode Interior>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::BuildInnerTree(
    BucketType idx, HyperPlaneSeq const& pivots,
    parlay::sequence<Node*> const& tree_nodes) {
  assert(idx < kPivotNum + kBucketNum + 1);

  if (idx > kPivotNum) {
    return tree_nodes[idx - kPivotNum - 1];
  }

  typename DerivedTree::OrthNodeArr multi_nodes;
  typename DerivedTree::Splitter split;
  for (DimsType i = 0; i < DerivedTree::kNodeRegions; ++i) {
    multi_nodes[i] = BuildInnerTree<Interior>(
        idx * DerivedTree::kNodeRegions + i, pivots, tree_nodes);
  }
  for (DimsType i = 0; i < DerivedTree::kSplitterNum; ++i) {
    split[i] = pivots[idx * (1 << i)];
    assert(i == 0 || pivots[idx * (1 << i)] == pivots[idx * (1 << i) + 1]);
  }

  return AllocInteriorNode<Interior>(multi_nodes, split,
                                     typename Interior::AT());
}
}  // namespace pstp

#endif  // PSTP_BASE_TREE_IMPL_TREE_OP_BUILD_INNER_TREE_HPP_
