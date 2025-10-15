#ifndef PSI_POINTER_VIEW_BASE_TREE_IMPL_TREE_OP_BUILD_INNER_TREE_HPP_
#define PSI_POINTER_VIEW_BASE_TREE_IMPL_TREE_OP_BUILD_INNER_TREE_HPP_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#include "../../base_tree.h"

namespace psi {
namespace pointer_view {
template <class TypeTrait, typename DerivedTree>
template <typename Leaf, IsBinaryNode Interior, typename ReturnType>
ReturnType BaseTree<TypeTrait, DerivedTree>::BuildInnerTree(
    BucketType idx, HyperPlaneSeq const& pivots,
    parlay::sequence<ReturnType> const& tree_nodes) {
  if (idx > kPivotNum) {
    assert(idx - kPivotNum - 1 < kBucketNum);
    return tree_nodes[idx - kPivotNum - 1];
  }

  ReturnType const L =
      BuildInnerTree<Leaf, Interior>(idx << 1, pivots, tree_nodes);
  ReturnType const R =
      BuildInnerTree<Leaf, Interior>(idx << 1 | 1, pivots, tree_nodes);

  return AllocInteriorNode<Interior>(L, R, pivots[idx]);
}

template <class TypeTrait, typename DerivedTree>
template <IsMultiNode Interior>
Node* BaseTree<TypeTrait, DerivedTree>::BuildInnerTree(
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

  return AllocInteriorNode<Interior>(multi_nodes, split);
}
}  // namespace pointer_view
}  // namespace psi

#endif  // PSI_POINTER_VIEW_BASE_TREE_IMPL_TREE_OP_BUILD_INNER_TREE_HPP_
