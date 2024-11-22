#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#include "../../base_tree.h"

namespace cpdd {

// NOTE: expand a multi node into a binary node recursively
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsBinaryNode BN, IsMultiNode MN>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::ExpandMultiNode(
    typename MN::ST const& split, BucketType idx, BucketType deep,
    parlay::sequence<Node*> const& tree_nodes) {
  if (idx >= MN::kRegions) {
    return tree_nodes[idx - MN::kRegions];
  }
  Node *L, *R;
  L = ExpandMultiNode<BN, MN>(split, idx * 2, deep + 1, tree_nodes);
  R = ExpandMultiNode<BN, MN>(split, idx * 2 + 1, deep + 1, tree_nodes);
  auto o = AllocInteriorNode<BN>(L, R, split[deep], typename BN::AT());
  assert(o->size == L->size + R->size);
  return o;
}

// NOTE: expand a tree with multi-way nodes into one with binary nodes
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsBinaryNode BN, IsMultiNode MN>
  requires std::same_as<typename BN::ST, typename MN::ST::value_type>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Expand2Binary(
    Node* T) {
  if (T->is_leaf) {
    return T;
  }
  MN* MI = static_cast<MN*>(T);
  parlay::sequence<Node*> tree_nodes(MN::kRegions);
  assert(MI->tree_nodes.size() == MN::kRegions);
  for (BucketType i = 0; i < MN::kRegions; ++i) {
    tree_nodes[i] = Expand2Binary<BN, MN>(MI->tree_nodes[i]);
  }
  assert(std::accumulate(tree_nodes.begin(), tree_nodes.end(), 0,
                         [](size_t acc, Node* n) -> size_t {
                           return acc + n->size;
                         }) == MI->size);
  auto split = MI->split;
  return ExpandMultiNode<BN, MN>(split, 1, 0, tree_nodes);
}
}  // namespace cpdd
