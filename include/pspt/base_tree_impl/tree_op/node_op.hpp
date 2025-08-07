#ifndef PSPT_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_
#define PSPT_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#include "../../base_tree.h"

namespace pspt {
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
static Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::MergeUp(
    Node* T) {
  if (T->is_leaf) {
    auto TL = static_cast<Leaf*>(T);
    TL->aug = GetBox(static_cast<Leaf*>(T)->pts.cut(0, T->size));
    return T;
  }
  auto TI = static_cast<Interior*>(T);
  auto left = MergeUp<Leaf, Interior>(TI->left);
  auto right = MergeUp<Leaf, Interior>(TI->right);
  auto left_box = left->is_leaf ? static_cast<Leaf*>(left)->aug
                                : static_cast<Interior*>(left)->aug.box;
  auto right_box = right->is_leaf ? static_cast<Leaf*>(right)->aug
                                  : static_cast<Interior*>(right)->aug.box;
  TI->aug.box = GetBox(left_box, right_box);
  return TI;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
typename Interior::ST const&
BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::GetSplit(Node const* node)
  requires std::same_as<
      typename BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::Box,
      typename Interior::ST>
{
  return node->is_leaf ? static_cast<Leaf const*>(node)->GetSplit()
                       : static_cast<Interior const*>(node)->GetSplit();
}

// NOTE: update the info of T by new children L and R
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsBinaryNode Interior, bool UpdateParFlag>
inline void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::UpdateInterior(
    Node* T, Node* L, Node* R) {
  assert(!T->is_leaf);
  Interior* TI = static_cast<Interior*>(T);
  if constexpr (UpdateParFlag) {
    TI->ResetParallelFlag();
  }
  TI->size = L->size + R->size;
  TI->left = L;
  TI->right = R;
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsBinaryNode Interior, bool UpdateParFlag>
inline void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::UpdateInterior(
    Node* T, NodeBox const& L, NodeBox const& R) {
  assert(!T->is_leaf);
  Interior* TI = static_cast<Interior*>(T);
  if constexpr (UpdateParFlag) {
    TI->ResetParallelFlag();
  }
  TI->size = L.first->size + R.first->size;
  TI->left = L.first;
  TI->right = R.first;
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <IsMultiNode Interior, bool UpdateParFlag>
inline void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::UpdateInterior(
    Node* T, typename Interior::NodeArr const& new_nodes) {
  assert(!T->is_leaf);
  Interior* TI = static_cast<Interior*>(T);
  if constexpr (UpdateParFlag) {
    TI->ResetParallelFlag();
  }
  TI->size = std::accumulate(
      new_nodes.begin(), new_nodes.end(), 0,
      [](size_t acc, Node* n) -> size_t { return acc + n->size; });
  TI->tree_nodes = new_nodes;
  return;
}
}  // namespace pspt

#endif  // PSPT_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_
