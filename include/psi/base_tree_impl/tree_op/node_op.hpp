#ifndef PSI_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_
#define PSI_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#include "../../base_tree.h"

#define BASETREE_TEMPLATE                                                 \
  template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight, \
            uint_fast8_t kImbaRatio>
#define BASETREE_CLASS BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>

namespace psi {
BASETREE_TEMPLATE
template <typename Leaf, typename Interior>
auto BASETREE_CLASS::RetriveBox(Node const* node)
  requires(HasBox<typename Leaf::AT> && HasBox<typename Interior::AT>)
{
  return node->is_leaf ? static_cast<Leaf const*>(node)->GetBox()
                       : static_cast<Interior const*>(node)->GetBox();
}

BASETREE_TEMPLATE
template <typename Leaf, typename Interior>
typename Interior::ST const& BASETREE_CLASS::GetSplit(Node const* node)
  requires std::same_as<typename BASETREE_CLASS::Box, typename Interior::ST>
{
  return node->is_leaf ? static_cast<Leaf const*>(node)->GetSplit()
                       : static_cast<Interior const*>(node)->GetSplit();
}

// NOTE: update the info of T by new children L and R
BASETREE_TEMPLATE
template <IsBinaryNode Interior, bool UpdateParFlag>
inline void BASETREE_CLASS::UpdateInterior(Node* T, Node* L, Node* R) {
  assert(!T->is_leaf);
  Interior* TI = static_cast<Interior*>(T);
  if constexpr (UpdateParFlag) {
    TI->ResetParallelFlag();
  }
  TI->size = L->size + R->size;
  TI->left = L;
  TI->right = R;
  TI->UpdateAug(L, R);
  return;
}

BASETREE_TEMPLATE
template <IsBinaryNode Interior, bool UpdateParFlag>
inline void BASETREE_CLASS::UpdateInterior(Node* T, NodeBox const& L,
                                           NodeBox const& R) {
  assert(!T->is_leaf);
  Interior* TI = static_cast<Interior*>(T);
  if constexpr (UpdateParFlag) {
    TI->ResetParallelFlag();
  }
  TI->size = L.first->size + R.first->size;
  TI->left = L.first;
  TI->right = R.first;
  TI->UpdateAug(L.first, R.first);
  return;
}

BASETREE_TEMPLATE
template <IsMultiNode Interior, bool UpdateParFlag>
inline void BASETREE_CLASS::UpdateInterior(
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
  TI->UpdateAug(new_nodes);
  return;
}
}  // namespace psi

#undef BASETREE_TEMPLATE
#undef BASETREE_CLASS

#endif  // PSI_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_
