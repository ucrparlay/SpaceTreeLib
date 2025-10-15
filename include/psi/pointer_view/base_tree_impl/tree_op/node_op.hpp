#ifndef PSI_POINTER_VIEW_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_
#define PSI_POINTER_VIEW_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_

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
template <typename Leaf, typename Interior>
auto BaseTree<TypeTrait, DerivedTree>::RetriveBox(Node const* node)
  requires(HasBox<typename Leaf::AT> && HasBox<typename Interior::AT>)
{
  return node->is_leaf ? static_cast<Leaf const*>(node)->GetBox()
                       : static_cast<Interior const*>(node)->GetBox();
}

template <class TypeTrait, typename DerivedTree>
template <typename Leaf, typename Interior>
typename Interior::ST const& BaseTree<TypeTrait, DerivedTree>::GetSplit(
    Node const* node)
  requires std::same_as<typename BaseTree<TypeTrait, DerivedTree>::Box,
                        typename Interior::ST>
{
  return node->is_leaf ? static_cast<Leaf const*>(node)->GetSplit()
                       : static_cast<Interior const*>(node)->GetSplit();
}

// NOTE: update the info of T by new children L and R
template <class TypeTrait, typename DerivedTree>
template <IsBinaryNode Interior, bool UpdateParFlag>
inline void BaseTree<TypeTrait, DerivedTree>::UpdateInterior(Node* T, Node* L,
                                                             Node* R) {
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

template <class TypeTrait, typename DerivedTree>
template <IsBinaryNode Interior, bool UpdateParFlag>
inline void BaseTree<TypeTrait, DerivedTree>::UpdateInterior(Node* T,
                                                             NodeBox const& L,
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

template <class TypeTrait, typename DerivedTree>
template <IsMultiNode Interior, bool UpdateParFlag>
inline void BaseTree<TypeTrait, DerivedTree>::UpdateInterior(
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
}  // namespace pointer_view
}  // namespace psi

#endif  // PSI_POINTER_VIEW_BASE_TREE_IMPL_TREE_OP_NODE_OP_HPP_
