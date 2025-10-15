#ifndef PSI_POINTER_VIEW_BASE_TREE_IMPL_TREE_OP_COMPRESS_EXPAND_TREE_HPP_
#define PSI_POINTER_VIEW_BASE_TREE_IMPL_TREE_OP_COMPRESS_EXPAND_TREE_HPP_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>

#include "dependence/tree_node.h"
#include "pointer_view/base_tree.h"

namespace psi {
namespace pointer_view {

// NOTE: expand a multi node into a binary node recursively
template <class TypeTrait, typename DerivedTree>
template <IsBinaryNode BN, IsMultiNode MN>
Node* BaseTree<TypeTrait, DerivedTree>::ExpandMultiNode(
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
template <class TypeTrait, typename DerivedTree>
template <IsBinaryNode BN, IsMultiNode MN>
Node* BaseTree<TypeTrait, DerivedTree>::Expand2Binary(Node* T)
  requires std::same_as<typename BN::ST, typename MN::ST::value_type>
{
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
  FreeNode<MN>(T);
  return ExpandMultiNode<BN, MN>(split, 1, 0, tree_nodes);
}

template <class TypeTrait, typename DerivedTree>
template <IsBinaryNode BN, IsMultiNode MN>
void BaseTree<TypeTrait, DerivedTree>::CompressBinaryNode(
    Node* bn_proto, typename MN::NodeArr& tree_nodes, typename MN::ST& split,
    BucketType idx) {
  if (idx >= MN::GetRegions()) {
    tree_nodes[idx - MN::GetRegions()] = bn_proto;
    return;
  }
  assert(!bn_proto->is_leaf);
  auto bn = static_cast<BN*>(bn_proto);
  split[idx] = bn->split;
  CompressBinaryNode<BN, MN>(bn->left, tree_nodes, split, idx * 2);
  CompressBinaryNode<BN, MN>(bn->right, tree_nodes, split, idx * 2 + 1);
  FreeNode<BN>(bn_proto);
  return;
}

template <class TypeTrait, typename DerivedTree>
template <IsBinaryNode BN, IsMultiNode MN>
Node* BaseTree<TypeTrait, DerivedTree>::Compress2Multi(Node* T)
  requires std::same_as<typename BN::ST, typename MN::ST::value_type>
{
  if (T->is_leaf) {
    return T;
  }

  auto bn = static_cast<BN*>(T);
  if (BN::template TestDepth<BN>(T, 0, MN::GetLevels())) {
    typename MN::NodeArr tree_nodes;
    typename MN::ST split;
    CompressBinaryNode<BN, MN>(T, tree_nodes, split, 1);
    for (BucketType i = 0; i < MN::GetRegions(); ++i) {
      tree_nodes[i] = Compress2Multi<BN, MN>(tree_nodes[i]);
    }
    auto o = AllocInteriorNode<MN>(tree_nodes, split, typename MN::AT());
    static_cast<MN*>(o)->SetParallelFlag(true);
    return o;
  } else {
    Node* L = Compress2Multi<BN, MN>(bn->left);
    Node* R = Compress2Multi<BN, MN>(bn->right);
    auto o = AllocInteriorNode<BN>(L, R, bn->split, typename BN::AT());
    static_cast<BN*>(o)->SetParallelFlag(false);
    assert(o->size == L->size + R->size);
    FreeNode<BN>(T);
    return o;
  }
}

}  // namespace pointer_view
}  // namespace psi

#endif  // PSI_POINTER_VIEW_BASE_TREE_IMPL_TREE_OP_COMPRESS_EXPAND_TREE_HPP_
