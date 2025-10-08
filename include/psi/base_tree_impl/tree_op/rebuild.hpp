#ifndef PSI_BASE_TREE_IMPL_REBUILD_HPP_
#define PSI_BASE_TREE_IMPL_REBUILD_HPP_

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

// NOTE: rebuild the tree
BASETREE_TEMPLATE
template <typename Leaf, typename Interior, bool granularity>
void BASETREE_CLASS::PrepareRebuild(Node* T, Points& wx, Points& wo) {
  wo = Points::uninitialized(T->size);
  wx = Points::uninitialized(T->size);
  FlattenRec<Leaf, Interior, Slice, granularity>(T, parlay::make_slice(wx));
  DeleteTreeRecursive<Leaf, Interior, granularity>(T);
  return;
}

// NOTE: rebuild with new input In
BASETREE_TEMPLATE
template <typename Leaf, typename Interior, bool granularity>
void BASETREE_CLASS::PrepareRebuild(Node* T, Slice In, Points& wx, Points& wo) {
  wo = Points::uninitialized(T->size + In.size());
  wx = Points::uninitialized(T->size + In.size());
  parlay::parallel_for(0, In.size(), [&](size_t j) { wx[j] = In[j]; });
  FlattenRec<Leaf, Interior, Slice, granularity>(T,
                                                 wx.cut(In.size(), wx.size()));
  DeleteTreeRecursive<Leaf, Interior, granularity>(T);
  return;
}

BASETREE_TEMPLATE
template <typename Leaf, typename Interior, typename PrepareFunc,
          typename... Args>
Node* BASETREE_CLASS::RebuildWithInsert(Node* T, PrepareFunc prepare_func,
                                        Slice In, Args&&... args) {
  Points w_in, w_out;
  PrepareRebuild<Leaf, Interior>(T, In, w_in, w_out);
  auto additional_arg = prepare_func(T, w_in, w_out);
  static_assert(
      std::is_invocable_v<decltype(&DerivedTree::BuildRecursive), DerivedTree*,
                          Slice, Slice, Args&&..., Box>);
  return static_cast<DerivedTree*>(this)->BuildRecursive(
      parlay::make_slice(w_in), parlay::make_slice(w_out),
      std::forward<Args>(args)..., additional_arg);
}

// PARA: when granularity set to false, it will disable the default value for
// granularity size, which is BT::kSerialBuildCutoff; instead, it will check
// whether the force parallel flag has been enabled
BASETREE_TEMPLATE
template <typename Leaf, typename Interior, bool granularity, typename... Args>
Node* BASETREE_CLASS::RebuildSingleTree(Node* T, Args&&... args) {
  Points wx, wo;
  PrepareRebuild<Leaf, Interior, granularity>(T, wx, wo);
  static_assert(std::is_invocable_v<decltype(&DerivedTree::BuildRecursive),
                                    DerivedTree*, Slice, Slice, Args&&...>);
  return static_cast<DerivedTree*>(this)->BuildRecursive(
      parlay::make_slice(wx), parlay::make_slice(wo),
      std::forward<Args>(args)...);
}

// NOTE: traverse a tree, if it satisfy the condition, then rebuild a binary
// tree
// PARA: if allow_enable_rebuild enabled, this method will re-balance the
// tree; otherwise, it flattens all sparcy node into leaf nodes
BASETREE_TEMPLATE
template <typename Leaf, IsBinaryNode Interior, bool granularity,
          typename PrepareFunc, typename... Args>
Node* BASETREE_CLASS::RebuildTreeRecursive(Node* T, PrepareFunc&& prepare_func,
                                           bool const allow_inba_rebuild,
                                           Args&&... args) {
  if (T->is_leaf) {
    return T;
  }

  Interior* TI = static_cast<Interior*>(T);
  // NOTE: rebuild the tree if it is sparcy or imbalance
  if (SparcyNode(0, TI->size) ||
      (allow_inba_rebuild && ImbalanceNode(TI->left->size, TI->size))) {
    return RebuildSingleTree<Leaf, Interior, granularity>(
        T, std::forward<Args>(args)...);
  }

  auto const [left_args, right_args] =
      prepare_func(T, std::forward<Args>(args)...);

  Node *L, *R;
  parlay::par_do_if(
      ForceParallelRecursion<Interior, granularity>(TI),
      [&] {
        L = std::apply(
            [&](auto&&... left_args) {
              return RebuildTreeRecursive<Leaf, Interior, granularity>(
                  TI->left, prepare_func, allow_inba_rebuild,
                  std::forward<decltype(left_args)>(left_args)...);
            },
            left_args);
      },
      [&] {
        R = std::apply(
            [&](auto&&... right_args) {
              return RebuildTreeRecursive<Leaf, Interior, granularity>(
                  TI->right, prepare_func, allow_inba_rebuild,
                  std::forward<decltype(right_args)>(right_args)...);
            },
            right_args);
      });

  UpdateInterior<Interior>(T, L, R);
  return T;
}

// NOTE: rebuild a multi-node tree
BASETREE_TEMPLATE
template <typename Leaf, IsMultiNode Interior, bool granularity,
          typename PrepareFunc, typename... Args>
Node* BASETREE_CLASS::RebuildTreeRecursive(
    Node* T, PrepareFunc&& prepare_func,
    [[maybe_unused]] bool const allow_inba_rebuild, Args&&... args) {
  if (T->is_leaf) {
    return T;
  }

  Interior* TI = static_cast<Interior*>(T);
  if (SparcyNode(0, TI->size)) {
    return RebuildSingleTree<Leaf, Interior, granularity>(
        T, std::forward<Args>(args)...);
  }

  typename Interior::NodeArr new_nodes;
  parlay::parallel_for(
      0, Interior::GetRegions(),
      [&](BucketType i) {
        auto const new_args = prepare_func(T, i, std::forward<Args>(args)...);
        std::apply(
            [&](auto&&... new_args) {
              new_nodes[i] = RebuildTreeRecursive<Leaf, Interior, granularity>(
                  TI->tree_nodes[i], prepare_func, allow_inba_rebuild,
                  std::forward<decltype(new_args)>(new_args)...);
            },
            new_args);
      },
      ForceParallelRecursion<Interior, granularity>(TI)
          ? 1
          : Interior::GetRegions());

  UpdateInterior<Interior>(T, new_nodes);
  return T;
}

}  // namespace psi

#undef BASETREE_TEMPLATE
#undef BASETREE_CLASS

#endif  // PSI_BASE_TREE_IMPL_REBUILD_HPP_
