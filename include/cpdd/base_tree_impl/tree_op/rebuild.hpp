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

// NOTE: rebuild the tree
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior, bool granularity>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::PrepareRebuild(
    Node* T, Points& wx, Points& wo) {
  wo = Points::uninitialized(T->size);
  wx = Points::uninitialized(T->size);
  FlattenRec<Leaf, Interior, Slice, granularity>(T, parlay::make_slice(wx));
  DeleteTreeRecursive<Leaf, Interior, granularity>(T);
  return;
}

// NOTE: rebuild with new input In
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior, bool granularity>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::PrepareRebuild(
    Node* T, Slice In, Points& wx, Points& wo) {
  wo = Points::uninitialized(T->size + In.size());
  wx = Points::uninitialized(T->size + In.size());
  parlay::parallel_for(0, In.size(), [&](size_t j) { wx[j] = In[j]; });
  FlattenRec<Leaf, Interior, Slice, granularity>(T,
                                                 wx.cut(In.size(), wx.size()));
  DeleteTreeRecursive<Leaf, Interior, granularity>(T);
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior, typename... Args>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RebuildWithInsert(
    Node* T, Slice In, Args&&... args) {
  Points wx, wo;
  PrepareRebuild<Leaf, Interior>(T, In, wx, wo);
  static_assert(
      std::is_invocable_v<decltype(&DerivedTree::BuildRecursive), DerivedTree*,
                          Slice, Slice, Args&&..., Box>);
  return static_cast<DerivedTree*>(this)->BuildRecursive(
      parlay::make_slice(wx), parlay::make_slice(wo),
      std::forward<Args>(args)..., GetBox(parlay::make_slice(wx)));
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior, bool granularity, typename... Args>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RebuildSingleTree(
    Node* T, Args&&... args) {
  Points wx, wo;
  PrepareRebuild<Leaf, Interior, granularity>(T, wx, wo);
  static_assert(std::is_invocable_v<decltype(&DerivedTree::BuildRecursive),
                                    DerivedTree*, Slice, Slice, Args&&...>);
  return static_cast<DerivedTree*>(this)->BuildRecursive(
      parlay::make_slice(wx), parlay::make_slice(wo),
      std::forward<Args>(args)...);
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsBinaryNode Interior, bool granularity,
          typename PrepareFunc, typename... Args>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RebuildTreeRecursive(
    Node* T, PrepareFunc&& prepare_func, Args&&... args) {
  if (T->is_leaf) {
    return T;
  }

  Interior* TI = static_cast<Interior*>(T);
  if (ImbalanceNode(TI->left->size, TI->size)) {
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
                  TI->left, prepare_func,
                  std::forward<decltype(left_args)>(left_args)...);
            },
            left_args);
      },
      [&] {
        R = std::apply(
            [&](auto&&... right_args) {
              return RebuildTreeRecursive<Leaf, Interior, granularity>(
                  TI->right, prepare_func,
                  std::forward<decltype(right_args)>(right_args)...);
            },
            right_args);
      });

  TI->ResetParallelFlag();
  UpdateInterior<Interior>(T, L, R);
  return T;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsMultiNode Interior, bool granularity,
          typename PrepareFunc, typename... Args>
Node* BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::RebuildTreeRecursive(
    Node* T, PrepareFunc&& prepare_func, Args&&... args) {
  if (T->is_leaf) {
    return T;
  }

  Interior* TI = static_cast<Interior*>(T);
  if (SparcyNode(0,
                 TI->size)) {  // NOTE: after diff points from the tree, all
                               // points has been removed, so there the remove
                               // points are 0, simply use it to check whether
                               // the number of points are below kLeaveWrap
    return RebuildSingleTree<Leaf, Interior, granularity>(
        T, std::forward<Args>(args)...);
  }

  typename Interior::NodeArr new_nodes;
  if (ForceParallelRecursion<Interior, granularity>(TI)) {
    parlay::parallel_for(0, TI->tree_nodes.size(), [&](BucketType i) {
      auto const new_args = prepare_func(T, i, std::forward<Args>(args)...);
      std::apply(
          [&](auto&&... new_args) {
            new_nodes[i] = RebuildTreeRecursive<Leaf, Interior, granularity>(
                TI->tree_nodes[i], prepare_func,
                std::forward<decltype(new_args)>(new_args)...);
          },
          new_args);
    });
  } else {
    for (BucketType i = 0; i < TI->tree_nodes.size(); ++i) {
      auto const new_args = prepare_func(T, i, std::forward<Args>(args)...);
      std::apply(
          [&](auto&&... new_args) {
            new_nodes[i] = RebuildTreeRecursive<Leaf, Interior, granularity>(
                TI->tree_nodes[i], prepare_func,
                std::forward<decltype(new_args)>(new_args)...);
          },
          new_args);
    }
  }

  TI->ResetParallelFlag();
  UpdateInterior<Interior>(T, new_nodes);
  return T;
}

}  // namespace cpdd
