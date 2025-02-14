#ifndef PSPT_BASE_TREE_IMPL_TREE_OP_FLATTEN_HPP_
#define PSPT_BASE_TREE_IMPL_TREE_OP_FLATTEN_HPP_

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
template <SupportsForceParallel Interior, bool granularity>
inline bool BaseTree<Point, DerivedTree, kSkHeight,
                     kImbaRatio>::ForceParallelRecursion(Interior const* TI) {
  return (granularity && TI->size > kSerialBuildCutoff) ||
         (!granularity && TI->ForceParallel());
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsBinaryNode Interior, typename Range,
          bool granularity>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::FlattenRec(
    Node* T, Range Out) {
  assert(T->size == Out.size());

  if (T->size == 0) return;

  if (T->is_leaf) {
    ExtractPointsInLeaf<Leaf>(T, Out);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);
  assert(TI->size == TI->left->size + TI->right->size);
  parlay::par_do_if(
      // WARN: check parallelisim using node size can be biased
      ForceParallelRecursion<Interior, granularity>(TI),
      [&]() {
        FlattenRec<Leaf, Interior>(TI->left, Out.cut(0, TI->left->size));
      },
      [&]() {
        FlattenRec<Leaf, Interior>(TI->right,
                                   Out.cut(TI->left->size, TI->size));
      });

  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsMultiNode Interior, typename Range, bool granularity>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::FlattenRec(
    Node* T, Range Out) {
  assert(T->size == Out.size());

  if (T->size == 0) return;

  if (T->is_leaf) {
    ExtractPointsInLeaf<Leaf>(T, Out);
    return;
  }

  Interior* TI = static_cast<Interior*>(T);

  assert(TI->size ==
         std::accumulate(
             TI->tree_nodes.begin(), TI->tree_nodes.end(),
             static_cast<size_t>(0),
             [](size_t acc, Node* n) -> size_t { return acc + n->size; }));

  parlay::parallel_for(
      0, TI->tree_nodes.size(),
      [&](BucketType i) {
        size_t start = 0;
        for (BucketType j = 0; j < i; ++j) {
          start += TI->tree_nodes[j]->size;
        }
        FlattenRec<Leaf, Interior, Range>(
            TI->tree_nodes[i], Out.cut(start, start + TI->tree_nodes[i]->size));
      },
      ForceParallelRecursion<Interior, granularity>(TI)
          ? 1
          : Interior::GetRegions());

  return;
}

// NOTE: for multi node @T, it only flatten the subtree with id @idx to @Out
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, IsMultiNode Interior, typename Range, bool granularity>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::PartialFlatten(
    Node* T, Range Out, BucketType idx) {
  if (idx == 1) {
    assert(T->size == Out.size());
    FlattenRec<Leaf, Interior>(T, Out.cut(0, T->size));
    return;
  } else if (idx >= Interior::GetRegions()) {
    Node* ns =
        static_cast<Interior*>(T)->tree_nodes[idx - Interior::GetRegions()];
    assert(ns->size == Out.size());
    FlattenRec<Leaf, Interior>(ns, Out.cut(0, ns->size));
    return;
  }

  Interior* TI = static_cast<Interior*>(T);
  size_t l_size = TI->MergeSize(idx << 1), r_size = TI->MergeSize(idx << 1 | 1);
  assert(l_size + r_size == Out.size());
  parlay::par_do_if(
      ForceParallelRecursion<Interior, granularity>(static_cast<Interior*>(T)),
      [&]() {
        PartialFlatten<Leaf, Interior>(T, Out.cut(0, l_size), idx << 1);
      },
      [&]() {
        PartialFlatten<Leaf, Interior>(T, Out.cut(l_size, l_size + r_size),
                                       idx << 1 | 1);
      });
  return;
}
}  // namespace pspt

#endif  // PSPT_BASE_TREE_IMPL_TREE_OP_FLATTEN_HPP_
