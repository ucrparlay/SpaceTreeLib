#ifndef PSPT_BASE_TREE_IMPL_DELETE_TREE_HPP_
#define PSPT_BASE_TREE_IMPL_DELETE_TREE_HPP_

#include "../base_tree.h"
#include "dependence/concepts.h"

namespace pspt {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::DeleteTreeWrapper() {
  if (this->root_ == nullptr) {
    return;
  }
  DeleteTreeRecursive<Leaf, Interior>(this->root_);
  this->tree_box_ = GetEmptyBox();
  this->root_ = nullptr;
  return;
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>  //* delete tree in parallel
template <typename Leaf, IsBinaryNode Interior, bool granularity>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::DeleteTreeRecursive(
    Node* T) {
  if (T == nullptr) {
    return;
  }
  if (T->is_leaf) {
    FreeNode<Leaf>(T);
  } else {
    Interior* TI = static_cast<Interior*>(T);
    // NOTE: enable granularity control by default, if it is disabled,
    // always delete in parallel
    parlay::par_do_if(
        ForceParallelRecursion<Interior, granularity>(TI),
        [&] { DeleteTreeRecursive<Leaf, Interior>(TI->left); },
        [&] { DeleteTreeRecursive<Leaf, Interior>(TI->right); });
    FreeNode<Interior>(T);
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>  //* delete tree in parallel
template <typename Leaf, IsMultiNode Interior, bool granularity>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::DeleteTreeRecursive(
    Node* T) {
  if (T == nullptr) {
    return;
  }
  if (T->is_leaf) {
    FreeNode<Leaf>(T);
  } else {
    Interior* TI = static_cast<Interior*>(T);

    // NOTE: enable granularity control by default, if it is disabled,
    // always delete in parallel
    parlay::parallel_for(
        0, TI->tree_nodes.size(),
        [&](size_t i) {
          DeleteTreeRecursive<Leaf, Interior>(TI->tree_nodes[i]);
        },
        ForceParallelRecursion<Interior, granularity>(TI)
            ? 1
            : Interior::GetRegions());

    FreeNode<Interior>(T);
  }
}

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>  //* delete tree in parallel
template <typename Leaf, IsDynamicNode Interior, bool granularity>
void BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::DeleteTreeRecursive(
    Node* T) {
  if (T == nullptr) {
    return;
  }
  if (T->is_leaf) {
    FreeNode<Leaf>(T);
  } else {
    Interior* TI = static_cast<Interior*>(T);

    // NOTE: enable granularity control by default, if it is disabled,
    // always delete in parallel
    parlay::parallel_for(
        0, TI->tree_nodes.size(),
        [&](size_t i) {
          DeleteTreeRecursive<Leaf, Interior>(TI->tree_nodes[i]);
        },
        ForceParallelRecursion<Interior, granularity>(TI)
            ? 1
            : TI->tree_nodes.size());

    FreeNode<Interior>(T);
  }
}
}  // namespace pspt

#endif  // PSPT_BASE_TREE_IMPL_DELETE_TREE_HPP_
