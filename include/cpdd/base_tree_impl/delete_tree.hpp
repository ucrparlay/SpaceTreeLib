#pragma once

#include "../base_tree.h"

namespace cpdd {

template <typename Point, typename DerivedTree, uint_fast8_t kBDO>
template <typename Leaf, typename Interior>
void BaseTree<Point, DerivedTree, kBDO>::DeleteTreeWrapper() {
  if (this->root_ == nullptr) {
    return;
  }
  DeleteTreeRecursive<Leaf, Interior>(this->root_);
  this->root_ = nullptr;
  return;
}

template <typename Point, typename DerivedTree,
          uint_fast8_t kBDO>  //* delete tree in parallel
template <typename Leaf, IsBinaryNode Interior, bool granularity>
void BaseTree<Point, DerivedTree, kBDO>::DeleteTreeRecursive(Node* T) {
  if (T == nullptr) {
    std::cout << "empty ptr" << std::endl;
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

template <typename Point, typename DerivedTree,
          uint_fast8_t kBDO>  //* delete tree in parallel
template <typename Leaf, IsMultiNode Interior, bool granularity>
void BaseTree<Point, DerivedTree, kBDO>::DeleteTreeRecursive(Node* T) {
  if (T == nullptr) return;
  if (T->is_leaf) {
    FreeNode<Leaf>(T);
  } else {
    Interior* TI = static_cast<Interior*>(T);

    // NOTE: enable granularity control by default, if it is disabled,
    // always delete in parallel
    if (ForceParallelRecursion<Interior, granularity>(TI)) {
      parlay::parallel_for(
          0, TI->tree_nodes.size(),
          [&](size_t i) {
            DeleteTreeRecursive<Leaf, Interior>(TI->tree_nodes[i]);
          },
          1);
    } else {
      for (size_t i = 0; i < TI->tree_nodes.size(); ++i) {
        DeleteTreeRecursive<Leaf, Interior>(TI->tree_nodes[i]);
      }
    }
    FreeNode<Interior>(T);
  }
}
}  // namespace cpdd
