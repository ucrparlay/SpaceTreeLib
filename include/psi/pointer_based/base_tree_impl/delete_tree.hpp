#ifndef PSI_BASE_TREE_IMPL_DELETE_TREE_HPP_
#define PSI_BASE_TREE_IMPL_DELETE_TREE_HPP_

#include "../base_tree.h"
#include "../../dependence/concepts.h"

#define BASETREE_TEMPLATE                                                 \
  template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight, \
            uint_fast8_t kImbaRatio>
#define BASETREE_CLASS BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>

namespace psi {

BASETREE_TEMPLATE
template <typename Leaf, typename Interior>
void BASETREE_CLASS::DeleteTreeWrapper() {
  if (this->root_ == nullptr) {
    return;
  }
  DeleteTreeRecursive<Leaf, Interior>(this->root_);
  this->tree_box_ = GetEmptyBox();
  this->root_ = nullptr;
  return;
}

BASETREE_TEMPLATE  //* delete tree in parallel
    template <typename Leaf, IsBinaryNode Interior, bool granularity>
    void BASETREE_CLASS::DeleteTreeRecursive(Node* T) {
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

BASETREE_TEMPLATE  //* delete tree in parallel
    template <typename Leaf, IsMultiNode Interior, bool granularity>
    void BASETREE_CLASS::DeleteTreeRecursive(Node* T) {
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

BASETREE_TEMPLATE  //* delete tree in parallel
    template <typename Leaf, IsDynamicNode Interior, bool granularity>
    void BASETREE_CLASS::DeleteTreeRecursive(Node* T) {
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
}  // namespace psi

#undef BASETREE_TEMPLATE
#undef BASETREE_CLASS

#endif  // PSI_BASE_TREE_IMPL_DELETE_TREE_HPP_
