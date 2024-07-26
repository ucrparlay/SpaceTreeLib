#include "../base_tree.h"

namespace cpdd {

template<typename Point>
template<typename leaf_node_type, typename interior_node_type>
void BaseTree<Point>::DeleteTreeWrapper() {
    if (this->root_ == nullptr) {
        return;
    }
    DeleteTreeRecursive<leaf_node_type, interior_node_type>(this->root_);
    this->root_ = nullptr;
    return;
}

template<typename Point>  //* delete tree in parallel
template<typename leaf_node_type, typename interior_node_type>
void BaseTree<Point>::DeleteTreeRecursive(Node* T, bool granularity) {
    if (T == nullptr) return;
    if (T->is_leaf) {
        FreeNode<leaf_node_type>(T);
    } else {
        interior_node_type* TI = static_cast<interior_node_type*>(T);
        // NOTE: enable granularity control by default, if it is disabled,
        // always delete in parallel
        parlay::par_do_if(
            (granularity && T->size > kSerialBuildCutoff) ||
                (!granularity && TI->aug),
            [&] {
                DeleteTreeRecursive<leaf_node_type, interior_node_type>(
                    TI->left, granularity);
            },
            [&] {
                DeleteTreeRecursive<leaf_node_type, interior_node_type>(
                    TI->right, granularity);
            });
        FreeNode<interior_node_type>(T);
    }
}

}  // namespace cpdd
