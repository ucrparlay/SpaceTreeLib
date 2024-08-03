#include "../base_tree.h"
#include "cpdd/dependence/tree_node.h"
#include "parlay/parallel.h"

namespace cpdd {

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior, typename MultiWayTag>
void BaseTree<Point, kBDO>::DeleteTreeWrapper() {
    if (this->root_ == nullptr) {
        return;
    }
    DeleteTreeRecursive<Leaf, Interior>(MultiWayTag(), this->root_);
    this->root_ = nullptr;
    return;
}

template<typename Point, uint8_t kBDO>  //* delete tree in parallel
template<typename Leaf, typename Interior>
void BaseTree<Point, kBDO>::DeleteTreeRecursive(BinaryInteriorTag, Node* T,
                                                bool granularity) {
    if (T == nullptr) {
        LOG << "empty ptr" << ENDL;
        return;
    }
    if (T->is_leaf) {
        FreeNode<Leaf>(T);
    } else {
        Interior* TI = static_cast<Interior*>(T);
        // NOTE: enable granularity control by default, if it is disabled,
        // always delete in parallel
        parlay::par_do_if((granularity && T->size > kSerialBuildCutoff) ||
                              (!granularity && TI->ForceParallel()),
                          [&] {
                              DeleteTreeRecursive<Leaf, Interior>(
                                  BinaryInteriorTag(), TI->left, granularity);
                          },
                          [&] {
                              DeleteTreeRecursive<Leaf, Interior>(
                                  BinaryInteriorTag(), TI->right, granularity);
                          });
        FreeNode<Interior>(T);
    }
}

template<typename Point, uint8_t kBDO>  //* delete tree in parallel
template<typename Leaf, typename Interior>
void BaseTree<Point, kBDO>::DeleteTreeRecursive(MultiWayInteriorTag, Node* T,
                                                bool granularity) {
    if (T == nullptr) return;
    if (T->is_leaf) {
        FreeNode<Leaf>(T);
    } else {
        Interior* TI = static_cast<Interior*>(T);

        // NOTE: enable granularity control by default, if it is disabled,
        // always delete in parallel
        if ((granularity && T->size > kSerialBuildCutoff) ||
            (!granularity && TI->ForceParallel())) {
            parlay::parallel_for(
                0, TI->tree_nodes.size(),
                [&](size_t i) {
                    DeleteTreeRecursive<Leaf, Interior>(
                        MultiWayInteriorTag(), TI->tree_nodes[i], granularity);
                },
                1);
        } else {
            for (size_t i = 0; i < TI->tree_nodes.size(); ++i) {
                DeleteTreeRecursive<Leaf, Interior>(
                    MultiWayInteriorTag(), TI->tree_nodes[i], granularity);
            }
        }
        FreeNode<Interior>(T);
    }
}
}  // namespace cpdd
