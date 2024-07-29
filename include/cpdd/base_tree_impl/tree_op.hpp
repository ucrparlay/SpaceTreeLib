
#pragma once

#include "../base_tree.h"

namespace cpdd {
template<typename Point, uint8_t kBDO>
template<typename Interior, typename SplitterSeq>
Node* BaseTree<Point, kBDO>::BuildInnerTree(
    BucketType idx, SplitterSeq& pivots, parlay::sequence<Node*>& tree_nodes) {
    if (idx > kPivotNum) {
        assert(idx - kPivotNum - 1 < kBucketNum);
        return tree_nodes[idx - kPivotNum - 1];
    }
    Node *L, *R;
    L = BuildInnerTree<Interior, SplitterSeq>(idx << 1, pivots, tree_nodes);
    R = BuildInnerTree<Interior, SplitterSeq>(idx << 1 | 1, pivots, tree_nodes);
    return AllocInteriorNode<Interior>(L, R, pivots[idx],
                                       typename Interior::AT());
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, typename Interior, typename Range>
void BaseTree<Point, kBDO>::FlattenRec(Node* T, Range Out, bool granularity) {
    assert(T->size == Out.size());

    if (T->size == 0) return;

    if (T->is_leaf) {
        Leaf* TL = static_cast<Leaf*>(T);
        for (int i = 0; i < TL->size; i++) {
            Out[i] = TL->pts[(!TL->is_dummy) * i];
        }
        return;
    }

    Interior* TI = static_cast<Interior*>(T);
    assert(TI->size == TI->left->size + TI->right->size);
    parlay::par_do_if(
        // WARN: check parallelisim using node size can be biased
        (granularity && TI->size > kSerialBuildCutoff) ||
            (!granularity && TI->ForceParallel()),
        [&]() {
            FlattenRec<Leaf, Interior, Range>(
                TI->left, Out.cut(0, TI->left->size), granularity);
        },
        [&]() {
            FlattenRec<Leaf, Interior, Range>(
                TI->right, Out.cut(TI->left->size, TI->size), granularity);
        });

    return;
}
}  // namespace cpdd
