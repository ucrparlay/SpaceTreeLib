
#pragma once

#include "../base_tree.h"

namespace cpdd {
template<typename Point>
template<typename Interior, typename SplitterSeq>
Node* BaseTree<Point>::BuildInnerTree(BucketType idx, SplitterSeq& pivots,
                                      parlay::sequence<Node*>& tree_nodes) {
    if (idx > kPivotNum) {
        assert(idx - kPivotNum - 1 < kBucketNum);
        return tree_nodes[idx - kPivotNum - 1];
    }
    Node *L, *R;
    L = BuildInnerTree<Interior, SplitterSeq>(idx << 1, pivots, tree_nodes);
    R = BuildInnerTree<Interior, SplitterSeq>(idx << 1 | 1, pivots, tree_nodes);
    return AllocInteriorNode<Point, typename Interior::ST,
                             typename Interior::AT>(L, R, pivots[idx],
                                                    typename Interior::AT());
}

template<typename Point>
template<typename Leaf, typename Interior, typename Range, typename UnaryPred>
void BaseTree<Point>::FlattenRec(Node* T, Range Out, UnaryPred&& F,
                                 bool granularity) {
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
            (!granularity && F(TI)),
        [&]() {
            FlattenRec<Leaf, Interior, Range, UnaryPred>(
                TI->left, Out.cut(0, TI->left->size),
                std::forward<UnaryPred>(F), granularity);
        },
        [&]() {
            FlattenRec<Leaf, Interior, Range, UnaryPred>(
                TI->right, Out.cut(TI->left->size, TI->size),
                std::forward<UnaryPred>(F), granularity);
        });

    return;
}
}  // namespace cpdd
