
#pragma once

#include <numeric>
#include "../base_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point, uint8_t kBDO>
template<typename Interior>
Node* BaseTree<Point, kBDO>::BuildInnerTree(
    BucketType idx, HyperPlaneSeq& pivots,
    parlay::sequence<Node*>& tree_nodes) {
    if (idx > kPivotNum) {
        assert(idx - kPivotNum - 1 < kBucketNum);
        return tree_nodes[idx - kPivotNum - 1];
    }
    Node *L, *R;
    L = BuildInnerTree<Interior>(idx << 1, pivots, tree_nodes);
    R = BuildInnerTree<Interior>(idx << 1 | 1, pivots, tree_nodes);
    return AllocInteriorNode<Interior>(L, R, pivots[idx],
                                       typename Interior::AT());
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, IsBinaryNode Interior, typename Range, bool granularity>
void BaseTree<Point, kBDO>::FlattenRec(Node* T, Range Out) {
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
            FlattenRec<Leaf, Interior>(TI->left, Out.cut(0, TI->left->size));
        },
        [&]() {
            FlattenRec<Leaf, Interior>(TI->right,
                                       Out.cut(TI->left->size, TI->size));
        });

    return;
}

template<typename Point, uint8_t kBDO>
template<typename Leaf, IsMultiNode Interior, typename Range, bool granularity>
void BaseTree<Point, kBDO>::FlattenRec(Node* T, Range Out) {
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

    assert(TI->size ==
           std::accumulate(
               TI->tree_nodes.begin(), TI->tree_nodes.end(),
               static_cast<size_t>(0),
               [](size_t acc, Node* n) -> size_t { return acc + n->size; }));

    if ((granularity && TI->size > kSerialBuildCutoff) ||
        (!granularity && TI->ForceParallel())) {
        parlay::parallel_for(0, TI->tree_nodes.size(), [&](BucketType i) {
            size_t start = 0;
            for (BucketType j = 0; j < i; ++j) {
                start += TI->tree_nodes[j]->size;
            }
            FlattenRec<Leaf, Interior, Range>(
                TI->tree_nodes[i], Out.cut(start, TI->tree_nodes[i]->size));
        });
    } else {
        size_t start = 0;
        for (BucketType i = 0; i < TI->tree_nodes.size(); ++i) {
            FlattenRec<Leaf, Interior, Range>(
                TI->tree_nodes[i], Out.cut(start, TI->tree_nodes[i]->size));
            start += TI->tree_nodes[i]->size;
        }
    }

    return;
}
}  // namespace cpdd
