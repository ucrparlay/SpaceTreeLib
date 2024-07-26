
#pragma once

#include "../base_tree.h"

namespace cpdd {
template<typename Point>
template<typename SplitType, typename AugType>
Node* BaseTree<Point>::BuildInnerTree(BucketType idx, SplitterSeq& pivots,
                                      parlay::sequence<Node*>& tree_nodes,
                                      const AugType& aug) {
    if (idx > kPivotNum) {
        assert(idx - kPivotNum - 1 < kBucketNum);
        return tree_nodes[idx - kPivotNum - 1];
    }
    Node *L, *R;
    L = BuildInnerTree<SplitType, AugType>(idx << 1, pivots, tree_nodes, aug);
    R = BuildInnerTree<SplitType, AugType>(idx << 1 | 1, pivots, tree_nodes,
                                           aug);
    return AllocInteriorNode<Point, SplitType, AugType>(L, R, pivots[idx], aug);
}
}  // namespace cpdd
