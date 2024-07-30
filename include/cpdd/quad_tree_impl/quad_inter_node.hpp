#pragma once

#include "../quad_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint8_t kBDO>
struct QuadTree<Point, SplitRule, kBDO>::QuadInteriorNode :
    InteriorNode<Point, SplitType, AugType> {
    using ST = SplitType;
    using AT = AugType;

    QuadInteriorNode(Node* _left, Node* _right, const ST& _split,
                     const AT& _aug) :
        InteriorNode<Point, SplitType, AugType>(_left, _right, _split, _aug) {}

    inline bool ForceParallel() const { return this->aug; }
};

}  // namespace cpdd
