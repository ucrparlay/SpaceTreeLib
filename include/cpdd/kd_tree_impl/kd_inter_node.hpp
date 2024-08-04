#pragma once

#include "../kd_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint8_t kBDO>
struct KdTree<Point, SplitRule, kBDO>::KdInteriorNode :
    BinaryNode<Point, Splitter, AugType> {
    using PT = Point;
    using ST = Splitter;
    using AT = AugType;

    KdInteriorNode(Node* _left, Node* _right, const ST& _split,
                   const AT& _aug) :
        BinaryNode<Point, Splitter, AugType>(_left, _right, _split, _aug) {}

    inline bool ForceParallel() const { return this->aug; }
};

}  // namespace cpdd
