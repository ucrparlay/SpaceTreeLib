#pragma once

#include <array>
#include "../quad_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint8_t kMD, uint8_t kBDO>
struct QuadTree<Point, SplitRule, kMD, kBDO>::QuadInteriorNode :
    MultiWayInteriorNode<Point, kMD, Splitter, AugType> {
    using Nodes = std::array<Node*, kMD>;
    using ST = Splitter;
    using AT = AugType;

    QuadInteriorNode(const Nodes& _tree_nodes, const ST& _split,
                     const AT& _aug) :
        MultiWayInteriorNode<Point, kMD, Splitter, AugType>(_tree_nodes, _split,
                                                            _aug) {}

    inline bool ForceParallel() const { return this->aug; }
};

}  // namespace cpdd