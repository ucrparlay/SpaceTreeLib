#pragma once

#include <array>
#include "../orth_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint8_t kMD, uint8_t kBDO>
struct OrthTree<Point, SplitRule, kMD, kBDO>::QuadInteriorNode :
    MultiNode<Point, kMD, Splitter, AugType> {
    using Nodes = std::array<Node*, 1 << kMD>;
    using PT = Point;
    using ST = Splitter;
    using AT = AugType;

    static_assert(
        Nodes().size() ==
        typename MultiNode<Point, kMD, Splitter, AugType>::Nodes().size());

    QuadInteriorNode(const Nodes& _tree_nodes, const ST& _split,
                     const AT& _aug) :
        MultiNode<Point, kMD, Splitter, AugType>(_tree_nodes, _split, _aug) {}

    inline bool ForceParallel() const { return this->aug; }
};

}  // namespace cpdd
