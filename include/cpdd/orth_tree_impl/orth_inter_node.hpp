#pragma once

#include <array>
#include "../orth_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint_fast8_t kMD,
         uint_fast8_t kBDO>
struct OrthTree<Point, SplitRule, kMD, kBDO>::OrthInteriorNode :
    MultiNode<Point, kMD, Splitter, AugType> {
    using MNode = MultiNode<Point, kMD, Splitter, AugType>;
    using OrthNodeSeq = std::array<Node*, MNode::kRegions>;
    using PT = Point;
    using ST = Splitter;
    using AT = AugType;

    static_assert(
        OrthNodeSeq().size() ==
        typename MultiNode<Point, kMD, Splitter, AugType>::OrthNodeSeq()
            .size());

    OrthInteriorNode(const OrthNodeSeq& _tree_nodes, const ST& _split,
                     const AT& _aug) :
        MNode(_tree_nodes, _split, _aug) {}

    inline void SetParallelFlag(bool flag) { this->aug = AT(flag); }

    inline void ResetParallelFlag() { this->aug = false; }

    inline bool ForceParallel() const {
        return this->aug ? *(this->aug) : this->size > BT::kSerialBuildCutoff;
    }
};

// NOTE: To expand as a kdtree node
template<typename Point, typename SplitRule, uint_fast8_t kMD,
         uint_fast8_t kBDO>
struct OrthTree<Point, SplitRule, kMD, kBDO>::KdInteriorNode :
    BinaryNode<Point, HyperPlane, AugType> {
    using PT = Point;
    using ST = HyperPlane;
    using AT = AugType;

    KdInteriorNode(Node* _left, Node* _right, const ST& _split,
                   const AT& _aug) :
        BinaryNode<Point, HyperPlane, AugType>(_left, _right, _split, _aug) {}

    inline void SetParallelFlag(bool flag) { this->aug = AT(flag); }

    inline void ResetParallelFlag() { this->aug = false; }

    inline bool ForceParallel() const {
        return this->aug ? *(this->aug) : this->size > BT::kSerialBuildCutoff;
    }
};
}  // namespace cpdd
