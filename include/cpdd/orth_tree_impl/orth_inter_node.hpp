#pragma once

#include "../orth_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
struct OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::OrthInteriorNode
    : MultiNode<Point, kMD, Splitter, AugType> {
  using BaseNode = MultiNode<Point, kMD, Splitter, AugType>;
  using OrthNodeArr = typename BaseNode::NodeArr;
  using PT = Point;
  using ST = Splitter;
  using AT = AugType;

  OrthInteriorNode(OrthNodeArr const& _tree_nodes, const ST& _split,
                   const AT& _aug)
      : BaseNode(_tree_nodes, _split, _aug) {}

  inline void SetParallelFlag(bool const flag) { this->aug.emplace(flag); }

  inline void ResetParallelFlag() { this->aug.reset(); }

  inline bool GetParallelFlagIniStatus() { return this->aug.has_value(); }

  // NOTE: use a tri-state bool to indicate whether a subtree needs to be
  // rebuilt. If aug is not INITIALIZED, then it means there is no need to
  // rebuild; otherwise, the value depends on the initial tree size before
  // rebuilding.
  inline bool ForceParallel() const {
    return this->aug.has_value() ? this->aug.value()
                                 : this->size > BT::kSerialBuildCutoff;
  }
};

// // NOTE: To expand as a kdtree node
// template<typename Point, typename SplitRule, uint_fast8_t kMD,
//          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
// struct OrthTree<Point, SplitRule, kMD, kSkHeight,kImbaRatio>::KdInteriorNode
// :
//     BinaryNode<Point, HyperPlane, AugType> {
//     using PT = Point;
//     using ST = HyperPlane;
//     using AT = AugType;
//
//     KdInteriorNode(Node* _left, Node* _right, const ST& _split,
//                    const AT& _aug) :
//         BinaryNode<Point, HyperPlane, AugType>(_left, _right, _split, _aug)
//         {}
//
//     inline void SetParallelFlag(bool flag) { this->aug = AT(flag); }
//
//     inline void ResetParallelFlag() { this->aug = false; }
//
//     inline bool ForceParallel() const {
//         return this->aug ? *(this->aug) : this->size >
//         BT::kSerialBuildCutoff;
//     }
// };
}  // namespace cpdd
