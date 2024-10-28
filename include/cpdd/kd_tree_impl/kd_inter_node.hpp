#pragma once

#include "../kd_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template <typename Point, typename SplitRule, uint_fast8_t kBDO>
struct KdTree<Point, SplitRule, kBDO>::KdInteriorNode
    : BinaryNode<Point, Splitter, AugType> {
  using PT = Point;
  using ST = Splitter;
  using AT = AugType;

  KdInteriorNode(Node* _left, Node* _right, const ST& _split, const AT& _aug)
      : BinaryNode<Point, Splitter, AugType>(_left, _right, _split, _aug) {}

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

}  // namespace cpdd
