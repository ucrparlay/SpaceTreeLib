#ifndef PSPT_R_TREE_IMPL_R_INTER_NODE_HPP_
#define PSPT_R_TREE_IMPL_R_INTER_NODE_HPP_

#include "../r_tree.h"
#include "pspt/dependence/tree_node.h"

namespace pspt {
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
struct RTree<Point, SplitRule, kSkHeight, kImbaRatio>::RInteriorNode
    : BinaryNode<Point, Splitter, AugType> {
  using PT = Point;
  using ST = Splitter;
  using AT = AugType;

  RInteriorNode(Node* _left, Node* _right, const ST& _split, const AT& _aug)
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

}  // namespace pspt

#endif  // PSPT_R_TREE_IMPL_R_INTER_NODE_HPP_
