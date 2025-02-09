#ifndef PSTP_COVER_TREE_IMPL_COVER_INTER_NODE_HPP_
#define PSTP_COVER_TREE_IMPL_COVER_INTER_NODE_HPP_

#include "../cover_tree.h"
#include "pstp/dependence/tree_node.h"

namespace pstp {
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
struct CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::CoverInteriorNode
    : DynamicNode<Point, Splitter, AugType> {
  using BaseNode = DynamicNode<Point, Splitter, AugType>;
  using CoverNodeArr = typename BaseNode::NodeArr;
  using PT = Point;
  using ST = Splitter;
  using AT = AugType;

  CoverInteriorNode(CoverNodeArr const& _tree_nodes, const ST& _split,
                    const AT& _aug)
      : BaseNode(_tree_nodes, _split, _aug) {}

  inline void SetParallelFlag(bool const flag) { this->aug.emplace(flag); }

  inline void ResetParallelFlag() { this->aug.reset(); }

  inline bool GetParallelFlagIniStatus() { return this->aug.has_value(); }

  // NOTE: use a tri-state bool to indicate whether a subtree needs to be
  // rebuilt. If aug is not INITIALIZED, then it means there is no need to
  // rebuild; otherwise, the value depends on the initial tree size before
  // rebuilding.
  bool ForceParallel() const {
    return this->aug.has_value() ? this->aug.value()
                                 : this->size > BT::kSerialBuildCutoff;
  }
};

}  // namespace pstp

#endif  // PSTP_COVER_TREE_IMPL_COVER_INTER_NODE_HPP_
