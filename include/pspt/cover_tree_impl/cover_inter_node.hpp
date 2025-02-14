#ifndef PSPT_COVER_TREE_IMPL_COVER_INTER_NODE_HPP_
#define PSPT_COVER_TREE_IMPL_COVER_INTER_NODE_HPP_

#include "../cover_tree.h"
#include "pspt/dependence/tree_node.h"

namespace pspt {
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

  auto& GetCoverCircle() const { return this->aug.cover_circle; }
  auto& GetParallelFlag() const { return this->aug.parallel_flag; }

  inline void SetParallelFlag(bool const flag) {
    this->GetParallelFlag().emplace(flag);
  }

  inline void ResetParallelFlag() { this->GetParallelFlag().reset(); }

  inline bool GetParallelFlagIniStatus() {
    return this->GetParallelFlag().has_value();
  }

  // NOTE: use a tri-state bool to indicate whether a subtree needs to be
  // rebuilt. If aug is not INITIALIZED, then it means there is no need to
  // rebuild; otherwise, the value depends on the initial tree size before
  // rebuilding.
  bool ForceParallel() const {
    return this->GetParallelFlag().has_value()
               ? this->GetParallelFlag().value()
               : this->size > BT::kSerialBuildCutoff;
  }
};

}  // namespace pspt

#endif  // PSPT_COVER_TREE_IMPL_COVER_INTER_NODE_HPP_
