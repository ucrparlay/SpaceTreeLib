#ifndef PSI_COVER_TREE_IMPL_COVER_INTER_NODE_HPP_
#define PSI_COVER_TREE_IMPL_COVER_INTER_NODE_HPP_

#define COVERTREE_TEMPLATE                                              \
  template <typename Point, typename SplitRule, uint_fast8_t kSkHeight, \
            uint_fast8_t kImbaRatio>
#define COVERTREE_CLASS CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>

#include "../cover_tree.h"
#include "psi/dependence/tree_node.h"

namespace psi {
COVERTREE_TEMPLATE
struct COVERTREE_CLASS::CoverInteriorNode
    : DynamicNode<Point, Splitter, AugType> {
  using BaseNode = DynamicNode<Point, Splitter, AugType>;
  using CoverNodeArr = typename BaseNode::NodeArr;
  using PT = Point;
  using ST = Splitter;
  using AT = AugType;
  using CircleType = decltype(AT::cover_circle);

  CoverInteriorNode(CoverNodeArr const& _tree_nodes, const ST& _split,
                    const AT& _aug)
      : BaseNode(_tree_nodes, _split, _aug) {}

  auto const& GetCoverCircle() const { return this->aug.cover_circle; }

  auto const& GetParallelFlag() const { return this->aug.parallel_flag; }

  auto GetSubTreeNum() const { return this->tree_nodes.size(); }

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

}  // namespace psi

#undef COVERTREE_TEMPLATE
#undef COVERTREE_CLASS

#endif
