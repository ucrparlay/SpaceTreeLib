#ifndef PSI_R_TREE_IMPL_R_INTER_NODE_HPP_
#define PSI_R_TREE_IMPL_R_INTER_NODE_HPP_

#include "../r_tree.h"

#define RTREE_TEMPLATE template <typename Point, typename SplitRule, \
    uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
#define RTREE_CLASS RTree<Point, SplitRule, kSkHeight, kImbaRatio>

#include "psi/dependence/tree_node.h"

namespace psi {
RTREE_TEMPLATE

struct RTREE_CLASS::RInteriorNode
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


#undef RTREE_TEMPLATE
#undef RTREE_CLASS
}  // namespace psi

#endif  // PSI_R_TREE_IMPL_R_INTER_NODE_HPP_
