#ifndef PSPT_KD_TREE_IMPL_KD_INTER_NODE_HPP_
#define PSPT_KD_TREE_IMPL_KD_INTER_NODE_HPP_

#include "../kd_tree.h"
#include "pspt/dependence/tree_node.h"

namespace pspt {
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
struct KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::KdAugType {
  std::optional<bool> parallel_flag;
  Box box;

  inline void SetParallelFlag(bool const flag) {
    this->parallel_flag.emplace(flag);
  }

  inline void ResetParallelFlag() { this->parallel_flag.reset(); }

  inline bool GetParallelFlagIniStatus() {
    return this->parallel_flag.has_value();
  }

  // NOTE: use a tri-state bool to indicate whether a subtree needs to be
  // rebuilt. If aug is not INITIALIZED, then it means there is no need to
  // rebuild; otherwise, the value depends on the initial tree size before
  // rebuilding.
  inline bool ForceParallel() const {
    // BUG: this is not consistent as the original one
    return this->parallel_flag.has_value() ? this->parallel_flag.value() : true;
  }
};

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
struct KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::KdInteriorNode
    : BinaryNode<Point, Splitter, AugType> {
  using PT = Point;
  using ST = Splitter;
  using AT = AugType;

  KdInteriorNode(Node* _left, Node* _right, const ST& _split, const AT& _aug)
      : BinaryNode<Point, Splitter, AugType>(_left, _right, _split, _aug) {}

  // Adding a virtual destructor makes Node polymorphic
  virtual ~KdInteriorNode() = default;
  inline void SetParallelFlag(bool const flag) {
    this->aug.SetParallelFlag(flag);
  }

  inline void ResetParallelFlag() { this->aug.ResetParallelFlag(); }

  inline bool GetParallelFlagIniStatus() { return false; }

  // NOTE: use a tri-state bool to indicate whether a subtree needs to be
  // rebuilt. If aug is not INITIALIZED, then it means there is no need to
  // rebuild; otherwise, the value depends on the initial tree size before
  // rebuilding.
  inline bool ForceParallel() const { return true; }
};

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <uint_fast8_t kMD>
struct KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::KdCompressionNode
    : MultiNode<Point, kMD, CompressNodeSplitter, AugType> {
  using BaseNode = MultiNode<Point, kMD, CompressNodeSplitter, AugType>;
  using KdNodeArr = typename BaseNode::NodeArr;
  using PT = Point;
  using ST = CompressNodeSplitter;
  using AT = AugType;

  KdCompressionNode(KdNodeArr const& _tree_nodes, const ST& _split,
                    const AT& _aug)
      : BaseNode(_tree_nodes, _split, _aug) {}

  // Adding a virtual destructor makes Node polymorphic
  virtual ~KdCompressionNode() = default;

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

#endif  // PSPT_KD_TREE_IMPL_KD_INTER_NODE_HPP_
