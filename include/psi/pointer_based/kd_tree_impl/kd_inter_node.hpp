#ifndef PSI_KD_TREE_IMPL_KD_INTER_NODE_HPP_
#define PSI_KD_TREE_IMPL_KD_INTER_NODE_HPP_

#include <cstddef>

#include "../kd_tree.h"
#include "../../dependence/tree_node.h"

#define KDTREE_TEMPLATE                                               \
  template <typename Point, typename SplitRule, typename LeafAugType, \
            typename InteriorAugType, uint_fast8_t kSkHeight,         \
            uint_fast8_t kImbaRatio>
#define KDTREE_CLASS \
  KdTree<Point, SplitRule, LeafAugType, InteriorAugType, kSkHeight, kImbaRatio>

namespace psi {

KDTREE_TEMPLATE
struct KDTREE_CLASS::KdInteriorNode
    : BinaryNode<Point, Splitter, InteriorAugType> {
  using PT = Point;
  using ST = Splitter;
  using AT = InteriorAugType;

  KdInteriorNode(Node* _left, Node* _right, const ST& _split)
      : BinaryNode<Point, Splitter, AT>(
            _left, _right, _split,
            AT(AT::template Create<Leaf, Interior>(_left, _right))) {}

  KdInteriorNode(Node* _left, Node* _right, const ST& _split, const AT& _aug)
      : BinaryNode<Point, Splitter, AT>(_left, _right, _split, _aug) {}

  // Adding a virtual destructor makes Node polymorphic
  virtual ~KdInteriorNode() = default;

  inline void SetParallelFlag(bool const flag) {
    this->aug.SetParallelFlag(flag);
  }

  inline void ResetParallelFlag() { this->aug.ResetParallelFlag(); }

  inline bool GetParallelFlagIniStatus() {
    return this->aug.GetParallelFlagIniStatus();
  }

  inline bool ForceParallel() const {
    return this->aug.ForceParallel(this->size);
  }

  auto UpdateAug(Node* l, Node* r) {
    return this->aug.template Update<Leaf, Interior>(l, r);
  }

  auto GetBox()
    requires HasBox<AT>
  {
    return this->aug.GetBox();
  }

  auto GetBox() const
    requires HasBox<AT>
  {
    return this->aug.GetBox();
  }

  auto ResetAug() { return this->aug.Reset(); }
};

// template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
//           uint_fast8_t kImbaRatio>
// template <uint_fast8_t kMD>
// struct KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::KdCompressionNode
//     : MultiNode<Point, kMD, CompressNodeSplitter, AugType> {
//   using BaseNode = MultiNode<Point, kMD, CompressNodeSplitter, AugType>;
//   using KdNodeArr = typename BaseNode::NodeArr;
//   using PT = Point;
//   using ST = CompressNodeSplitter;
//   using AT = AugType;

//   KdCompressionNode(KdNodeArr const& _tree_nodes, const ST& _split,
//                     const AT& _aug)
//       : BaseNode(_tree_nodes, _split, _aug) {}

//   // Adding a virtual destructor makes Node polymorphic
//   virtual ~KdCompressionNode() = default;

//   inline void SetParallelFlag(bool const flag) { this->aug.emplace(flag); }

//   inline void ResetParallelFlag() { this->aug.reset(); }

//   inline bool GetParallelFlagIniStatus() { return this->aug.has_value(); }

//   // NOTE: use a tri-state bool to indicate whether a subtree needs to be
//   // rebuilt. If aug is not INITIALIZED, then it means there is no need to
//   // rebuild; otherwise, the value depends on the initial tree size before
//   // rebuilding.
//   inline bool ForceParallel() const {
//     return this->aug.has_value() ? this->aug.value()
//                                  : this->size > BT::kSerialBuildCutoff;
//   }
// };

}  // namespace psi

#undef KDTREE_TEMPLATE
#undef KDTREE_CLASS

#endif  // PSI_KD_TREE_IMPL_KD_INTER_NODE_HPP_
