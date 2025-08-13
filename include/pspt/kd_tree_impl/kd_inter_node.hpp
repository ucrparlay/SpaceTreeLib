#ifndef PSPT_KD_TREE_IMPL_KD_INTER_NODE_HPP_
#define PSPT_KD_TREE_IMPL_KD_INTER_NODE_HPP_

#include <cstddef>

#include "../kd_tree.h"
#include "pspt/dependence/tree_node.h"

namespace pspt {

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
struct KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::LeafAugType {
  LeafAugType() : box(BT::GetEmptyBox()) {};
  LeafAugType(Box const& _box) : box(_box) {};
  LeafAugType(Slice In) : box(BT::GetBox(In)) {};

  Box& GetBox() { return this->box; }
  Box const& GetBox() const { return this->box; }

  void UpdateAug(Slice In) {
    this->box = BT::GetBox(box, BT::GetBox(In));
    return;
  }

  void Reset() {
    this->box = BT::GetEmptyBox();
    return;
  }

  Box box;
};

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
struct KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::InteriorAugType {
  InteriorAugType() : box(BT::GetEmptyBox()) { force_par_indicator.reset(); }
  InteriorAugType(Box const& _box) : box(_box) { force_par_indicator.reset(); }
  InteriorAugType(Node* l, Node* r)
      : box(BT::GetBox(BT::template RetriveBox<Leaf, Interior>(l),
                       BT::template RetriveBox<Leaf, Interior>(r))) {
    force_par_indicator.reset();
  }

  void SetParallelFlag(bool const flag) {
    this->force_par_indicator.emplace(flag);
  }

  void ResetParallelFlag() { this->force_par_indicator.reset(); }

  bool GetParallelFlagIniStatus() {
    return this->force_par_indicator.has_value();
  }

  bool ForceParallel(size_t sz) const {
    return this->force_par_indicator.has_value()
               ? this->force_par_indicator.value()
               : sz > BT::kSerialBuildCutoff;
  }

  Box& GetBox() { return this->box; }
  Box const& GetBox() const { return this->box; }

  void Update(Node* l, Node* r) {
    this->box = BT::GetBox(BT::template RetriveBox<Leaf, Interior>(l),
                           BT::template RetriveBox<Leaf, Interior>(r));
    return;
  }

  void Reset() {
    this->box = BT::GetEmptyBox();
    this->force_par_indicator.reset();
  }

  // NOTE: use a tri-state bool to indicate whether a subtree needs to be
  // rebuilt. If aug is not INITIALIZED, then it means there is no need to
  // rebuild; otherwise, the value depends on the initial tree size before
  // rebuilding.
  std::optional<bool> force_par_indicator;
  Box box;
};

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
struct KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::KdInteriorNode
    : BinaryNode<Point, Splitter, InteriorAugType> {
  using PT = Point;
  using ST = Splitter;
  using AT = InteriorAugType;

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

  auto UpdateAug(Node* l, Node* r) { return this->aug.Update(l, r); }

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

}  // namespace pspt

#endif  // PSPT_KD_TREE_IMPL_KD_INTER_NODE_HPP_
