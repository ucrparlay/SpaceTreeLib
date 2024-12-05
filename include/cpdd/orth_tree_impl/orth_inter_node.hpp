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
  bool ForceParallel() const {
    return this->aug.has_value() ? this->aug.value()
                                 : this->size > BT::kSerialBuildCutoff;
  }

  // NOTE: specific part

  // NOTE: Generate the hyperplane sequences for the node
  template <typename HyperPlaneSeq>
  void GenerateHyperPlaneSeq(HyperPlaneSeq& hyper_seq, auto idx, auto deep) {
    if (idx >= BaseNode::GetRegions()) {
      return;
    }
    hyper_seq[idx] = this->split[deep];
    GenerateHyperPlaneSeq(hyper_seq, idx * 2, deep + 1);
    GenerateHyperPlaneSeq(hyper_seq, idx * 2 + 1, deep + 1);
    return;
  }

  // NOTE: given an idx in the tree skeleton, modify its value its position in
  // left or right of the splitter
  BucketType SeievePoint(Point const& p, BucketType idx) {
    if constexpr (kMD == 2) {
      idx = 2 * idx + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[this->split[0].second], this->split[0].first));
      idx = 2 * idx + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[this->split[1].second], this->split[1].first));
    } else if constexpr (kMD == 3) {
      idx = 2 * idx + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[this->split[0].second], this->split[0].first));
      idx = 2 * idx + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[this->split[1].second], this->split[1].first));
      idx = 2 * idx + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[this->split[2].second], this->split[2].first));
    } else {
      for (BucketType i = 0; i < kMD; ++i) {
        idx = 2 * idx + 1 -
              static_cast<BucketType>(
                  Num::Lt(p.pnt[this->split[i].second], this->split[i].first));
      }
    }
    return idx;
  }

  template <typename BoxSeqSlice>
  void ComputeSubregionsRec(BoxSeqSlice box_seq, BucketType deep) {
    if (deep == kMD) {
      assert(box_seq.size() == 1);
      return;
    }

    BucketType n = box_seq.size();
    assert(n - n / 2 == n / 2);

    std::for_each_n(box_seq.begin(), n / 2, [&](auto& bx) {
      bx.second.pnt[this->split[deep].second] = this->split[deep].first;
    });
    std::for_each_n(box_seq.begin() + n / 2, n - n / 2, [&](auto& bx) {
      bx.first.pnt[this->split[deep].second] = this->split[deep].first;
    });

    ComputeSubregionsRec(box_seq.cut(0, n / 2), deep + 1);
    ComputeSubregionsRec(box_seq.cut(n / 2, n), deep + 1);
    return;
  }

  // NOTE: Given a box, produce new sub-boxes equivalent to modify the existing
  // boxes
  template <typename BoxSeq, typename Box>
  BoxSeq ComputeSubregions(Box const& box) {
    auto box_seq = BoxSeq(BaseNode::GetRegions(), box);
    ComputeSubregionsRec(parlay::make_slice(box_seq), 0);
    return std::move(box_seq);
  }

  // NOTE: Given a box and a bucket id, construct new box for that bucket
  template <typename Box>
  Box GetBoxByRegionId(BucketType id, Box const& box) {
    Box bx(box);
    assert(id >= 0 && id < BaseNode::GetRegions());

    if constexpr (kMD == 2) {  // dim = 2
      auto& target1 = (id & (1 << 1)) ? bx.first : bx.second;
      target1.pnt[this->split[0].second] = this->split[0].first;

      auto& target2 = (id & (1 << 0)) ? bx.first : bx.second;
      target2.pnt[this->split[1].second] = this->split[1].first;
    } else if constexpr (kMD == 3) {  // dim = 3
      auto& target1 = (id & (1 << 2)) ? bx.first : bx.second;
      target1.pnt[this->split[0].second] = this->split[0].first;

      auto& target2 = (id & (1 << 1)) ? bx.first : bx.second;
      target2.pnt[this->split[1].second] = this->split[1].first;

      auto& target3 = (id & (1 << 0)) ? bx.first : bx.second;
      target3.pnt[this->split[2].second] = this->split[2].first;
    } else {  // higher dim
      for (BucketType i = kMD; i > 0; --i) {
        auto& target = (id & (1 << (i - 1))) ? bx.first : bx.second;
        target.pnt[this->split[kMD - i].second] = this->split[kMD - i].first;
      }
    }

    return std::move(bx);
  }

  // NOTE: Given a box and a bucket id, modify the box for that bucket
  template <typename Box>
  inline void ModifyBoxById(BucketType id, Box& box) {
    assert(id >= 0 && id <= BaseNode::GetRegions());

    if constexpr (kMD == 2) {  // dim =2
      auto& target1 = (id & (1 << 1)) ? box.first : box.second;
      target1.pnt[this->split[0].second] = this->split[0].first;

      auto& target2 = (id & (1 << 0)) ? box.first : box.second;
      target2.pnt[this->split[1].second] = this->split[1].first;
    } else if constexpr (kMD == 3) {  // dim =3
      auto& target1 = (id & (1 << 2)) ? box.first : box.second;
      target1.pnt[this->split[0].second] = this->split[0].first;

      auto& target2 = (id & (1 << 1)) ? box.first : box.second;
      target2.pnt[this->split[1].second] = this->split[1].first;

      auto& target3 = (id & (1 << 0)) ? box.first : box.second;
      target3.pnt[this->split[2].second] = this->split[2].first;
    } else {  // dim > 3
      for (BucketType i = kMD; i > 0; --i) {
        auto& target = (id & (1 << (i - 1))) ? box.first : box.second;
        target.pnt[this->split[kMD - i].second] = this->split[kMD - i].first;
      }
    }
    return;
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
