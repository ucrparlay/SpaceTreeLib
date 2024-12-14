#ifndef PSTP_ORTH_TREE_IMPL_ORTH_INTER_NODE_HPP_
#define PSTP_ORTH_TREE_IMPL_ORTH_INTER_NODE_HPP_

#include "../orth_tree.h"
#include "pstp/dependence/tree_node.h"

namespace pstp {
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
  // NOTE: compute the spliiter
  template <typename Box>
  static ST ComputeSplitter(Box const& box) {
    ST split;

    if constexpr (kMD == 2) {  // for dim = 2
      split[0] = HyperPlane(BT::GetBoxMid(0, box), 0);
      split[1] = HyperPlane(BT::GetBoxMid(1, box), 1);
    } else if constexpr (kMD == 3) {  // for dim = 3
      split[0] = HyperPlane(BT::GetBoxMid(0, box), 0);
      split[1] = HyperPlane(BT::GetBoxMid(1, box), 1);
      split[2] = HyperPlane(BT::GetBoxMid(2, box), 2);
    } else {
      for (DimsType i = 0; i < kMD; ++i) {  // for dim > 3
        split[i] = HyperPlane(BT::GetBoxMid(i, box), i);
      }
    }

    return std::move(split);
  }

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
  static BucketType SeievePoint(Point const& p, const ST& split,
                                BucketType idx) {
    if constexpr (kMD == 2) {  // for dim = 2
      idx = 2 * idx + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[split[0].second], split[0].first));
      idx = 2 * idx + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[split[1].second], split[1].first));
    } else if constexpr (kMD == 3) {  // for dim = 3
      idx = 2 * idx + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[split[0].second], split[0].first));
      idx = 2 * idx + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[split[1].second], split[1].first));
      idx = 2 * idx + 1 -
            static_cast<BucketType>(
                Num::Lt(p.pnt[split[2].second], split[2].first));
    } else {
      for (BucketType i = 0; i < kMD; ++i) {  // for dim > 3
        idx = 2 * idx + 1 -
              static_cast<BucketType>(
                  Num::Lt(p.pnt[split[i].second], split[i].first));
      }
    }
    return idx;
  }

  BucketType SeievePoint(Point const& p, BucketType idx) {
    return SeievePoint(p, this->split, idx);
  }

  // NOTE: Given a box, produce new sub-boxes equivalent to modify the existing
  // boxes
  template <typename BoxSeqSlice>
  static void ComputeSubregionsRec(BoxSeqSlice box_seq, const ST& split,
                                   BucketType deep) {
    if (deep == kMD) {
      assert(box_seq.size() == 1);
      return;
    }

    BucketType n = box_seq.size();
    assert(n - n / 2 == n / 2);

    std::for_each_n(box_seq.begin(), n / 2, [&](auto& bx) {
      bx.second.pnt[split[deep].second] = split[deep].first;
    });
    std::for_each_n(box_seq.begin() + n / 2, n - n / 2, [&](auto& bx) {
      bx.first.pnt[split[deep].second] = split[deep].first;
    });

    ComputeSubregionsRec(box_seq.cut(0, n / 2), split, deep + 1);
    ComputeSubregionsRec(box_seq.cut(n / 2, n), split, deep + 1);
    return;
  }

  template <typename BoxSeq, typename Box>
  BoxSeq ComputeSubregions(Box const& box) {
    auto box_seq = BoxSeq(BaseNode::GetRegions(), box);
    ComputeSubregionsRec(parlay::make_slice(box_seq), this->split, 0);
    return std::move(box_seq);
  }

  template <typename BoxSeq, typename Box>
  static BoxSeq ComputeSubregions(Box const& box, ST const& split) {
    auto box_seq = BoxSeq(BaseNode::GetRegions(), box);
    ComputeSubregionsRec(parlay::make_slice(box_seq), split, 0);
    return std::move(box_seq);
  }

  // NOTE: Given a box and a bucket id, modify the box for that bucket
  template <typename Box>
  static void ModifyBoxById(BucketType id, const ST& split, Box& box) {
    assert(id >= 0 && id <= BaseNode::GetRegions());

    if constexpr (kMD == 2) {  // dim =2
      auto& target1 = (id & (1 << 1)) ? box.first : box.second;
      target1.pnt[split[0].second] = split[0].first;

      auto& target2 = (id & (1 << 0)) ? box.first : box.second;
      target2.pnt[split[1].second] = split[1].first;
    } else if constexpr (kMD == 3) {  // dim =3
      auto& target1 = (id & (1 << 2)) ? box.first : box.second;
      target1.pnt[split[0].second] = split[0].first;

      auto& target2 = (id & (1 << 1)) ? box.first : box.second;
      target2.pnt[split[1].second] = split[1].first;

      auto& target3 = (id & (1 << 0)) ? box.first : box.second;
      target3.pnt[split[2].second] = split[2].first;
    } else {  // dim > 3
      for (BucketType i = kMD; i > 0; --i) {
        auto& target = (id & (1 << (i - 1))) ? box.first : box.second;
        target.pnt[split[kMD - i].second] = split[kMD - i].first;
      }
    }
    return;
  }

  template <typename Box>
  void ModifyBoxById(BucketType id, Box& box) {
    ModifyBoxById(id, this->split, box);
    return;
  }

  // NOTE: return a new box by modifying the box for a specific region
  template <typename Box>
  static Box GetBoxByRegionId(BucketType id, ST const& split, Box const& box) {
    Box new_box(box);
    ModifyBoxById(id, split, new_box);
    return std::move(new_box);
  }

  template <typename Box>
  Box GetBoxByRegionId(BucketType id, Box const& box) {
    Box new_box(box);
    ModifyBoxById(id, this->split, new_box);
    return std::move(new_box);
  }
};

}  // namespace pstp

#endif  // PSTP_ORTH_TREE_IMPL_ORTH_INTER_NODE_HPP_
