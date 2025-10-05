#ifndef PSI_BASE_TREE_IMPL_INNER_TREE_MULTI_HPP_
#define PSI_BASE_TREE_IMPL_INNER_TREE_MULTI_HPP_

#include "dependence/concepts.h"
#include "inner_tree_common.hpp"

namespace psi {
namespace inner_tree_detail {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio, typename Leaf, typename Interior>
  requires IsMultiNode<Interior>
struct MultiNodeOps {
  using BT = BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>;
  using InnerTree = typename BT::template InnerTree<Leaf, Interior>;
  using Box = typename BT::Box;
  using BoxSeq = typename BT::BoxSeq;
  using NodeTagSeq = typename BT::NodeTagSeq;
  using NodeTag = typename BT::NodeTag;
  using BucketType = typename BT::BucketType;
  using BallsType = typename BT::BallsType;
  using DimsType = typename BT::DimsType;

  static Box GetBoxByRegionIdx(InnerTree* self, int const idx, Box const& box) {
    assert(InnerTree::GetDepthByIndex(19) == 4);
    assert(BT::kDim ==
           static_cast<BucketType>(std::log2(Interior::GetRegions())));

    Box bx(box);
    BucketType h = InnerTree::GetDepthByIndex(idx);
    for (BucketType i = h, new_idx = 1; i > 0; i -= BT::kDim) {
      BucketType local_id = (idx >> (i - BT::kDim)) & ((1 << BT::kDim) - 1);
      auto TI = static_cast<Interior*>(self->tags[new_idx].first);
      TI->ModifyBoxById(local_id, bx);
      new_idx = new_idx << BT::kDim | local_id;
    }
    return bx;
  }

  static void AssignNodeTag(InnerTree* self, Node* T, BucketType idx) {
    if (T->is_leaf || idx > BT::kPivotNum) {
      assert(self->tags_num < BT::kBucketNum);
      self->tags[idx] = NodeTag(T, self->tags_num);
      self->rev_tag[self->tags_num++] = idx;
      return;
    }

    self->tags[idx] = NodeTag(T, BT::kBucketNum);
    Interior* TI = static_cast<Interior*>(T);
    for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
      AssignNodeTag(self, TI->tree_nodes[i], idx * Interior::GetRegions() + i);
    }
  }

  static void ReduceSums(InnerTree* self, BucketType idx) {
    if (idx > BT::kPivotNum || self->tags[idx].first->is_leaf) {
      assert(self->tags[idx].second < BT::kBucketNum);
      self->sums_tree[idx] = self->sums[self->tags[idx].second];
      return;
    }

    self->sums_tree[idx] = 0;
    for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
      ReduceSums(self, idx * Interior::GetRegions() + i);
      self->sums_tree[idx] += self->sums_tree[idx * Interior::GetRegions() + i];
    }
  }

  template <typename ViolateFunc>
  static void PickTag(InnerTree* self, BucketType idx,
                      ViolateFunc&& violate_func) {
    if (idx > BT::kPivotNum || self->tags[idx].first->is_leaf) {
      self->tags[idx].second = BT::kBucketNum + 1;
      self->rev_tag[self->tags_num++] = idx;
      return;
    }

    assert(self->tags[idx].second == BT::kBucketNum &&
           (!self->tags[idx].first->is_leaf));
    if (InvokeWithOptionalArg<bool>(violate_func, idx)) {
      self->tags[idx].second = BT::kBucketNum + 2;
      self->rev_tag[self->tags_num++] = idx;
      return;
    }

    for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
      PickTag(self, idx * Interior::GetRegions() + i, violate_func);
    }
  }

  template <bool kSetParallelFlag, typename ViolateFunc>
  static void MarkNodesForRebuild(InnerTree* self, BucketType idx,
                                  BucketType& re_num, size_t& tot_re_size,
                                  BoxSeq& box_seq, Box const& box,
                                  bool has_tomb, ViolateFunc&& violate_func) {
    if (idx > BT::kPivotNum || self->tags[idx].first->is_leaf) {
      assert(self->tags[idx].second >= 0 &&
             self->tags[idx].second < BT::kBucketNum);
      if (!has_tomb) {
        self->tags[idx].second = BT::kBucketNum + 2;
        if constexpr (kSetParallelFlag) {
          if (!self->tags[idx].first->is_leaf) {
            auto TI = static_cast<Interior*>(self->tags[idx].first);
            TI->SetParallelFlag(TI->size > BT::kSerialBuildCutoff);
          }
        }
      } else {
        self->tags[idx].second = BT::kBucketNum + 1;
      }

      box_seq[self->tags_num] = box;
      self->rev_tag[self->tags_num++] = idx;
      return;
    }

    assert(self->tags[idx].second == BT::kBucketNum &&
           (!self->tags[idx].first->is_leaf));
    Interior* TI = static_cast<Interior*>(self->tags[idx].first);
    if (has_tomb && InvokeWithOptionalArg<bool>(violate_func, idx)) {
      self->tags[idx].second = BT::kBucketNum + 3;
      has_tomb = false;
      TI->SetParallelFlag(TI->size > BT::kSerialBuildCutoff);
      re_num++, tot_re_size += TI->size;
    }

    for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
      MarkNodesForRebuild<kSetParallelFlag>(
          self, idx * Interior::GetRegions() + i, re_num, tot_re_size, box_seq,
          TI->GetBoxByRegionId(i, box), has_tomb, violate_func);
    }
  }

  static void TagOversizedNodesRecursive(InnerTree* self, BucketType idx) {
    if (idx > BT::kPivotNum || self->tags[idx].first->is_leaf) {
      self->tags[idx].second = BT::kBucketNum + 2;
      if (!self->tags[idx].first->is_leaf) {
        Interior* TI = static_cast<Interior*>(self->tags[idx].first);
        TI->SetParallelFlag(TI->size > BT::kSerialBuildCutoff);
      }
      return;
    }

    assert(!self->tags[idx].first->is_leaf);
    Interior* TI = static_cast<Interior*>(self->tags[idx].first);
    TI->SetParallelFlag(TI->size > BT::kSerialBuildCutoff);

    for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
      TagOversizedNodesRecursive(self, idx * Interior::GetRegions() + i);
    }
  }


  template <bool UpdateParFlag, typename NodeOrNodeBox>
  static NodeOrNodeBox UpdateInnerTreePointersRecursive(
      InnerTree* self, BucketType idx,
      parlay::sequence<NodeOrNodeBox> const& tree_nodes, BucketType& p) {
    if (self->tags[idx].second == BT::kBucketNum + 1 ||
        self->tags[idx].second == BT::kBucketNum + 2) {
      return tree_nodes[p++];
    }

    typename Interior::OrthNodeArr new_nodes;
    for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
      new_nodes[i] = UpdateInnerTreePointersRecursive<UpdateParFlag>(
          self, idx * Interior::GetRegions() + i, tree_nodes, p);
    }

    BT::template UpdateInterior<Interior, UpdateParFlag>(self->tags[idx].first,
                                                         new_nodes);
    return self->tags[idx].first;
  }

  template <bool UpdateParFlag, typename NodeOrNodeBox, typename Func>
  static NodeOrNodeBox TagNodesForRebuildRecursive(
      InnerTree* self, BucketType idx,
      parlay::sequence<NodeOrNodeBox> const& tree_nodes, BucketType& p,
      Func&& func) {
    if (self->tags[idx].second == BT::kBucketNum + 1 ||
        self->tags[idx].second == BT::kBucketNum + 2) {
      return tree_nodes[p++];
    }

    typename Interior::OrthNodeArr new_nodes;
    for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
      new_nodes[i] = TagNodesForRebuildRecursive<UpdateParFlag>(
          self, idx * Interior::GetRegions() + i, tree_nodes, p, func);
    }

    BT::template UpdateInterior<Interior, UpdateParFlag>(self->tags[idx].first,
                                                         new_nodes);
    if (self->tags[idx].second == BT::kBucketNum + 3) {
      func(idx);
    }
    return self->tags[idx].first;
  }

  template <bool UpdateParFlag, typename NodeOrNodeBox, typename Func>
  static NodeOrNodeBox UpdateAfterDeletionRecursive(
      InnerTree* self, BucketType idx,
      parlay::sequence<NodeOrNodeBox> const& tree_nodes, BucketType& p,
      Func&& func) {
    if (self->tags[idx].second == BT::kBucketNum + 1 ||
        self->tags[idx].second == BT::kBucketNum + 2) {
      return tree_nodes[p++];
    }

    if (self->tags[idx].second == BT::kBucketNum + 3) {
      func(0);
      assert(func(1) == true);
    }

    typename Interior::OrthNodeArr new_nodes;
    for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
      new_nodes[i] = UpdateAfterDeletionRecursive<UpdateParFlag>(
          self, idx * Interior::GetRegions() + i, tree_nodes, p, func);
    }

    if (!func(1)) {
      BT::template UpdateInterior<Interior, UpdateParFlag>(
          self->tags[idx].first, new_nodes);
      return self->tags[idx].first;
    } else if (self->tags[idx].second == BT::kBucketNum + 3) {
      func(0);
      assert(func(1) == false);
      if (!self->tags[idx].first->is_leaf) {
        static_cast<Interior*>(self->tags[idx].first)->ResetParallelFlag();
      }
      return self->tags[idx].first;
    } else {
      return nullptr;
    }
  }
};

}  // namespace inner_tree_detail
}  // namespace psi

#endif  // PSI_BASE_TREE_IMPL_INNER_TREE_MULTI_HPP_
