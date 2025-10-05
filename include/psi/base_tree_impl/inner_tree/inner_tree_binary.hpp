#ifndef PSI_BASE_TREE_IMPL_INNER_TREE_BINARY_HPP_
#define PSI_BASE_TREE_IMPL_INNER_TREE_BINARY_HPP_

#include "inner_tree_common.hpp"

namespace psi {
namespace inner_tree_detail {

template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio, typename Leaf, typename Interior>
  requires IsBinaryNode<Interior>
struct BinaryNodeOps {
  using BT = BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>;
  using InnerTree = typename BT::template InnerTree<Leaf, Interior>;
  using Box = typename BT::Box;
  using BoxSeq = typename BT::BoxSeq;
  using BoxCut = typename BT::BoxCut;
  using NodeTagSeq = typename BT::NodeTagSeq;
  using NodeTag = typename BT::NodeTag;
  using BucketType = typename BT::BucketType;
  using BallsType = typename BT::BallsType;

  static Box GetBoxByRegionIdx(InnerTree* self, int const idx, Box const& box) {
    assert(InnerTree::GetDepthByIndex(10) == 3);
    Box bx(box);
    int h = InnerTree::GetDepthByIndex(idx);
    for (int i = h - 1, new_idx = 1; i >= 0; i--) {
      int local_id = (idx >> i) & 1;
      auto TI = static_cast<Interior*>(self->tags[new_idx].first);
      auto& target = local_id ? bx.first : bx.second;
      target.pnt[TI->split.second] = TI->split.first;
      new_idx = new_idx << 1 | local_id;
      assert(new_idx <= idx);
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
    AssignNodeTag(self, TI->left, idx << 1);
    AssignNodeTag(self, TI->right, idx << 1 | 1);
  }

  static void ReduceSums(InnerTree* self, BucketType idx) {
    if (idx > BT::kPivotNum || self->tags[idx].first->is_leaf) {
      assert(self->tags[idx].second < BT::kBucketNum);
      self->sums_tree[idx] = self->sums[self->tags[idx].second];
      return;
    }

    ReduceSums(self, idx << 1);
    ReduceSums(self, idx << 1 | 1);
    self->sums_tree[idx] =
        self->sums_tree[idx << 1] + self->sums_tree[idx << 1 | 1];
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

    PickTag(self, idx << 1, violate_func);
    PickTag(self, idx << 1 | 1, violate_func);
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

    BoxCut box_cut(box, TI->split, true);
    MarkNodesForRebuild<kSetParallelFlag>(self, idx << 1, re_num, tot_re_size,
                                          box_seq, box_cut.GetFirstBoxCut(),
                                          has_tomb, violate_func);
    MarkNodesForRebuild<kSetParallelFlag>(
        self, idx << 1 | 1, re_num, tot_re_size, box_seq,
        box_cut.GetSecondBoxCut(), has_tomb, violate_func);
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

    TagOversizedNodesRecursive(self, idx << 1);
    TagOversizedNodesRecursive(self, idx << 1 | 1);
  }

  template <bool UpdateParFlag, typename NodeOrNodeBox>
  static NodeOrNodeBox UpdateInnerTreePointersRecursive(
      InnerTree* self, BucketType idx,
      parlay::sequence<NodeOrNodeBox> const& tree_nodes, BucketType& p) {
    using NodeBox = typename BT::NodeBox;

    if (self->tags[idx].second == BT::kBucketNum + 1 ||
        self->tags[idx].second == BT::kBucketNum + 2) {
      return tree_nodes[p++];
    }

    NodeOrNodeBox const& left = UpdateInnerTreePointersRecursive<UpdateParFlag>(
        self, idx << 1, tree_nodes, p);
    NodeOrNodeBox const& right = UpdateInnerTreePointersRecursive<UpdateParFlag>(
        self, idx << 1 | 1, tree_nodes, p);

    BT::template UpdateInterior<Interior, UpdateParFlag>(self->tags[idx].first,
                                                         left, right);
    if constexpr (IsPointerToNode<NodeOrNodeBox>) {
      return self->tags[idx].first;
    } else {
      return NodeBox(self->tags[idx].first, Box());
    }
  }

  template <bool UpdateParFlag, typename NodeOrNodeBox>
  static NodeOrNodeBox UpdateInnerTreePointersWithBoxRecursive(
      InnerTree* self, BucketType idx,
      parlay::sequence<NodeOrNodeBox> const& tree_nodes, BucketType& p) {
    using NodeBox = typename BT::NodeBox;

    if (self->tags[idx].second == BT::kBucketNum + 1 ||
        self->tags[idx].second == BT::kBucketNum + 2) {
      return tree_nodes[p++];
    }

    NodeOrNodeBox const& left =
        UpdateInnerTreePointersWithBoxRecursive<UpdateParFlag>(self, idx << 1,
                                                               tree_nodes, p);
    NodeOrNodeBox const& right =
        UpdateInnerTreePointersWithBoxRecursive<UpdateParFlag>(
            self, idx << 1 | 1, tree_nodes, p);

    BT::template UpdateInterior<Interior, UpdateParFlag>(self->tags[idx].first,
                                                         left, right);
    return NodeBox(self->tags[idx].first,
                   BT::GetBox(left.second, right.second));
  }

  template <bool UpdateParFlag, typename NodeOrNodeBox, typename Func>
  static NodeOrNodeBox TagNodesForRebuildRecursive(
      InnerTree* self, BucketType idx,
      parlay::sequence<NodeOrNodeBox> const& tree_nodes, BucketType& p,
      Func&& func) {
    using NodeBox = typename BT::NodeBox;

    if (self->tags[idx].second == BT::kBucketNum + 1 ||
        self->tags[idx].second == BT::kBucketNum + 2) {
      return tree_nodes[p++];
    }

    NodeOrNodeBox const& left = TagNodesForRebuildRecursive<UpdateParFlag>(
        self, idx << 1, tree_nodes, p, func);
    NodeOrNodeBox const& right = TagNodesForRebuildRecursive<UpdateParFlag>(
        self, idx << 1 | 1, tree_nodes, p, func);

    BT::template UpdateInterior<Interior, UpdateParFlag>(self->tags[idx].first,
                                                         left, right);
    auto new_box = BT::GetBox(left.second, right.second);
    if (self->tags[idx].second == BT::kBucketNum + 3) {
      func(new_box, idx);
    }
    return NodeBox(self->tags[idx].first, std::move(new_box));
  }

  template <bool UpdateParFlag, typename NodeOrNodeBox, typename Func>
  static NodeOrNodeBox UpdateAfterDeletionRecursive(
      InnerTree* self, BucketType idx,
      parlay::sequence<NodeOrNodeBox> const& tree_nodes, BucketType& p,
      Func&& func) {
    using NodeBox = typename BT::NodeBox;

    if (self->tags[idx].second == BT::kBucketNum + 1 ||
        self->tags[idx].second == BT::kBucketNum + 2) {
      return tree_nodes[p++];
    }

    if (self->tags[idx].second == BT::kBucketNum + 3) {
      func(0);
      assert(func(1) == true);
    }

    NodeOrNodeBox const& left = UpdateAfterDeletionRecursive<UpdateParFlag>(
        self, idx << 1, tree_nodes, p, func);
    NodeOrNodeBox const& right = UpdateAfterDeletionRecursive<UpdateParFlag>(
        self, idx << 1 | 1, tree_nodes, p, func);

    if (!func(1)) {
      BT::template UpdateInterior<Interior, UpdateParFlag>(
          self->tags[idx].first, left, right);
      return NodeBox(self->tags[idx].first, Box());
    } else if (self->tags[idx].second == BT::kBucketNum + 3) {
      func(0);
      if (!self->tags[idx].first->is_leaf) {
        static_cast<Interior*>(self->tags[idx].first)->ResetParallelFlag();
      }
      assert(func(1) == false);
      return NodeBox(self->tags[idx].first, Box());
    } else {
      return NodeBox(nullptr, Box());
    }
  }
};

}  // namespace inner_tree_detail
}  // namespace psi

#endif  // PSI_BASE_TREE_IMPL_INNER_TREE_BINARY_HPP_
