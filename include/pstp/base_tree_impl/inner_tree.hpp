#ifndef PSTP_BASE_TREE_IMPL_INNER_TREE_HPP_
#define PSTP_BASE_TREE_IMPL_INNER_TREE_HPP_

#include <utility>

#include "../base_tree.h"
#include "dependence/concepts.h"

namespace pstp {
template <typename Point, typename DerivedTree, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Leaf, typename Interior>
struct BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>::InnerTree {
  using BT = BaseTree<Point, DerivedTree, kSkHeight, kImbaRatio>;
  InnerTree(BT& _btref)
      : BTRef(_btref),
        tags_num(0),
        tags(NodeTagSeq::uninitialized(kPivotNum + kBucketNum + 1)),
        sums_tree(parlay::sequence<BallsType>(kPivotNum + kBucketNum + 1)),
        rev_tag(BucketSeq::uninitialized(kBucketNum)) {}

  // NOTE: helpers
  bool AssertSize(Node* T) const {
    if (T->is_leaf) {
      Leaf* TI = static_cast<Leaf*>(T);
      assert(T->size <= TI->pts.size() && T->size <= kLeaveWrap);
      return true;
    }
    Interior* TI = static_cast<Interior*>(T);
    if constexpr (IsBinaryNode<Interior>) {
      assert(TI->size == TI->left->size + TI->right->size);
    } else if constexpr (IsMultiNode<Interior>) {
      assert(std::cmp_equal(
          TI->size,
          std::accumulate(TI->tree_nodes.begin(), TI->tree_nodes.end(), 0,
                          [](size_t sum, Node* T) { return sum + T->size; })));
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
    return true;
  }

  inline BucketType GetNodeIdx(BucketType idx, Node* T) {
    if (tags[idx].first == T) return idx;
    if (idx > kPivotNum || tags[idx].first->is_leaf) return -1;
    auto pos = GetNodeIdx(idx << 1, T);
    if (pos != -1) return pos;
    return GetNodeIdx(idx << 1 | 1, T);
  }

  // NOTE: cores
  inline void ResetTagsNum() { tags_num = 0; }

  Box GetBoxByRegionIdx(int const idx, Box const& box) {
    if constexpr (IsBinaryNode<Interior>) {
      assert(GetDepthByIndex(10) == 3);
      Box bx(box);
      int h = GetDepthByIndex(idx);
      for (int i = h - 1, new_idx = 1; i >= 0; i--) {
        int local_id = (idx >> i) & 1;
        auto TI = static_cast<Interior*>(tags[new_idx].first);
        auto& target = local_id ? bx.first : bx.second;
        target.pnt[TI->split.second] = TI->split.first;
        new_idx = new_idx << 1 | local_id;
        assert(new_idx <= idx);
      }
      return std::move(bx);
    } else if constexpr (IsMultiNode<Interior>) {
      assert(GetDepthByIndex(19) == 4);
      assert(kDim ==
             static_cast<BucketType>(std::log2(Interior::GetRegions())));

      Box bx(box);
      BucketType h = GetDepthByIndex(idx);
      for (BucketType i = h, new_idx = 1; i > 0; i -= kDim) {
        BucketType local_id = (idx >> (i - kDim)) & ((1 << kDim) - 1);
        auto TI = static_cast<Interior*>(tags[new_idx].first);
        TI->ModifyBoxById(local_id, bx);
        new_idx = new_idx << kDim | local_id;
      }
      return std::move(bx);
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
  }

  // NOTE: Each Node in the skeleton receives a tag
  // NOTE: A Leaf Node receives the tag < BUCKETNUM
  // NOTE: All internal Node has tag == BUCKETNUM
  void AssignNodeTag(Node* T, BucketType idx) {
    if (T->is_leaf || idx > kPivotNum) {
      assert(tags_num < kBucketNum);
      tags[idx] = NodeTag(T, tags_num);
      rev_tag[tags_num++] = idx;  // WARN: cannot remove
      return;
    }

    // INFO: BUCKET ID in [0, kBucketNum)
    tags[idx] = NodeTag(T, kBucketNum);
    Interior* TI = static_cast<Interior*>(T);
    if constexpr (IsBinaryNode<Interior>) {
      AssignNodeTag(TI->left, idx << 1);
      AssignNodeTag(TI->right, idx << 1 | 1);
    } else if constexpr (IsMultiNode<Interior>) {
      for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
        AssignNodeTag(TI->tree_nodes[i], idx * Interior::GetRegions() + i);
      }
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
    return;
  }

  // NOTE: reduce sums is travsersal of the skeleton with counting the points
  // seieved onto every node, it is good to determine whether we need forcing
  // parallel in the following operations, e.g., flatten/rebuild the tree.
  void ReduceSums(BucketType idx) const {
    if (idx > kPivotNum || tags[idx].first->is_leaf) {
      assert(tags[idx].second < kBucketNum);
      sums_tree[idx] = sums[tags[idx].second];

      // PERF: no need to update the parallel flag here, as it is either a leaf
      // node or it will be handled by recursive calls
      return;
    }

    if constexpr (IsBinaryNode<Interior>) {
      ReduceSums(idx << 1);
      ReduceSums(idx << 1 | 1);
      sums_tree[idx] = sums_tree[idx << 1] + sums_tree[idx << 1 | 1];
    } else if constexpr (IsMultiNode<Interior>) {
      sums_tree[idx] = 0;
      for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
        ReduceSums(idx * Interior::GetRegions() + i);
        sums_tree[idx] += sums_tree[idx * Interior::GetRegions() + i];
      }
    }

    // PERF: Don't add force parallel here as it not precise: whether a node
    // should be rebuilt depends on only the tomb status, rather than the points
    // sieved to that node
    return;
  }

  static inline BucketType GetDepthByIndex(BucketType idx) {
    BucketType h = 0;
    while (idx > 1) {
      idx >>= 1;
      ++h;
    }
    return h;
  }

  // NOTE: a bucket/leaf has id kBucketNum+1
  // a node needs to be rebuilt has id kBucketNum+2
  // otherwise, it has id kBucketNum
  template <typename ViolateFunc>
  void PickTag(BucketType idx, ViolateFunc&& violate_func) {
    if (idx > kPivotNum || tags[idx].first->is_leaf) {
      tags[idx].second = kBucketNum + 1;
      rev_tag[tags_num++] = idx;
      return;
    }

    assert(tags[idx].second == kBucketNum && (!tags[idx].first->is_leaf));
    // if (violate_func(idx)) {
    if (InvokeWithOptionalArg<bool>(violate_func, idx)) {
      tags[idx].second = kBucketNum + 2;
      rev_tag[tags_num++] = idx;
      return;
    }

    if constexpr (IsBinaryNode<Interior>) {
      PickTag(idx << 1, violate_func);
      PickTag(idx << 1 | 1, violate_func);
    } else if constexpr (IsMultiNode<Interior>) {
      for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
        PickTag(idx * Interior::GetRegions() + i, violate_func);
      }
    }
    return;
  }

  template <typename... Args>
  void TagInbalanceNode(Args&&... args) {
    ReduceSums(1);
    ResetTagsNum();
    PickTag(1, std::forward<Args>(args)...);
    assert(AssertSize(tags[1].first));
    return;
  }

  // NOTE: the Node which needs to be rebuilt has tag kBucketNum+3
  // the *bucket* Node whose ancestor has been rebuilt has tag kBucketNum+2
  // the *bucket* Node whose ancestor has not been ... has kBucketNum+1
  // otherwise, it's kBucketNum
  template <bool kSetParallelFlag, typename ViolateFunc>
  void MarkTomb(BucketType idx, BucketType& re_num, size_t& tot_re_size,
                BoxSeq& box_seq, Box const& box, bool has_tomb,
                ViolateFunc&& violate_func) {
    if (idx > kPivotNum || tags[idx].first->is_leaf) {
      assert(tags[idx].second >= 0 && tags[idx].second < kBucketNum);
      if (!has_tomb) {
        tags[idx].second = kBucketNum + 2;
        if constexpr (kSetParallelFlag) {  // INFO: the sub-tree will be rebuilt
                                           // in the future and need to force
                                           // the parallisim
          if (!tags[idx].first->is_leaf) {
            auto TI = static_cast<Interior*>(tags[idx].first);
            TI->SetParallelFlag(TI->size > BT::kSerialBuildCutoff);
          }
        }
      } else {
        tags[idx].second = kBucketNum + 1;
      }

      box_seq[tags_num] = box;
      rev_tag[tags_num++] = idx;
      return;
    }

    // NOTE: no need to mark the internal nodes with tag kBucketNum
    assert(tags[idx].second == kBucketNum && (!tags[idx].first->is_leaf));
    Interior* TI = static_cast<Interior*>(tags[idx].first);
    if (has_tomb && InvokeWithOptionalArg<bool>(violate_func, idx)) {
      tags[idx].second = kBucketNum + 3;
      has_tomb = false;
      TI->SetParallelFlag(TI->size > BT::kSerialBuildCutoff);
      re_num++, tot_re_size += TI->size;
    }

    if constexpr (IsBinaryNode<Interior>) {
      BoxCut box_cut(box, TI->split, true);
      MarkTomb<kSetParallelFlag>(idx << 1, re_num, tot_re_size, box_seq,
                                 box_cut.GetFirstBoxCut(), has_tomb,
                                 violate_func);
      MarkTomb<kSetParallelFlag>(idx << 1 | 1, re_num, tot_re_size, box_seq,
                                 box_cut.GetSecondBoxCut(), has_tomb,
                                 violate_func);
    } else if constexpr (IsMultiNode<Interior>) {
      // BoxSeq new_box(TI->template ComputeSubregions<BoxSeq>(box));
      for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
        MarkTomb<kSetParallelFlag>(
            idx * Interior::GetRegions() + i, re_num, tot_re_size, box_seq,
            TI->GetBoxByRegionId(i, box), has_tomb, violate_func);
      }
    }
    return;
  }

  // NOTE: the Node which needs to be rebuilt has tag kBucketNum+3
  // the *bucket* Node whose ancestor has been rebuilt has tag kBucketNum+2
  // the *bucket* Node whose ancestor has not been ... has kBucketNum+1
  // otherwise, it's kBucketNum

  template <bool kSetParallelFlag, typename... Args>
  auto TagInbalanceNodeDeletion(Args&&... args) {
    ReduceSums(1);
    ResetTagsNum();
    BucketType re_num = 0;
    size_t tot_re_size = 0;
    MarkTomb<kSetParallelFlag>(1, re_num, tot_re_size,
                               std::forward<Args>(args)...);
    return std::make_pair(re_num, tot_re_size);
  }

  void TagPuffyNodesRecursive(BucketType idx) {
    if (idx > kPivotNum || tags[idx].first->is_leaf) {
      tags[idx].second =
          kBucketNum +
          2;  // PERF: ensure the following UpdateInterior meets the base case
      if (!tags[idx].first->is_leaf) {
        Interior* TI = static_cast<Interior*>(tags[idx].first);
        TI->SetParallelFlag(TI->size > BT::kSerialBuildCutoff);
      }
      return;
    }

    assert(!tags[idx].first->is_leaf);
    Interior* TI = static_cast<Interior*>(tags[idx].first);
    TI->SetParallelFlag(TI->size > BT::kSerialBuildCutoff);

    if constexpr (IsBinaryNode<Interior>) {
      TagPuffyNodesRecursive(idx << 1);
      TagPuffyNodesRecursive(idx << 1 | 1);
    } else if constexpr (IsMultiNode<Interior>) {
      for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
        TagPuffyNodesRecursive(idx * Interior::GetRegions() + i);
      }
    }
    return;
  }

  void TagPuffyNodes() {
    TagPuffyNodesRecursive(1);
    return;
  }

  // NOTE: kUpdatePointer: Update the pointer only, if it contains box, return
  //   empty box
  // kUpdatePointerBox: Update pointer and box
  // kTagRebuildNode: Update the pointer and box, meanwhile it assign the
  //   imbalance node with a new tag
  // kPostDelUpdate: Update the skeleton after rebuild (e.g., size, children),
  //   which needs to avoid touch the deleted nodes
  enum UpdateType {
    kUpdatePointer,
    kUpdatePointerBox,
    kTagRebuildNode,
    kPostDelUpdate
  };

  // NOTE: update the skeleton based on the @UpdateType
  template <UpdateType kUT, bool UpdateParFlag = true, typename Base,
            typename... Args>
    requires IsPointerToNode<Base> || IsNodeBox<Base, Point>
  Base UpdateInnerTree(parlay::sequence<Base> const& tree_nodes,
                       Args&&... args) {
    BucketType p = 0;
    if constexpr (kUT == kUpdatePointer ||
                  kUT == kUpdatePointerBox) {  // NOTE: update the inner tree
                                               // nodes or box
      return UpdateInnerTreeRecursive<kUT, UpdateParFlag>(1, tree_nodes, p,
                                                          [&]() {});
    } else if constexpr (kUT == kTagRebuildNode) {  // NOTE: tag the node that
                                                    // needs to be rebuild
      this->ResetTagsNum();

      auto func_2_rebuild_node = [&](auto const&... params) -> void {
        if constexpr (IsBinaryNode<Interior>) {  // NOTE: needs to save the box
          auto const& new_box = std::get<0>(std::forward_as_tuple(params...));
          auto const& idx = std::get<1>(std::forward_as_tuple(params...));
          rev_tag[tags_num] = idx;
          FindVar<BoxSeq>(std::forward<Args>(args)...)[tags_num++] = new_box;
        } else if constexpr (IsMultiNode<Interior>) {  // NOTE: the box is fixed
                                                       // in orth node, no need
                                                       // to save
          auto const& idx = std::get<0>(std::forward_as_tuple(params...));
          rev_tag[tags_num++] = idx;
        }
      };

      return UpdateInnerTreeRecursive<kUT, UpdateParFlag>(1, tree_nodes, p,
                                                          func_2_rebuild_node);
    } else if constexpr (kUT == kPostDelUpdate) {  // NOTE: avoid touch the node
                                                   // that has been deleted
      // PARA: op == 0 -> toggle whether under a rebuild tree
      // op == 1 -> query current status
      bool under_rebuild_tree = false;
      return UpdateInnerTreeRecursive<kUT, UpdateParFlag>(
          1, tree_nodes, p, [&](bool op) -> bool {
            return op == 0 ? (under_rebuild_tree = !under_rebuild_tree)
                           : under_rebuild_tree;
          });
    } else {
      static_assert(0);
    }
  }

  // NOTE: udpate inner tree for binary nodes
  template <UpdateType kUT, bool UpdateParFlag, typename Base, typename Func>
    requires IsBinaryNode<Interior>
  Base UpdateInnerTreeRecursive(BucketType idx,
                                parlay::sequence<Base> const& tree_nodes,
                                BucketType& p, Func&& func) {
    // WARN: needs to ensure this success for both insert and delete
    if (this->tags[idx].second == kBucketNum + 1 ||
        this->tags[idx].second == kBucketNum + 2) {
      return tree_nodes[p++];
    }

    if constexpr (kUT == kPostDelUpdate) {
      if (this->tags[idx].second == kBucketNum + 3) {
        func(0);  // close the under_rebuild_tree flag
        assert(func(1) == true);
      }
    }

    Base const& left = UpdateInnerTreeRecursive<kUT, UpdateParFlag>(
        idx << 1, tree_nodes, p, func);
    Base const& right = UpdateInnerTreeRecursive<kUT, UpdateParFlag>(
        idx << 1 | 1, tree_nodes, p, func);

    if constexpr (kUT == kUpdatePointer) {  // only update the pointers
      UpdateInterior<Interior, UpdateParFlag>(this->tags[idx].first, left,
                                              right);
      if constexpr (IsPointerToNode<Base>) {
        return this->tags[idx].first;
      } else {  // WARN: if only update pointer, then avoid update box
        return NodeBox(this->tags[idx].first, Box());
      }
    } else if constexpr (kUT == kUpdatePointerBox) {  // update pointer and box
      UpdateInterior<Interior, UpdateParFlag>(this->tags[idx].first, left,
                                              right);
      return NodeBox(this->tags[idx].first,
                     BT::GetBox(left.second, right.second));
    } else if constexpr (kUT == kTagRebuildNode) {  // retag and save
                                                    // box for rebuild
      UpdateInterior<Interior, UpdateParFlag>(this->tags[idx].first, left,
                                              right);
      auto new_box = BT::GetBox(left.second, right.second);
      if (this->tags[idx].second == BT::kBucketNum + 3) {
        func(new_box, idx);
      }
      return NodeBox(this->tags[idx].first, std::move(new_box));
    } else if constexpr (kUT == kPostDelUpdate) {  // avoid update pointers
                                                   // for deleted trees
      if (!func(1)) {  // query whether under the rebuild_tree
        UpdateInterior<Interior, UpdateParFlag>(this->tags[idx].first, left,
                                                right);
        return NodeBox(this->tags[idx].first,
                       Box());  // box has been computed before
      } else if (this->tags[idx].second == kBucketNum + 3) {  // recurse back
        func(0);  // disable the under_rebuild_tree flag
        if (!this->tags[idx].first->is_leaf) {
          static_cast<Interior*>(this->tags[idx].first)->ResetParallelFlag();
        }
        assert(func(1) == false);
        return NodeBox(this->tags[idx].first, Box());
      } else {  // the tree has been deleted
        return NodeBox(nullptr, Box());
      }
    } else {
      static_assert(0);
    }
  }

  // NOTE: update inner tree for multi nodes
  template <UpdateType kUT, bool UpdateParFlag, typename Base, typename Func>
    requires IsMultiNode<Interior>
  Base UpdateInnerTreeRecursive(BucketType idx,
                                parlay::sequence<Base> const& tree_nodes,
                                BucketType& p, Func&& func) {
    if (tags[idx].second == BT::kBucketNum + 1 ||
        tags[idx].second == BT::kBucketNum + 2) {
      return tree_nodes[p++];
    }

    if constexpr (kUT == kPostDelUpdate) {
      if (this->tags[idx].second == kBucketNum + 3) {
        func(0);  // close the under_rebuild_tree flag
        assert(func(1) == true);
      }
    }

    typename Interior::OrthNodeArr new_nodes;
    for (BucketType i = 0; i < Interior::GetRegions(); ++i) {
      new_nodes[i] = UpdateInnerTreeRecursive<kUT, UpdateParFlag>(
          idx * Interior::GetRegions() + i, tree_nodes, p, func);
    }

    if constexpr (kUT == kUpdatePointer) {
      BT::template UpdateInterior<Interior, UpdateParFlag>(tags[idx].first,
                                                           new_nodes);
      return this->tags[idx].first;
    } else if constexpr (kUT == kTagRebuildNode) {
      UpdateInterior<Interior, UpdateParFlag>(this->tags[idx].first, new_nodes);
      if (this->tags[idx].second == BT::kBucketNum + 3) {
        func(idx);
      }
      return this->tags[idx].first;
    } else if constexpr (kUT == kPostDelUpdate) {  // NOTE: avoid update
                                                   // pointers for deleted trees
      if (!func(1)) {                              // not under rebuild tree
        UpdateInterior<Interior, UpdateParFlag>(this->tags[idx].first,
                                                new_nodes);
        return this->tags[idx].first;  // box has been computed before
      } else if (this->tags[idx].second == kBucketNum + 3) {  // back
        func(0);  // disable the under_rebuild_tree flag
        assert(func(1) == false);
        if (!this->tags[idx].first->is_leaf) {
          static_cast<Interior*>(this->tags[idx].first)->ResetParallelFlag();
        }
        return this->tags[idx].first;
      } else {  // the tree has been deleted
        return nullptr;
      }
    } else {
      static_assert(false);
    }
  }

  void Reset() {
    ResetTagsNum();
    tags = NodeTagSeq::uninitialized(kPivotNum + kBucketNum + 1);
    sums_tree = parlay::sequence<BallsType>(kPivotNum + kBucketNum + 1);
    rev_tag = BucketSeq::uninitialized(kBucketNum);
  }

  // NOTE: variables
  BT& BTRef;
  BucketType tags_num;
  NodeTagSeq tags;  // PARA: Assign each node a tag, aka skeleton
  mutable parlay::sequence<BallsType> sums_tree;
  mutable BucketSeq rev_tag;  // PARA: maps tag to the position in skeleton
  parlay::sequence<BallsType> sums;
};
};  // namespace pstp

#endif  // PSTP_BASE_TREE_IMPL_INNER_TREE_HPP_
