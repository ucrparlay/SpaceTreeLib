#ifndef PSI_POINTER_BASED_BASE_TREE_IMPL_INNER_TREE_INNER_TREE_HPP_
#define PSI_POINTER_BASED_BASE_TREE_IMPL_INNER_TREE_INNER_TREE_HPP_

#include <utility>

#include "dependence/concepts.h"
#include "inner_tree_binary.hpp"
#include "inner_tree_multi.hpp"
#include "pointer_based/base_tree.h"

// Macro to simplify BinaryNodeOps template usage
#define BINARY_NODE_OPS \
  inner_tree_detail::BinaryNodeOps<TypeTrait, DerivedTree, Leaf, Interior>

// Macro to simplify MultiNodeOps template usage
#define MULTI_NODE_OPS \
  inner_tree_detail::MultiNodeOps<TypeTrait, DerivedTree, Leaf, Interior>

namespace psi {
template <class TypeTrait, typename DerivedTree>
template <typename Leaf, typename Interior>
struct BaseTree<TypeTrait, DerivedTree>::InnerTree {
  using BT = BaseTree<TypeTrait, DerivedTree>;
  using Geo = GeoBase<TypeTrait>;

  InnerTree()
      : tags_num(0),
        tags(NodeTagSeq::uninitialized(kPivotNum + kBucketNum + 1)),
        sums_tree(parlay::sequence<BallsType>(kPivotNum + kBucketNum + 1)),
        rev_tag(BucketSeq::uninitialized(kBucketNum)) {}

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
    if (idx > kPivotNum || tags[idx].first->is_leaf)
      return inner_tree_detail::kInvalidBucket;
    auto pos = GetNodeIdx(idx << 1, T);
    if (pos != inner_tree_detail::kInvalidBucket) return pos;
    return GetNodeIdx(idx << 1 | 1, T);
  }

  inline void ResetTagsNum() { tags_num = 0; }

  Box GetBoxByRegionIdx(int const idx, Box const& box) {
    if constexpr (IsBinaryNode<Interior>) {
      return BINARY_NODE_OPS::GetBoxByRegionIdx(this, idx, box);
    } else if constexpr (IsMultiNode<Interior>) {
      return MULTI_NODE_OPS::GetBoxByRegionIdx(this, idx, box);
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
  }

  void AssignNodeTag(Node* T, BucketType idx) {
    if constexpr (IsBinaryNode<Interior>) {
      BINARY_NODE_OPS::AssignNodeTag(this, T, idx);
    } else if constexpr (IsMultiNode<Interior>) {
      MULTI_NODE_OPS::AssignNodeTag(this, T, idx);
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
  }

  void ReduceSums(BucketType idx) const {
    if constexpr (IsBinaryNode<Interior>) {
      BINARY_NODE_OPS::ReduceSums(const_cast<InnerTree*>(this), idx);
    } else if constexpr (IsMultiNode<Interior>) {
      MULTI_NODE_OPS::ReduceSums(const_cast<InnerTree*>(this), idx);
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
  }

  static inline BucketType GetDepthByIndex(BucketType idx) {
    BucketType h = 0;
    while (idx > 1) {
      idx >>= 1;
      ++h;
    }
    return h;
  }

  template <typename ViolateFunc>
  void PickTag(BucketType idx, ViolateFunc&& violate_func) {
    if constexpr (IsBinaryNode<Interior>) {
      BINARY_NODE_OPS::PickTag(this, idx,
                               std::forward<ViolateFunc>(violate_func));
    } else if constexpr (IsMultiNode<Interior>) {
      MULTI_NODE_OPS::PickTag(this, idx,
                              std::forward<ViolateFunc>(violate_func));
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
  }

  template <typename... Args>
  void TagInbalanceNode(Args&&... args) {
    ReduceSums(1);
    ResetTagsNum();
    PickTag(1, std::forward<Args>(args)...);
    assert(AssertSize(tags[1].first));
  }

  template <bool kSetParallelFlag, typename ViolateFunc>
  void MarkNodesForRebuild(BucketType idx, BucketType& re_num,
                           size_t& tot_re_size, BoxSeq& box_seq, Box const& box,
                           bool has_tomb, ViolateFunc&& violate_func) {
    if constexpr (IsBinaryNode<Interior>) {
      BINARY_NODE_OPS::template MarkNodesForRebuild<kSetParallelFlag>(
          this, idx, re_num, tot_re_size, box_seq, box, has_tomb,
          std::forward<ViolateFunc>(violate_func));
    } else if constexpr (IsMultiNode<Interior>) {
      MULTI_NODE_OPS::template MarkNodesForRebuild<kSetParallelFlag>(
          this, idx, re_num, tot_re_size, box_seq, box, has_tomb,
          std::forward<ViolateFunc>(violate_func));
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
  }

  template <bool kSetParallelFlag, typename... Args>
  auto TagInbalanceNodeDeletion(Args&&... args) {
    ReduceSums(1);
    ResetTagsNum();
    BucketType re_num = 0;
    size_t tot_re_size = 0;
    MarkNodesForRebuild<kSetParallelFlag>(1, re_num, tot_re_size,
                                          std::forward<Args>(args)...);
    return std::make_pair(re_num, tot_re_size);
  }

  void TagOversizedNodesRecursive(BucketType idx) {
    if constexpr (IsBinaryNode<Interior>) {
      BINARY_NODE_OPS::TagOversizedNodesRecursive(this, idx);
    } else if constexpr (IsMultiNode<Interior>) {
      MULTI_NODE_OPS::TagOversizedNodesRecursive(this, idx);
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
  }

  void TagOversizedNodes() { TagOversizedNodesRecursive(1); }

  template <bool UpdateParFlag = true, typename NodeOrNodeBox>
    requires IsPointerToNode<NodeOrNodeBox> || IsNodeBox<NodeOrNodeBox, Point>
  NodeOrNodeBox UpdateInnerTreePointers(
      parlay::sequence<NodeOrNodeBox> const& tree_nodes) {
    BucketType p = 0;
    return UpdateInnerTreePointersRecursive<UpdateParFlag>(1, tree_nodes, p);
  }

  template <bool UpdateParFlag, typename NodeOrNodeBox>
  NodeOrNodeBox UpdateInnerTreePointersRecursive(
      BucketType idx, parlay::sequence<NodeOrNodeBox> const& tree_nodes,
      BucketType& p) {
    if constexpr (IsBinaryNode<Interior>) {
      return BINARY_NODE_OPS::template UpdateInnerTreePointersRecursive<
          UpdateParFlag>(this, idx, tree_nodes, p);
    } else if constexpr (IsMultiNode<Interior>) {
      return MULTI_NODE_OPS::template UpdateInnerTreePointersRecursive<
          UpdateParFlag>(this, idx, tree_nodes, p);
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
  }

  template <bool UpdateParFlag = true, typename NodeOrNodeBox>
    requires IsNodeBox<NodeOrNodeBox, Point>
  NodeOrNodeBox UpdateInnerTreePointersWithBox(
      parlay::sequence<NodeOrNodeBox> const& tree_nodes) {
    BucketType p = 0;
    return UpdateInnerTreePointersWithBoxRecursive<UpdateParFlag>(1, tree_nodes,
                                                                  p);
  }

  template <bool UpdateParFlag, typename NodeOrNodeBox>
  NodeOrNodeBox UpdateInnerTreePointersWithBoxRecursive(
      BucketType idx, parlay::sequence<NodeOrNodeBox> const& tree_nodes,
      BucketType& p) {
    if constexpr (IsBinaryNode<Interior>) {
      return BINARY_NODE_OPS::template UpdateInnerTreePointersWithBoxRecursive<
          UpdateParFlag>(this, idx, tree_nodes, p);
    } else if constexpr (IsMultiNode<Interior>) {
      return MULTI_NODE_OPS::template UpdateInnerTreePointersWithBoxRecursive<
          UpdateParFlag>(this, idx, tree_nodes, p);
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
  }

  template <bool UpdateParFlag = true, typename NodeOrNodeBox, typename... Args>
    requires IsPointerToNode<NodeOrNodeBox> || IsNodeBox<NodeOrNodeBox, Point>
  NodeOrNodeBox TagNodesForRebuild(
      parlay::sequence<NodeOrNodeBox> const& tree_nodes, Args&&... args) {
    this->ResetTagsNum();

    auto func_2_rebuild_node = [&](auto const&... params) -> void {
      if constexpr (IsBinaryNode<Interior>) {
        auto const& new_box = std::get<0>(std::forward_as_tuple(params...));
        auto const& idx = std::get<1>(std::forward_as_tuple(params...));
        rev_tag[tags_num] = idx;
        FindVar<BoxSeq>(std::forward<Args>(args)...)[tags_num++] = new_box;
      } else if constexpr (IsMultiNode<Interior>) {
        auto const& idx = std::get<0>(std::forward_as_tuple(params...));
        rev_tag[tags_num++] = idx;
      }
    };

    BucketType p = 0;
    return TagNodesForRebuildRecursive<UpdateParFlag>(1, tree_nodes, p,
                                                      func_2_rebuild_node);
  }

  template <bool UpdateParFlag, typename NodeOrNodeBox, typename Func>
  NodeOrNodeBox TagNodesForRebuildRecursive(
      BucketType idx, parlay::sequence<NodeOrNodeBox> const& tree_nodes,
      BucketType& p, Func&& func) {
    if constexpr (IsBinaryNode<Interior>) {
      return BINARY_NODE_OPS::template TagNodesForRebuildRecursive<
          UpdateParFlag>(this, idx, tree_nodes, p, std::forward<Func>(func));
    } else if constexpr (IsMultiNode<Interior>) {
      return MULTI_NODE_OPS::template TagNodesForRebuildRecursive<
          UpdateParFlag>(this, idx, tree_nodes, p, std::forward<Func>(func));
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
  }

  template <bool UpdateParFlag = true, typename NodeOrNodeBox>
    requires IsPointerToNode<NodeOrNodeBox> || IsNodeBox<NodeOrNodeBox, Point>
  NodeOrNodeBox UpdateAfterDeletion(
      parlay::sequence<NodeOrNodeBox> const& tree_nodes) {
    bool under_rebuild_tree = false;
    BucketType p = 0;
    return UpdateAfterDeletionRecursive<UpdateParFlag>(
        1, tree_nodes, p, [&](bool op) -> bool {
          return op == 0 ? (under_rebuild_tree = !under_rebuild_tree)
                         : under_rebuild_tree;
        });
  }

  template <bool UpdateParFlag, typename NodeOrNodeBox, typename Func>
  NodeOrNodeBox UpdateAfterDeletionRecursive(
      BucketType idx, parlay::sequence<NodeOrNodeBox> const& tree_nodes,
      BucketType& p, Func&& func) {
    if constexpr (IsBinaryNode<Interior>) {
      return BINARY_NODE_OPS::template UpdateAfterDeletionRecursive<
          UpdateParFlag>(this, idx, tree_nodes, p, std::forward<Func>(func));
    } else if constexpr (IsMultiNode<Interior>) {
      return MULTI_NODE_OPS::template UpdateAfterDeletionRecursive<
          UpdateParFlag>(this, idx, tree_nodes, p, std::forward<Func>(func));
    } else {
      static_assert(IsBinaryNode<Interior> || IsMultiNode<Interior>);
    }
  }

  void Reset() {
    ResetTagsNum();
    tags = NodeTagSeq::uninitialized(kPivotNum + kBucketNum + 1);
    sums_tree = parlay::sequence<BallsType>(kPivotNum + kBucketNum + 1);
    rev_tag = BucketSeq::uninitialized(kBucketNum);
  }

  BucketType tags_num;
  NodeTagSeq tags;
  mutable parlay::sequence<BallsType> sums_tree;
  mutable BucketSeq rev_tag;
  parlay::sequence<BallsType> sums;
};

// Clean up macros
#undef BINARY_NODE_OPS
#undef MULTI_NODE_OPS

};  // namespace psi

#endif  // PSI_POINTER_BASED_BASE_TREE_IMPL_INNER_TREE_INNER_TREE_HPP_