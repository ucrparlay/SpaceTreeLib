#pragma once

#include "../kd_tree.h"

namespace cpdd {

// NOTE: default batch delete
template <typename Point, typename SplitRule, uint_fast8_t kBDO>
template <typename Range>
void KdTree<Point, SplitRule, kBDO>::BatchDelete(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  BatchDelete_(A);
  return;
}

// NOTE: assume all Points are fully covered in the tree
template <typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::BatchDelete_(Slice A) {
  Points B = Points::uninitialized(A.size());
  Node* T = this->root_;
  Box bx = this->tree_box_;
  DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
  std::tie(this->root_, this->tree_box_) =
      BatchDeleteRecursive(T, bx, A, parlay::make_slice(B), d, 1);
  return;
}

// NOTE: the Node which needs to be rebuilt has tag BT::kBucketNum+3
// NOTE: the bucket Node whose ancestor has been rebuilt has tag
// BT::kBucketNum+2
// NOTE: the bucket Node whose ancestor has not been ... has
// BT::kBucketNum+1
// NOTE: otherwise, it's BT::kBucketNum
template <typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::DeleteInnerTree(BucketType idx,
                                                NodeTagSeq const& tags,
                                                NodeBoxSeq& tree_nodes,
                                                BucketType& p, DimsType d) {
  if (tags[idx].second == BT::kBucketNum + 1 ||
      tags[idx].second == BT::kBucketNum + 2) {
    return tree_nodes[p++];
  }

  auto& [L, Lbox] =
      DeleteInnerTree(idx << 1, tags, tree_nodes, p, (d + 1) % BT::kDim);
  auto& [R, Rbox] =
      DeleteInnerTree(idx << 1 | 1, tags, tree_nodes, p, (d + 1) % BT::kDim);

  BT::template UpdateInterior<Interior>(tags[idx].first, L, R);

  // WARN: this blocks the parallelsim as it will rebuild the tree one-by-one
  if (tags[idx].second == BT::kBucketNum + 3) {  // NOTE: launch rebuild
    Interior const* TI = static_cast<Interior*>(tags[idx].first);
    assert(BT::ImbalanceNode(TI->left->size, TI->size) ||
           TI->size < BT::kThinLeaveWrap);

    if (tags[idx].first->size == 0) {  // NOTE: special judge for empty tree
      BT::template DeleteTreeRecursive<Leaf, Interior, false>(tags[idx].first);
      return NodeBox(AllocEmptyLeafNode<Slice, Leaf>(), BT::GetEmptyBox());
    }

    return BT::template RebuildSingleTree<Leaf, Interior, false>(
        tags[idx].first, d, BT::GetBox(Lbox, Rbox));
  }

  return NodeBox(tags[idx].first, BT::GetBox(Lbox, Rbox));
}

// NOTE: delete with rebuild, with the assumption that all Points are in the
// tree
// WARN: the param d can be only used when rotate cutting is applied
template <typename Point, typename SplitRule, uint_fast8_t kBDO>
typename KdTree<Point, SplitRule, kBDO>::NodeBox
KdTree<Point, SplitRule, kBDO>::BatchDeleteRecursive(
    Node* T, typename KdTree<Point, SplitRule, kBDO>::Box const& bx, Slice In,
    Slice Out, DimsType d, bool has_tomb) {
  size_t n = In.size();

  if (n == 0) {
    assert(BT::WithinBox(BT::template GetBox<Leaf, Interior>(T), bx));
    return NodeBox(T, bx);
  }

  // INFO: may can be used to accelerate the whole deletion process
  if (n == T->size) {
    if (has_tomb) {  // rebuild this subtree
      BT::template DeleteTreeRecursive<Leaf, Interior>(T);
      return NodeBox(AllocEmptyLeafNode<Slice, Leaf>(), BT::GetEmptyBox());
    }
    // wihtin a rebuild tree
    if (!T->is_leaf) {
      auto TI = static_cast<Interior*>(T);
      // WARN: only set the flag for root ,rather than the whole sub-tree
      TI->SetParallelFlag(T->size > BT::kSerialBuildCutoff);
    }
    T->size = 0;
    return NodeBox(T, bx);
  }

  if (T->is_leaf) {
    return BT::template DeletePoints4Leaf<Leaf, NodeBox>(T, In);
  }

  // if (1) {
  if (In.size() <= BT::kSerialBuildCutoff) {
    Interior* TI = static_cast<Interior*>(T);
    PointsIter split_iter =
        std::ranges::partition(In, [&](Point const& p) {
          return Num::Lt(p.pnt[TI->split.second], TI->split.first);
        }).begin();

    bool putTomb = has_tomb && (BT::ImbalanceNode(
                                    TI->left->size - (split_iter - In.begin()),
                                    TI->size - In.size()) ||
                                TI->size - In.size() < BT::kThinLeaveWrap);
    has_tomb = putTomb ? false : has_tomb;
    assert(putTomb ? (!has_tomb) : true);

    DimsType nextDim = (d + 1) % BT::kDim;
    BoxCut box_cut(bx, TI->split, true);

    auto [L, Lbox] = BatchDeleteRecursive(
        TI->left, box_cut.GetFirstBoxCut(), In.cut(0, split_iter - In.begin()),
        Out.cut(0, split_iter - In.begin()), nextDim, has_tomb);
    auto [R, Rbox] = BatchDeleteRecursive(TI->right, box_cut.GetSecondBoxCut(),
                                          In.cut(split_iter - In.begin(), n),
                                          Out.cut(split_iter - In.begin(), n),
                                          nextDim, has_tomb);

    bool par_flag = TI->size > BT::kSerialBuildCutoff;
    BT::template UpdateInterior<Interior>(T, L, R);
    if (!has_tomb) {  // WARN: Above update will reset parallel flag
      TI->SetParallelFlag(par_flag);
    }

    assert(T->size == L->size + R->size && TI->split.second >= 0 &&
           TI->is_leaf == false);

    // NOTE: rebuild
    if (putTomb) {
      assert(BT::ImbalanceNode(TI->left->size, TI->size) ||
             TI->size < BT::kThinLeaveWrap);

      auto const new_box = BT::GetBox(Lbox, Rbox);
      assert(BT::WithinBox(BT::template GetBox<Leaf, Interior>(T), new_box));

      return NodeBox(
          BT::template RebuildSingleTree<Leaf, Interior, false>(T, d, new_box),
          new_box);
    }

    return NodeBox(T, BT::GetBox(Lbox, Rbox));
  }

  InnerTree IT(*this);
  IT.AssignNodeTag(T, 1);
  assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
  BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                      IT.tags_num);

  auto tree_nodes = NodeBoxSeq::uninitialized(IT.tags_num);
  auto box_seq = parlay::sequence<Box>::uninitialized(IT.tags_num);

  auto [re_num, tot_re_size] = IT.TagInbalanceNodeDeletion(
      box_seq, bx, has_tomb, [&](BucketType idx) -> bool {
        Interior* TI = static_cast<Interior*>(IT.tags[idx].first);
        return BT::ImbalanceNode(TI->left->size - IT.sums_tree[idx << 1],
                                 TI->size - IT.sums_tree[idx]) ||
               (TI->size - IT.sums_tree[idx] < BT::kThinLeaveWrap);
      });

  assert(re_num <= IT.tags_num);

  parlay::parallel_for(
      0, IT.tags_num,
      [&](decltype(IT.tags_num) i) {
        size_t start = 0;
        for (decltype(IT.tags_num) j = 0; j < i; j++) {
          start += IT.sums[j];
        }

        assert(IT.sums_tree[IT.rev_tag[i]] == IT.sums[i]);
        assert(IT.tags[IT.rev_tag[i]].first->size >= IT.sums[i]);
        assert(BT::WithinBox(
            BT::GetBox(Out.cut(start, start + IT.sums[i])),
            BT::template GetBox<Leaf, Interior>(IT.tags[IT.rev_tag[i]].first)));

        DimsType nextDim = (d + IT.GetDepthByIndex(IT.rev_tag[i])) % BT::kDim;

        tree_nodes[i] = BatchDeleteRecursive(
            IT.tags[IT.rev_tag[i]].first, box_seq[i],
            Out.cut(start, start + IT.sums[i]),
            In.cut(start, start + IT.sums[i]), nextDim,
            IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1);
      },
      1);

  // NOTE: handling of rebuild
  // NOTE: get new box for skeleton root and rebuild nodes
  Box const new_box =
      std::get<1>(IT.template UpdateInnerTree<InnerTree::kTagRebuildNode>(
          tree_nodes, box_seq));
  assert(IT.tags_num == re_num);

  parlay::parallel_for(0, IT.tags_num, [&](size_t i) {
    assert(IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 3);

    if (IT.tags[IT.rev_tag[i]].first->size == 0) {  // NOTE: empty
      BT::template DeleteTreeRecursive<Leaf, Interior, false>(
          IT.tags[IT.rev_tag[i]].first);
      IT.tags[IT.rev_tag[i]].first = AllocEmptyLeafNode<Slice, Leaf>();
    } else {  // NOTE: rebuild
      assert(BT::WithinBox(
          BT::template GetBox<Leaf, Interior>(IT.tags[IT.rev_tag[i]].first),
          box_seq[i]));

      auto next_dim = (d + IT.GetDepthByIndex(IT.rev_tag[i])) % BT::kDim;
      IT.tags[IT.rev_tag[i]].first =
          BT::template RebuildSingleTree<Leaf, Interior, false>(
              IT.tags[IT.rev_tag[i]].first, next_dim, box_seq[i]);
    }
  });

  // PARA: op == 0 -> toggle whether under a rebuild tree
  // op == 1 -> query current status
  auto const new_root = std::get<0>(
      IT.template UpdateInnerTree<InnerTree::kPostRebuild>(tree_nodes));
  return NodeBox(new_root, new_box);
}

}  // namespace cpdd
