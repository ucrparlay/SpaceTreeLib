#pragma once

#include "../orth_tree.h"

namespace cpdd {

// NOTE: default batch delete
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
template <typename Range>
void OrthTree<Point, SplitRule, kMD, kBDO>::BatchDiff(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  BatchDiff_(A);
  return;
}

// NOTE: assume all Points are fully covered in the tree
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
void OrthTree<Point, SplitRule, kMD, kBDO>::BatchDiff_(Slice A) {
  Points B = Points::uninitialized(A.size());
  this->root_ = BatchDiffRecursive(this->root_, A, parlay::make_slice(B),
                                   this->tree_box_, 1);
  return;
}

// NOTE: delete with rebuild, with the assumption that all Points are in the
// tree
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
Node* OrthTree<Point, SplitRule, kMD, kBDO>::BatchDiffRecursive(
    Node* T, Slice In, Slice Out, Box const& box) {
  size_t n = In.size();

  if (n == 0) {
    return T;
  }

  if (T->is_leaf) {
    return BT::template DiffPoints4Leaf<Leaf, Node*>(T, In);
  }

  if (In.size() <= BT::kSerialBuildCutoff) {
    parlay::sequence<BallsType> sums(kNodeRegions, 0);
    SerialSplitSkeleton(T, In, 0, 1, sums);
    assert(std::cmp_equal(std::accumulate(sums.begin(), sums.end(), 0), n));

    auto TI = static_cast<Interior*>(T);
    OrthNodeArr new_nodes;
    BoxSeq new_box(TI->template ComputeSubregions<BoxSeq>(box));

    size_t start = 0;
    for (DimsType i = 0; i < kNodeRegions; ++i) {
      new_nodes[i] =
          BatchDiffRecursive(TI->tree_nodes[i], In.cut(start, start + sums[i]),
                             Out.cut(start, start + sums[i]), new_box[i]);
      start += sums[i];
    }

    // bool par_flag = TI->size > BT::kSerialBuildCutoff;
    BT::template UpdateInterior<Interior>(T, new_nodes);
    assert(T->is_leaf == false);

    return T;
  }

  InnerTree IT(*this);
  IT.AssignNodeTag(T, 1);
  assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
  BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                      IT.tags_num);

  BoxSeq box_seq(IT.tags_num);  // PARA: the box for bucket nodes
  auto [re_num, tot_re_size] = IT.TagInbalanceNodeDeletion(
      box_seq, box, false, [&]() -> bool { return false; });

  assert(re_num <= IT.tags_num);

  auto tree_nodes = parlay::sequence<Node*>::uninitialized(IT.tags_num);
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

        tree_nodes[i] = BatchDiffRecursive(
            IT.tags[IT.rev_tag[i]].first, Out.cut(start, start + IT.sums[i]),
            In.cut(start, start + IT.sums[i]), box_seq[i],
            IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1);
      },
      1);

  // NOTE: handling of rebuild
  // WARN: the rebuild node is on top
  // NOTE: retag the inba-nodes and save the bounding boxes
  // TODO: here
  [[maybe_unused]] Node* new_node =
      IT.template UpdateInnerTree<InnerTree::kTagRebuildNode>(tree_nodes);
  assert(IT.tags_num == re_num);

  return IT.template UpdateInnerTree<InnerTree::kPostDelUpdate>(tree_nodes);
}

}  // namespace cpdd
