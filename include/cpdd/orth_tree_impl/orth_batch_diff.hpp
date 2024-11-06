#pragma once

#include <tuple>

#include "../orth_tree.h"

namespace cpdd {

// NOTE: default batch delete
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
template <typename Range>
void OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::BatchDiff(
    Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  BatchDiff_(A);
  return;
}

// NOTE: assume points are partially covered in the tree
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
void OrthTree<Point, SplitRule, kMD, kSkHeight, kImbaRatio>::BatchDiff_(
    Slice A) {
  // NOTE: diff points from the tree
  Points B = Points::uninitialized(A.size());
  this->root_ = BatchDiffRecursive(this->root_, A, parlay::make_slice(B));

  // NOTE: launch the rebuild
  auto prepare_func = [&](Node* T, size_t i, Box const& box) {
    auto new_box = static_cast<Interior*>(T)->GetBoxByRegionId(i, box);
    assert(BT::WithinBox(new_box, box));
    return std::make_tuple(std::move(new_box));
  };
  this->root_ = BT::template RebuildTreeRecursive<Leaf, Interior, true>(
      this->root_, prepare_func, this->tree_box_);
  return;
}

// NOTE: the orth does not need box since the box will never change
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kSkHeight, uint_fast8_t kImbaRatio>
Node* OrthTree<Point, SplitRule, kMD, kSkHeight,
               kImbaRatio>::BatchDiffRecursive(Node* T, Slice In, Slice Out) {
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
    // BoxSeq new_box(TI->template ComputeSubregions<BoxSeq>(box));

    size_t start = 0;
    for (DimsType i = 0; i < kNodeRegions; ++i) {
      new_nodes[i] =
          BatchDiffRecursive(TI->tree_nodes[i], In.cut(start, start + sums[i]),
                             Out.cut(start, start + sums[i]));
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

  // BoxSeq box_seq(IT.tags_num);  // PARA: the box for bucket nodes
  // [[maybe_unused]] auto [re_num, tot_re_size] = IT.TagInbalanceNodeDeletion(
  //     box_seq, box, false, [&]() -> bool { return false; });

  // assert(re_num == 0 && tot_re_size == 0);

  auto tree_nodes = parlay::sequence<Node*>::uninitialized(IT.tags_num);
  parlay::parallel_for(
      0, IT.tags_num,
      [&](decltype(IT.tags_num) i) {
        size_t start = 0;
        for (decltype(IT.tags_num) j = 0; j < i; j++) {
          start += IT.sums[j];
        }

        assert(IT.sums_tree[IT.rev_tag[i]] == IT.sums[i]);

        tree_nodes[i] = BatchDiffRecursive(IT.tags[IT.rev_tag[i]].first,
                                           Out.cut(start, start + IT.sums[i]),
                                           In.cut(start, start + IT.sums[i]));
      },
      1);

  return IT.template UpdateInnerTree<InnerTree::kUpdatePointer>(tree_nodes);
}

}  // namespace cpdd
