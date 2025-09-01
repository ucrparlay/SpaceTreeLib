#ifndef PSI_ORTH_TREE_IMPL_ORTH_BATCH_DIFF_HPP_
#define PSI_ORTH_TREE_IMPL_ORTH_BATCH_DIFF_HPP_

#include <tuple>

#include "../orth_tree.h"

namespace psi {

// NOTE: default batch delete
template <typename Point, typename SplitRule, typename LeafAugType,
          typename InteriorAugType, uint_fast8_t kMD, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void OrthTree<Point, SplitRule, LeafAugType, InteriorAugType, kMD, kSkHeight,
              kImbaRatio>::BatchDiff(Range&& In) {
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
template <typename Point, typename SplitRule, typename LeafAugType,
          typename InteriorAugType, uint_fast8_t kMD, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void OrthTree<Point, SplitRule, LeafAugType, InteriorAugType, kMD, kSkHeight,
              kImbaRatio>::BatchDiff_(Slice A) {
  // NOTE: diff points from the tree
  Points B = Points::uninitialized(A.size());
  this->root_ = BatchDiffRecursive(this->root_, A, parlay::make_slice(B));

  // NOTE: launch the rebuild
  // PARA: @prepare_func: function that computes the new parameters before the
  // rebuildtree recursive
  auto prepare_func = [&]([[maybe_unused]] Node* T, [[maybe_unused]] size_t i,
                          Box const& box) {
    auto new_box = static_cast<Interior*>(T)->GetBoxByRegionId(i, box);
    assert(BT::WithinBox(new_box, box));
    return std::make_tuple(std::move(new_box));
  };
  this->root_ = BT::template RebuildTreeRecursive<Leaf, Interior, false>(
      this->root_, prepare_func, split_rule_.AllowRebuild(), this->tree_box_);
  return;
}

// NOTE: the orth does not need box since the box will never change
template <typename Point, typename SplitRule, typename LeafAugType,
          typename InteriorAugType, uint_fast8_t kMD, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
Node* OrthTree<Point, SplitRule, LeafAugType, InteriorAugType, kMD, kSkHeight,
               kImbaRatio>::BatchDiffRecursive(Node* T, Slice In, Slice Out) {
  size_t n = In.size();

  if (n == 0) {
    return T;
  }

  if (T->is_leaf) {
    return BT::template DiffPoints4Leaf<Leaf, Node*>(T, In);
  }

  // if (In.size()) {
  if (In.size() <= BT::kSerialBuildCutoff) {
    parlay::sequence<BallsType> sums(kNodeRegions, 0);
    SerialSplitSkeleton(T, In, 0, 1, sums);
    assert(std::cmp_equal(std::accumulate(sums.begin(), sums.end(), 0), n));

    auto TI = static_cast<Interior*>(T);
    OrthNodeArr new_nodes;

    size_t start = 0;
    for (DimsType i = 0; i < kNodeRegions; ++i) {
      new_nodes[i] =
          BatchDiffRecursive(TI->tree_nodes[i], In.cut(start, start + sums[i]),
                             Out.cut(start, start + sums[i]));
      start += sums[i];
    }

    bool const force_parallel_flag = TI->size > BT::kSerialBuildCutoff;
    BT::template UpdateInterior<Interior>(T, new_nodes);
    assert(T->is_leaf == false);

    if (BT::SparcyNode(0, TI->size)) {
      TI->SetParallelFlag(force_parallel_flag);
    }

    return T;
  }

  InnerTree IT(*this);
  IT.AssignNodeTag(T, 1);
  assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
  BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                      IT.tags_num);
  IT.TagPuffyNodes();

  // PERF: no need to call tag inbalance node here, as the bounding box for
  // orth-tree is fixed

  auto tree_nodes = parlay::sequence<Node*>::uninitialized(IT.tags_num);
  parlay::parallel_for(
      0, IT.tags_num,
      [&](decltype(IT.tags_num) i) {
        size_t start = 0;
        for (decltype(IT.tags_num) j = 0; j < i; j++) {
          start += IT.sums[j];
        }

        tree_nodes[i] = BatchDiffRecursive(IT.tags[IT.rev_tag[i]].first,
                                           Out.cut(start, start + IT.sums[i]),
                                           In.cut(start, start + IT.sums[i]));

        // NOTE: after pick the tag, the tag id is same as the bucket id.
        // In order to match the base case in UpdateInnerTree, we need to
        // manually change it to kBucketNum+2, i.e., none-of its ancestor has
        // been rebuilt
        IT.tags[IT.rev_tag[i]].second = BT::kBucketNum + 2;
      },
      1);

  return IT.template UpdateInnerTree<InnerTree::kUpdatePointer, false>(
      tree_nodes);
}

}  // namespace psi

#endif  // PSI_ORTH_TREE_IMPL_ORTH_BATCH_DIFF_HPP_
