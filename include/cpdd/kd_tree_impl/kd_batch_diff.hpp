#pragma once
#include "../kd_tree.h"

namespace cpdd {
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDiff(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  BatchDiff_(A);
  return;
}

// NOTE: batch delete suitable for Points that are pratially covered in the tree
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDiff_(Slice A) {
  Points B = Points::uninitialized(A.size());
  Node* T = this->root_;
  Box box = this->tree_box_;

  // NOTE: diff poitns from the tree
  DimsType d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
  std::tie(T, this->tree_box_) =
      BatchDiffRecursive(T, box, A, parlay::make_slice(B), d);

  // NOTE: launch rebuild
  d = T->is_leaf ? 0 : static_cast<Interior*>(T)->split.second;
  auto prepare_rebuild_func = [&](Node* T, DimsType d, Box const& box) {
    DimsType new_dim = (d + 1) % BT::kDim;
    BoxCut box_cut(box, static_cast<Interior*>(T)->split, true);
    auto left_args = std::make_pair(new_dim, box_cut.GetFirstBoxCut());
    auto right_args = std::make_pair(new_dim, box_cut.GetSecondBoxCut());
    return std::make_pair(std::move(left_args), std::move(right_args));
  };
  this->root_ = BT::template RebuildTreeRecursive<Leaf, Interior>(
      T, prepare_rebuild_func, d, this->tree_box_);

  return;
}

// NOTE: only sieve the Points, without rebuilding the tree
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
typename KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::NodeBox
KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::BatchDiffRecursive(
    Node* T,
    typename KdTree<Point, SplitRule, kSkHeight, kImbaRatio>::Box const& box,
    Slice In, Slice Out, DimsType d) {
  size_t n = In.size();

  if (n == 0) return NodeBox(T, box);

  if (T->is_leaf) {
    return BT::template DiffPoints4Leaf<Leaf, NodeBox>(T, In);
  }

  if (In.size() <= BT::kSerialBuildCutoff) {
    // if (In.size()) {
    Interior* TI = static_cast<Interior*>(T);
    PointsIter split_iter =
        std::ranges::partition(In, [&](Point const& p) {
          return Num::Lt(p.pnt[TI->split.second], TI->split.first);
        }).begin();

    DimsType nextDim = (d + 1) % BT::kDim;

    BoxCut box_cut(box, TI->split, true);
    auto [L, Lbox] = BatchDiffRecursive(
        TI->left, box_cut.GetFirstBoxCut(), In.cut(0, split_iter - In.begin()),
        Out.cut(0, split_iter - In.begin()), nextDim);
    auto [R, Rbox] =
        BatchDiffRecursive(TI->right, box_cut.GetSecondBoxCut(),
                           In.cut(split_iter - In.begin(), n),
                           Out.cut(split_iter - In.begin(), n), nextDim);

    BT::template UpdateInterior<Interior>(T, L, R);
    assert(T->size == L->size + R->size && TI->split.second >= 0 &&
           TI->is_leaf == false);

    return NodeBox(T, BT::GetBox(Lbox, Rbox));
  }

  InnerTree IT(*this);
  IT.AssignNodeTag(T, 1);
  assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
  BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                      IT.tags_num);

  auto tree_nodes = NodeBoxSeq::uninitialized(IT.tags_num);
  auto box_seq = parlay::sequence<Box>::uninitialized(IT.tags_num);

  // NOTE: never set tomb, this equivalent to only compute the bounding box,
  [[maybe_unused]] auto [re_num, tot_re_size] =
      IT.TagInbalanceNodeDeletion(box_seq, box, false, []() { return false; });
  assert(re_num == 0 && tot_re_size == 0);

  parlay::parallel_for(
      0, IT.tags_num,
      // NOTE: i is the index of the tags
      [&](BucketType i) {
        size_t start = 0;
        for (BucketType j = 0; j < i; j++) {
          // NOTE: should have same effect as using sums_tree
          // if using sums_tree then it should be sums_tree[rev_tag[j]]
          start += IT.sums[j];
        }

        DimsType nextDim = (d + IT.GetDepthByIndex(IT.rev_tag[i])) % BT::kDim;
        tree_nodes[i] =
            BatchDiffRecursive(IT.tags[IT.rev_tag[i]].first, box_seq[i],
                               Out.cut(start, start + IT.sums[i]),
                               In.cut(start, start + IT.sums[i]), nextDim);
      },
      1);

  return IT.template UpdateInnerTree<InnerTree::kUpdatePointerBox>(tree_nodes);
}

}  // namespace cpdd
