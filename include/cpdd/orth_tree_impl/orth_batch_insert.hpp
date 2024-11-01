#pragma once

#include <utility>

#include "../orth_tree.h"

namespace cpdd {

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
template <typename Range>
void OrthTree<Point, SplitRule, kMD, kBDO>::BatchInsert(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);
  static_assert(BT::kBuildDepthOnce % kMD == 0);
  assert(kMD == BT::kDim);
  // TODO: handling the case that insert box is no in the tree box
  assert(BT::WithinBox(BT::GetBox(In), this->tree_box_));

  Slice A = parlay::make_slice(In);
  BatchInsert_(A);
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
void OrthTree<Point, SplitRule, kMD, kBDO>::BatchInsert_(Slice A) {
  if (this->root_ == nullptr) {  // TODO: may check using explicity tag
    return Build_(A);
  }

  Points B = Points::uninitialized(A.size());
  Node* T = this->root_;
  // this->tree_box_ = BT::GetBox(this->tree_box_, BT::GetBox(A));
  // PERF: no need to compute bounding box here, checked previously
  this->root_ = BatchInsertRecursive(T, A, B.cut(0, A.size()));
  assert(this->root_ != NULL);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
void OrthTree<Point, SplitRule, kMD, kBDO>::SerialSplitSkeleton(
    Node* T, Slice In, DimsType dim, DimsType idx,
    parlay::sequence<BallsType>& sums) {
  if (dim == BT::kDim) {
    sums[idx - kNodeRegions] = In.size();
    return;
  }

  auto mid = static_cast<Interior*>(T)->split[dim].first;
  assert(dim == static_cast<Interior*>(T)->split[dim].second);

  PointsIter split_iter = std::ranges::partition(In, [&](Point const& p) {
                            return Num::Lt(p.pnt[dim], mid);
                          }).begin();

  SerialSplitSkeleton(T, In.cut(0, split_iter - In.begin()), dim + 1, 2 * idx,
                      sums);
  SerialSplitSkeleton(T, In.cut(split_iter - In.begin(), In.size()), dim + 1,
                      2 * idx + 1, sums);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
Node* OrthTree<Point, SplitRule, kMD, kBDO>::BatchInsertRecursive(Node* T,
                                                                  Slice In,
                                                                  Slice Out) {
  size_t n = In.size();

  if (n == 0) return T;

  if (T->is_leaf) {
    Leaf* TL = static_cast<Leaf*>(T);
    if ((!TL->is_dummy && n + T->size <= BT::kLeaveWrap) ||
        (TL->is_dummy &&
         parlay::all_of(In, [&](Point const& p) { return p == TL->pts[0]; }))) {
      return BT::template InsertPoints2Leaf<Leaf>(T, In);
    } else {
      return BT::template RebuildWithInsert<Leaf, Interior>(T, In);
    }
  }

  // if (n) {
  if (n <= BT::kSerialBuildCutoff) {
    parlay::sequence<BallsType> sums(kNodeRegions, 0);
    SerialSplitSkeleton(T, In, 0, 1, sums);
    assert(std::cmp_equal(std::accumulate(sums.begin(), sums.end(), 0), n));

    OrthNodeArr new_nodes;
    size_t start = 0;
    for (DimsType i = 0; i < kNodeRegions; ++i) {
      new_nodes[i] = BatchInsertRecursive(
          static_cast<Interior*>(T)->tree_nodes[i],
          In.cut(start, start + sums[i]), Out.cut(start, start + sums[i]));
      start += sums[i];
    }
    BT::template UpdateInterior<Interior>(T, new_nodes);
    assert(T->is_leaf == false);
    return T;
  }

  // NOTE: assign each Node a tag
  InnerTree IT(*this);
  // IT.Init();
  // assert(IT.rev_tag.size() == BT::kBucketNum);
  IT.AssignNodeTag(T, 1);
  assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);

  BT::template SeievePoints<Interior>(In, Out, n, IT.tags, IT.sums,
                                      IT.tags_num);

  // NOTE: no need to tag imbalance node in orth tree as it never rebuilds the
  // tree, used to remap the bucket node tag to kBucketNum+1
  // TODO: this may incur extra work e.g., reduce sum
  IT.TagInbalanceNode([]() -> bool { return false; });

  auto tree_nodes = parlay::sequence<Node*>::uninitialized(IT.tags_num);

  parlay::parallel_for(
      0, IT.tags_num,
      [&](decltype(IT.tags_num) i) {
        size_t s = 0;
        for (decltype(IT.tags_num) j = 0; j < i; j++) {
          s += IT.sums_tree[IT.rev_tag[j]];
        }

        assert(IT.tags[IT.rev_tag[i]].second == BT::kBucketNum + 1);
        tree_nodes[i] =
            BatchInsertRecursive(IT.tags[IT.rev_tag[i]].first,
                                 Out.cut(s, s + IT.sums_tree[IT.rev_tag[i]]),
                                 In.cut(s, s + IT.sums_tree[IT.rev_tag[i]]));
      },
      1);

  return IT.template UpdateInnerTree<InnerTree::kUpdatePointer>(tree_nodes);
}
}  // namespace cpdd
