#pragma once

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include <algorithm>

#include "cpdd/dependence/tree_node.h"
#include "cpdd/orth_tree.h"

namespace cpdd {
template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
template <typename Range, typename... Args>
void OrthTree<Point, SplitRule, kMD, kBDO>::Build(Range&& In, Args&&... args) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);
  static_assert(BT::kBuildDepthOnce % kMD == 0);
  assert(kMD == BT::kDim);

  Slice A = parlay::make_slice(In);
  Build_(A, std::forward<Args>(args)...);
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
void OrthTree<Point, SplitRule, kMD, kBDO>::DivideRotate(HyperPlaneSeq& pivots,
                                                         DimsType dim,
                                                         BucketType idx,
                                                         BoxSeq& box_seq,
                                                         Box const& box) {
  if (idx > BT::kPivotNum) {
    // WARN: sometimes cut dimension can be -1
    //  never use pivots[idx].first == -1 to check whether it is in bucket;
    //  instead, use idx > PIVOT_NUM
    box_seq[idx - BT::kBucketNum] = box;
    pivots[idx] = HyperPlane(0, idx - BT::kBucketNum);
    return;
  }

  pivots[idx] = HyperPlane(
      static_cast<Coord>((box.first.pnt[dim] + box.second.pnt[dim]) / 2), dim);

  Box lbox(box), rbox(box);
  lbox.second.pnt[dim] = pivots[idx].first;
  rbox.first.pnt[dim] = pivots[idx].first;

  dim = (dim + 1) % BT::kDim;
  DivideRotate(pivots, dim, 2 * idx, box_seq, lbox);
  DivideRotate(pivots, dim, 2 * idx + 1, box_seq, rbox);

  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
void OrthTree<Point, SplitRule, kMD, kBDO>::SerialSplit(
    Slice In, DimsType dim, DimsType idx, Box const& box, Splitter const& split,
    parlay::sequence<BallsType>& sums, BoxSeq& box_seq) {
  assert(dim <= BT::kDim);

  if (dim == BT::kDim) {
    sums[idx - kNodeRegions] = In.size();
    box_seq[idx - kNodeRegions] = box;
    return;
  }

  auto mid = split[dim].first;
  assert(dim == split[dim].second);

  PointsIter split_iter = std::ranges::partition(In, [&](Point const& p) {
                            return Num::Lt(p.pnt[dim], mid);
                          }).begin();

  Box lbox(box), rbox(box);
  lbox.second.pnt[dim] = mid;
  rbox.first.pnt[dim] = mid;
  SerialSplit(In.cut(0, split_iter - In.begin()), dim + 1, idx << 1, lbox,
              split, sums, box_seq);
  SerialSplit(In.cut(split_iter - In.begin(), In.size()), dim + 1, idx << 1 | 1,
              rbox, split, sums, box_seq);
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
Node* OrthTree<Point, SplitRule, kMD, kBDO>::SerialBuildRecursive(
    Slice In, Slice Out, Box const& box, bool checked_duplicate) {
  assert(In.size() == 0 || BT::WithinBox(BT::GetBox(In), box));
  size_t n = In.size();

  if (n == 0) {
    return AllocEmptyLeafNode<Slice, Leaf>();
  }

  if (n <= BT::kLeaveWrap) {
    return AllocNormalLeafNode<Slice, Leaf>(In);
  }

  assert(kSplitterNum == BT::kDim);

  DimsType dim = 0;
  Splitter split;
  for (DimsType i = 0; i < kSplitterNum; ++i) {
    split[i] = HyperPlane(
        static_cast<Coord>((box.first.pnt[i] + box.second.pnt[i]) / 2), i);
  }

  parlay::sequence<BallsType> sums(kNodeRegions, 0);
  BoxSeq box_seq(kNodeRegions);
  SerialSplit(In, dim, 1, box, split, sums, box_seq);
  assert(std::cmp_equal(std::accumulate(sums.begin(), sums.end(), 0), n));

  if (std::ranges::count(sums, 0) == kNodeRegions - 1) {
    // NOTE: avoid the repeat check as the last
    if (!checked_duplicate) {
      if (std::ranges::find_if_not(In, [&](Point const& p) {
            return p.SameDimension(In[0]);
          }) == In.end()) {
        // WARN: Need to pass full range
        return AllocDummyLeafNode<Slice, Leaf>(In);
      }
      checked_duplicate = true;
    }
  } else {
    checked_duplicate = false;
  }

  OrthNodeArr tree_nodes;
  size_t start = 0;
  for (DimsType i = 0; i < kNodeRegions; ++i) {
    // NOTE: iterate through non-empty partitions, put them into the
    // position identified by non_empty_node
    tree_nodes[i] = SerialBuildRecursive(In.cut(start, start + sums[i]),
                                         Out.cut(start, start + sums[i]),
                                         box_seq[i], checked_duplicate);
    start += sums[i];
  }

  return AllocInteriorNode<Interior>(tree_nodes, split, AugType());
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
Node* OrthTree<Point, SplitRule, kMD, kBDO>::OrthBuildInnerTree(
    BucketType idx, HyperPlaneSeq const& pivots,
    parlay::sequence<Node*> const& tree_nodes) {
  assert(idx < BT::kPivotNum + BT::kBucketNum + 1);

  if (idx > BT::kPivotNum) {
    return tree_nodes[idx - BT::kPivotNum - 1];
  }

  OrthNodeArr multi_nodes;
  Splitter split;
  for (DimsType i = 0; i < kNodeRegions; ++i) {
    multi_nodes[i] =
        OrthBuildInnerTree(idx * kNodeRegions + i, pivots, tree_nodes);
  }
  for (DimsType i = 0; i < kSplitterNum; ++i) {
    split[i] = pivots[idx * (1 << i)];
    assert(i == 0 || pivots[idx * (1 << i)] == pivots[idx * (1 << i) + 1]);
  }

  return AllocInteriorNode<Interior>(multi_nodes, split, AugType());
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
Node* OrthTree<Point, SplitRule, kMD, kBDO>::BuildRecursive(Slice In, Slice Out,
                                                            Box const& box) {
  // TODO: may ensure the bucket is corresponding the the splitter
  assert(In.size() == 0 || BT::WithinBox(BT::GetBox(In), box));

  // if (In.size()) {
  if (In.size() <= BT::kSerialBuildCutoff) {
    return SerialBuildRecursive(In, Out, box, false);
  }

  auto pivots =
      HyperPlaneSeq::uninitialized(BT::kPivotNum + BT::kBucketNum + 1);
  auto box_seq = BoxSeq::uninitialized(BT::kBucketNum);
  parlay::sequence<BallsType> sums;

  DivideRotate(pivots, 0, 1, box_seq, box);
  BT::Partition(In, Out, In.size(), pivots, sums);

  // NOTE: if random sampling failed to split points, re-partitions using
  // serail approach
  auto tree_nodes = parlay::sequence<Node*>::uninitialized(BT::kBucketNum);
  auto nodes_map = BucketSeq::uninitialized(BT::kBucketNum);
  BucketType zeros = 0, cnt = 0;
  for (BucketType i = 0; i < BT::kBucketNum; ++i) {
    if (!sums[i]) {
      ++zeros;
      tree_nodes[i] = AllocEmptyLeafNode<Slice, Leaf>();
    } else {
      nodes_map[cnt++] = i;
    }
  }

  parlay::parallel_for(
      0, BT::kBucketNum - zeros,
      [&](BucketType i) {
        size_t start = 0;
        for (BucketType j = 0; j < nodes_map[i]; ++j) {
          start += sums[j];
        }

        tree_nodes[nodes_map[i]] = BuildRecursive(
            Out.cut(start, start + sums[nodes_map[i]]),
            In.cut(start, start + sums[nodes_map[i]]), box_seq[nodes_map[i]]);
      },
      1);

  return OrthBuildInnerTree(1, pivots, tree_nodes);
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
void OrthTree<Point, SplitRule, kMD, kBDO>::Build_(Slice A) {
  Points B = Points::uninitialized(A.size());
  this->tree_box_ = BT::GetBox(A);
  this->root_ = BuildRecursive(A, B.cut(0, A.size()), this->tree_box_);
  // this->root_ = SerialBuildRecursive(A, B.cut(0, A.size()), BT::kDim,
  //                                    this->tree_box_, false);
  assert(this->root_ != nullptr);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kMD,
          uint_fast8_t kBDO>
void OrthTree<Point, SplitRule, kMD, kBDO>::Build_(Slice A, Box const& box) {
  assert(BT::WithinBox(BT::GetBox(A), box));

  Points B = Points::uninitialized(A.size());
  this->tree_box_ = box;
  this->root_ = BuildRecursive(A, B.cut(0, A.size()), this->tree_box_);
  // this->root_ = SerialBuildRecursive(A, B.cut(0, A.size()), BT::kDim,
  //                                    this->tree_box_, false);
  assert(this->root_ != nullptr);
  return;
}

}  // namespace cpdd
