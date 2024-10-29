#pragma once

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include "../kd_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template <typename Point, typename SplitRule, uint_fast8_t kBDO>
template <typename Range>
void KdTree<Point, SplitRule, kBDO>::Build(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  Build_(A);
}

template <typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::DivideRotate(Slice In, SplitterSeq& pivots,
                                                  DimsType dim, BucketType idx,
                                                  BoxSeq& boxs, Box const& bx) {
  if (idx > BT::kPivotNum) {
    // WARN: sometimes cut dimension can be -1
    //  never use pivots[idx].first to check whether it is in bucket;
    //  instead, use idx > PIVOT_NUM
    boxs[idx - BT::kBucketNum] = bx;
    pivots[idx] = Splitter(0, idx - BT::kBucketNum);
    return;
  }
  size_t n = In.size();
  DimsType d = split_rule_.FindCuttingDimension(bx, dim);
  assert(d < BT::kDim);

  std::ranges::nth_element(In, In.begin() + n / 2,
                           [&](Point const& p1, Point const& p2) {
                             return Num::Lt(p1.pnt[d], p2.pnt[d]);
                           });

  pivots[idx] = Splitter(In[n / 2].pnt[d], d);

  Box lbox(bx), rbox(bx);
  lbox.second.pnt[d] = pivots[idx].first;  // PERF: loose
  rbox.first.pnt[d] = pivots[idx].first;

  d = (d + 1) % BT::kDim;
  DivideRotate(In.cut(0, n / 2), pivots, d, 2 * idx, boxs, lbox);
  DivideRotate(In.cut(n / 2, n), pivots, d, 2 * idx + 1, boxs, rbox);
  return;
}

// NOTE: starting at dimesion dim and pick pivots in a rotation manner
template <typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::PickPivots(Slice In, size_t const& n,
                                                SplitterSeq& pivots,
                                                DimsType const dim,
                                                BoxSeq& boxs, Box const& bx) {
  size_t size = std::min(n, static_cast<size_t>(32 * BT::kBucketNum));
  assert(size <= n);

  Points arr = Points::uninitialized(size);
  BT::SamplePoints(In, arr);

  // NOTE: pick pivots
  DivideRotate(arr.cut(0, size), pivots, dim, 1, boxs, bx);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kBDO>
Node* KdTree<Point, SplitRule, kBDO>::SerialBuildRecursive(Slice In, Slice Out,
                                                           DimsType dim,
                                                           Box const& bx) {
  size_t n = In.size();

  if (n == 0) return AllocEmptyLeafNode<Slice, Leaf>();

  if (n <= BT::kLeaveWrap) return AllocNormalLeafNode<Slice, Leaf>(In);

  DimsType d = split_rule_.FindCuttingDimension(bx, dim);
  PointsIter splitIter = BT::SerialPartition(In, d);

  Splitter split;

  if (splitIter <= In.begin() + n / 2) {  // NOTE: split is on left half
    split = Splitter(In[n / 2].pnt[d], d);
  } else if (splitIter != In.end()) {  // NOTE: split is on right half
    auto minEleIter = std::ranges::min_element(
        splitIter, In.end(), [&](Point const& p1, Point const& p2) {
          return Num::Lt(p1.pnt[d], p2.pnt[d]);
        });
    split = Splitter(minEleIter->pnt[d], d);
  } else if (In.end() == std::ranges::find_if_not(In, [&](Point const& p) {
               return p.SameDimension(In[0]);
             })) {  // NOTE: check whether all elements are identical
    return AllocDummyLeafNode<Slice, Leaf>(In);
  } else {  // NOTE: current dim d is same but other dims are not
    // WARN: this will break the rotate dimension mannar
    auto [new_box, new_dim] = split_rule_.SwitchDimension(In, d, bx);
    assert(IsMaxStretchSplit<SplitRule> || new_dim != d);
    return SerialBuildRecursive(In, Out, new_dim, new_box);
  }

  assert(std::ranges::all_of(In.begin(), splitIter, [&](Point& p) {
    return Num::Lt(p.pnt[split.second], split.first);
  }));
  assert(std::ranges::all_of(splitIter, In.end(), [&](Point& p) {
    return Num::Geq(p.pnt[split.second], split.first);
  }));

  Box lbox(bx), rbox(bx);
  lbox.second.pnt[d] = split.first;  //* loose
  rbox.first.pnt[d] = split.first;

  d = (d + 1) % BT::kDim;
  Node *L, *R;

  L = SerialBuildRecursive(In.cut(0, splitIter - In.begin()),
                           Out.cut(0, splitIter - In.begin()), d, lbox);
  R = SerialBuildRecursive(In.cut(splitIter - In.begin(), n),
                           Out.cut(splitIter - In.begin(), n), d, rbox);
  return AllocInteriorNode<Interior>(L, R, split, AugType());
}

template <typename Point, typename SplitRule, uint_fast8_t kBDO>
Node* KdTree<Point, SplitRule, kBDO>::BuildRecursive(Slice In, Slice Out,
                                                     DimsType dim,
                                                     Box const& bx) {
  assert(In.size() == 0 || BT::WithinBox(BT::GetBox(In), bx));

  // if (In.size()) {
  if (In.size() <= BT::kSerialBuildCutoff) {
    return SerialBuildRecursive(In, Out, dim, bx);
  }

  auto pivots = SplitterSeq::uninitialized(BT::kPivotNum + BT::kBucketNum + 1);
  auto boxs = BoxSeq::uninitialized(BT::kBucketNum);
  parlay::sequence<BallsType> sums;

  PickPivots(In, In.size(), pivots, dim, boxs, bx);
  BT::Partition(In, Out, In.size(), pivots, sums);

  // NOTE: if random sampling failed to split points, re-partitions using
  // serail approach
  auto tree_nodes = parlay::sequence<Node*>::uninitialized(BT::kBucketNum);
  auto nodes_map = BucketSeq::uninitialized(BT::kBucketNum);
  BucketType zeros = std::ranges::count(sums, 0), cnt = 0;

  if (zeros == BT::kBucketNum - 1) {  // NOTE: switch to seral
    // TODO: add parallelsim within this call
    // see parallel kth element
    return SerialBuildRecursive(In, Out, dim, bx);
  }

  // NOTE: alloc empty leaf beforehand to avoid spawn threads
  for (BucketType i = 0; i < BT::kBucketNum; ++i) {
    if (!sums[i]) {
      tree_nodes[i] = AllocEmptyLeafNode<Slice, Leaf>();
    } else {
      nodes_map[cnt++] = i;
    }
  }

  dim = (dim + BT::kBuildDepthOnce) % BT::kDim;

  parlay::parallel_for(
      0, BT::kBucketNum - zeros,
      [&](BucketType i) {
        size_t start = 0;
        for (BucketType j = 0; j < nodes_map[i]; ++j) {
          start += sums[j];
        }

        tree_nodes[nodes_map[i]] = BuildRecursive(
            Out.cut(start, start + sums[nodes_map[i]]),
            In.cut(start, start + sums[nodes_map[i]]), dim, boxs[nodes_map[i]]);
      },
      1);

  return BT::template BuildInnerTree<Interior>(1, pivots, tree_nodes);
}

template <typename Point, typename SplitRule, uint_fast8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::Build_(Slice A) {
  Points B = Points::uninitialized(A.size());
  this->tree_box_ = BT::GetBox(A);
  this->root_ = BuildRecursive(A, B.cut(0, A.size()), 0, this->tree_box_);
  assert(this->root_ != nullptr);
  return;
}

}  // namespace cpdd
