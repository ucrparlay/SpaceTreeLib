#ifndef PSPT_R_TREE_IMPL_ORTH_BUILD_TREE_HPP_
#define PSPT_R_TREE_IMPL_ORTH_BUILD_TREE_HPP_

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include "../r_tree.h"
#include "pspt/dependence/tree_node.h"

namespace pspt {
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range>
void RTree<Point, SplitRule, kSkHeight, kImbaRatio>::Build(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  Build_(A);
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void RTree<Point, SplitRule, kSkHeight, kImbaRatio>::Build_(Slice A) {
  Points B = Points::uninitialized(A.size());
  this->tree_box_ = BT::GetBox(A.cut(0, A.size()));
  this->root_ = BuildRecursive(A, B.cut(0, A.size()), 0, this->tree_box_);
  assert(this->root_ != nullptr);
  assert(BT::SameBox(this->tree_box_,
                     BT::template GetSplit<Leaf, Interior>(this->root_)));
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void RTree<Point, SplitRule, kSkHeight, kImbaRatio>::DivideRotate(
    Slice In, HyperPlaneSeq& pivots, DimsType dim, BucketType idx,
    BoxSeq& box_seq, Box const& bx) {
  if (idx > BT::kPivotNum) {
    // WARN: sometimes cut dimension can be -1
    //  never use pivots[idx].first to check whether it is in bucket;
    //  instead, use idx > PIVOT_NUM
    box_seq[idx - BT::kBucketNum] = bx;
    pivots[idx] = HyperPlane(0, idx - BT::kBucketNum);
    return;
  }
  size_t n = In.size();
  DimsType d = split_rule_.FindCuttingDimension(bx, dim);
  assert(d < BT::kDim);

  std::ranges::nth_element(In, In.begin() + n / 2,
                           [&](Point const& p1, Point const& p2) {
                             return Num::Lt(p1.pnt[d], p2.pnt[d]);
                           });

  pivots[idx] = HyperPlane(In[n / 2].pnt[d], d);

  BoxCut box_cut(bx, pivots[idx], true);

  d = (d + 1) % BT::kDim;
  DivideRotate(In.cut(0, n / 2), pivots, d, 2 * idx, box_seq,
               box_cut.GetFirstBoxCut());
  DivideRotate(In.cut(n / 2, n), pivots, d, 2 * idx + 1, box_seq,
               box_cut.GetSecondBoxCut());
  return;
}

// NOTE: starting at dimesion dim and pick pivots in a rotation manner
template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void RTree<Point, SplitRule, kSkHeight, kImbaRatio>::PickPivots(
    Slice In, size_t const& n, HyperPlaneSeq& pivots, DimsType const dim,
    BoxSeq& box_seq, Box const& bx) {
  size_t size = std::min(n, static_cast<size_t>(32 * BT::kBucketNum));
  assert(size <= n);

  Points arr = Points::uninitialized(size);
  BT::SamplePoints(In, arr);

  // NOTE: pick pivots
  DivideRotate(arr.cut(0, size), pivots, dim, 1, box_seq, bx);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
Node* RTree<Point, SplitRule, kSkHeight, kImbaRatio>::SerialBuildRecursive(
    Slice In, Slice Out, DimsType dim, Box const& bx) {
  size_t n = In.size();

  if (n == 0) {
    return AllocEmptyLeafNode<Slice, Leaf>(BT::GetEmptyBox());
  }

  if (n <= BT::kLeaveWrap) {
    return AllocNormalLeafNode<Slice, Leaf>(In, BT::GetBox(In));
  }

  DimsType d = split_rule_.FindCuttingDimension(bx, dim);
  PointsIter split_iter = BT::SerialPartition(In, d);

  HyperPlane split;

  if (split_iter <= In.begin() + n / 2) {  // split is on left half
    split = HyperPlane(In[n / 2].pnt[d], d);
  } else if (split_iter != In.end()) {  // split is on right half
    auto min_elem_iter = std::ranges::min_element(
        split_iter, In.end(), [&](Point const& p1, Point const& p2) {
          return Num::Lt(p1.pnt[d], p2.pnt[d]);
        });
    split = HyperPlane(min_elem_iter->pnt[d], d);
  } else if (In.end() == std::ranges::find_if_not(In, [&](Point const& p) {
               return p.SameDimension(In[0]);
             })) {  //  check whether all elements are identical
    return AllocDummyLeafNode<Slice, Leaf>(In);
  } else {  //  current dim d is same but other dims are not
    // WARN: this will break the rotate dimension mannar
    // auto [new_box, new_dim] = split_rule_.SwitchDimension(In, d, bx);
    // BUG: handling rtree switch
    auto new_dim = 0;
    auto new_box = BT::GetBox(In);
    assert(IsMaxStretchDim<SplitRule> || new_dim != d);
    return SerialBuildRecursive(In, Out, new_dim, new_box);
  }

  assert(std::ranges::all_of(In.begin(), split_iter, [&](Point& p) {
    return Num::Lt(p.pnt[split.second], split.first);
  }));
  assert(std::ranges::all_of(split_iter, In.end(), [&](Point& p) {
    return Num::Geq(p.pnt[split.second], split.first);
  }));

  BoxCut box_cut(bx, split, true);

  d = (d + 1) % BT::kDim;
  Node *L, *R;

  L = SerialBuildRecursive(In.cut(0, split_iter - In.begin()),
                           Out.cut(0, split_iter - In.begin()), d,
                           box_cut.GetFirstBoxCut());
  R = SerialBuildRecursive(In.cut(split_iter - In.begin(), n),
                           Out.cut(split_iter - In.begin(), n), d,
                           box_cut.GetSecondBoxCut());

  return AllocInteriorNode<Interior>(
      L, R,
      BT::GetBox(BT::template GetSplit<Leaf, Interior>(L),
                 BT::template GetSplit<Leaf, Interior>(R)),
      AugType());
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
Node* RTree<Point, SplitRule, kSkHeight, kImbaRatio>::BuildRecursive(
    Slice In, Slice Out, DimsType dim, Box const& bx) {
  assert(In.size() == 0 || BT::WithinBox(BT::GetBox(In), bx));

  // if (In.size()) {
  if (In.size() <= BT::kSerialBuildCutoff) {
    return SerialBuildRecursive(In, Out, dim, bx);
  }

  auto pivots =
      HyperPlaneSeq::uninitialized(BT::kPivotNum + BT::kBucketNum + 1);
  auto box_seq = BoxSeq::uninitialized(BT::kBucketNum);
  parlay::sequence<BallsType> sums;

  PickPivots(In, In.size(), pivots, dim, box_seq, bx);
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
      tree_nodes[i] = AllocEmptyLeafNode<Slice, Leaf>(BT::GetEmptyBox());
    } else {
      nodes_map[cnt++] = i;
    }
  }

  dim = (dim + BT::kBuildDepthOnce) % BT::kDim;

  parlay::parallel_for(
      0, BT::kBucketNum - zeros,
      [&](BucketType i) {
        size_t start =
            std::accumulate(sums.begin(), sums.begin() + nodes_map[i], 0);

        tree_nodes[nodes_map[i]] =
            BuildRecursive(Out.cut(start, start + sums[nodes_map[i]]),
                           In.cut(start, start + sums[nodes_map[i]]), dim,
                           box_seq[nodes_map[i]]);
      },
      1);

  return BT::template BuildInnerTree<Leaf, Interior>(1, pivots, tree_nodes);
}

}  // namespace pspt

#endif  // PSPT_R_TREE_IMPL_ORTH_BUILD_TREE_HPP_
