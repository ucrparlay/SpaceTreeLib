#ifndef PSTP_COVER_TREE_IMPL_COVER_BUILD_TREE_HPP_
#define PSTP_COVER_TREE_IMPL_COVER_BUILD_TREE_HPP_

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include <algorithm>

#include "pstp/cover_tree.h"
#include "pstp/dependence/tree_node.h"

namespace pstp {

// // TODO: maybe we don't need this function, it can be directly computed by
// value template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
//           uint_fast8_t kImbaRatio>
// void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::DivideRotate(
//     HyperPlaneSeq& pivots, DimsType dim, BucketType idx, BoxSeq& box_seq,
//     Box const& box) {
//   if (idx > BT::kPivotNum) {
//     // WARN: sometimes cut dimension can be -1, never use pivots[idx].first
//     ==
//     // -1 to check whether it is in bucket; instead, use idx > PIVOT_NUM
//     box_seq[idx - BT::kBucketNum] = box;
//     pivots[idx] = HyperPlane(0, idx - BT::kBucketNum);
//     return;
//   }
//
//   pivots[idx] = split_rule_.SplitSample(Slice(nullptr, nullptr), dim, box);
//
//   BoxCut box_cut(box, pivots[idx], true);
//   // dim = (dim + 1) % BT::kDim;
//   dim = split_rule_.NextDimension(dim);
//   DivideRotate(pivots, dim, 2 * idx, box_seq, box_cut.GetFirstBoxCut());
//   DivideRotate(pivots, dim, 2 * idx + 1, box_seq, box_cut.GetSecondBoxCut());
//
//   return;
// }
//
// template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
//           uint_fast8_t kImbaRatio>
// void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::SerialSplit(
//     Slice In, DimsType dim, DimsType idx, Box const& box,
//     parlay::sequence<BallsType>& sums) {
//   assert(dim <= BT::kDim);
//
//   if (dim == BT::kDim) {
//     sums[idx - kNodeRegions] = In.size();
//     return;
//   }
//
//   auto [split_iter, split] = split_rule_.SplitInput(In, dim, box);
//
//   SerialSplit(In.cut(0, split_iter - In.begin()), dim + 1, idx << 1, box,
//   sums); SerialSplit(In.cut(split_iter - In.begin(), In.size()), dim + 1, idx
//   << 1 | 1,
//               box, sums);
// }
//
// template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
//           uint_fast8_t kImbaRatio>
// Node* CoverTree<Point, SplitRule, kMD, kSkHeight,
//                 kImbaRatio>::SerialBuildRecursive(Slice In, Slice Out,
//                                                   Box const& box,
//                                                   bool checked_duplicate) {
//   assert(In.size() == 0 || BT::WithinBox(BT::GetBox(In), box));
//   size_t n = In.size();
//
//   if (n == 0) {
//     return AllocEmptyLeafNode<Slice, Leaf>();
//   }
//
//   if (n <= BT::kLeaveWrap) {
//     return AllocNormalLeafNode<Slice, Leaf>(In);
//   }
//
//   assert(kSplitterNum == BT::kDim);
//
//   DimsType dim = 0;
//   parlay::sequence<BallsType> sums(kNodeRegions, 0);
//   auto splitter = Interior::ComputeSplitter(box);
//   SerialSplit(In, dim, 1, box, sums);
//   assert(std::cmp_equal(std::accumulate(sums.begin(), sums.end(), 0), n));
//
//   if (std::ranges::count(sums, 0) == kNodeRegions - 1) {    // split fails
//     if (std::ranges::find_if_not(In, [&](Point const& p) {  // early return
//           return p.SameDimension(In[0]);
//         }) == In.end()) {
//       // WARN: Need to pass full range, since it needs to compute the size
//       return AllocDummyLeafNode<Slice, Leaf>(In);
//     } else {
//       return split_rule_.HandlingUndivide(*this, In, Out, box);
//     }
//   }
//
//   CoverNodeArr tree_nodes;
//   size_t start = 0;
//   for (DimsType i = 0; i < kNodeRegions; ++i) {
//     // NOTE: iterate through non-empty partitions, put them into the
//     // position identified by non_empty_node
//     tree_nodes[i] = SerialBuildRecursive(
//         In.cut(start, start + sums[i]), Out.cut(start, start + sums[i]),
//         Interior::GetBoxByRegionId(i, splitter, box), checked_duplicate);
//     start += sums[i];
//   }
//
//   return AllocInteriorNode<Interior>(tree_nodes, splitter, AugType());
// }
//
// template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
//           uint_fast8_t kImbaRatio>
// Node* CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::BuildRecursive(
//     Slice In, Slice Out, Box const& box) {
//   // TODO: may ensure the bucket is corresponding the the splitter
//   assert(In.size() == 0 || BT::WithinBox(BT::GetBox(In), box));
//   size_t n = In.size();
//
//   // if (In.size()) {
//   if (n <= BT::kSerialBuildCutoff) {
//     return SerialBuildRecursive(In, Out, box, false);
//   }
//
//   auto pivots =
//       HyperPlaneSeq::uninitialized(BT::kPivotNum + BT::kBucketNum + 1);
//   auto box_seq = BoxSeq::uninitialized(BT::kBucketNum);
//   parlay::sequence<BallsType> sums;
//
//   DivideRotate(pivots, 0, 1, box_seq, box);
//   BT::Partition(In, Out, In.size(), pivots, sums);
//
//   auto tree_nodes = parlay::sequence<Node*>::uninitialized(BT::kBucketNum);
//   auto nodes_map = BucketSeq::uninitialized(BT::kBucketNum);
//   BucketType zeros = std::ranges::count(sums, 0), cnt = 0;
//
//   if (zeros == BT::kBucketNum - 1) {  // NOTE: switch to seral
//     // TODO: add parallelsim within this call
//     // see parallel kth element
//     return SerialBuildRecursive(In, Out, box, false);
//   }
//
//   for (BucketType i = 0; i < BT::kBucketNum; ++i) {
//     if (sums[i] == 0) {
//       tree_nodes[i] = AllocEmptyLeafNode<Slice, Leaf>();
//     } else {
//       nodes_map[cnt++] = i;
//     }
//   }
//
//   parlay::parallel_for(
//       0, BT::kBucketNum - zeros,
//       [&](BucketType i) {
//         size_t start = 0;
//         for (BucketType j = 0; j < nodes_map[i]; ++j) {
//           start += sums[j];
//         }
//
//         tree_nodes[nodes_map[i]] = BuildRecursive(
//             Out.cut(start, start + sums[nodes_map[i]]),
//             In.cut(start, start + sums[nodes_map[i]]),
//             box_seq[nodes_map[i]]);
//       },
//       1);
//
//   return BT::template BuildInnerTree<Interior>(1, pivots, tree_nodes);
// }

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::Build_(Slice A) {
  Points B = Points::uninitialized(A.size());
  this->tree_box_ = BT::GetBox(A);
  // this->root_ = BuildRecursive(A, B.cut(0, A.size()), this->tree_box_);
  // this->root_ = SerialBuildRecursive(A, B.cut(0, A.size()), BT::kDim,
  //                                    this->tree_box_, false);
  assert(this->root_ != nullptr);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::Build_(
    Slice A, Box const& box) {
  assert(BT::WithinBox(BT::GetBox(A), box));

  Points B = Points::uninitialized(A.size());
  this->tree_box_ = box;
  // this->root_ = BuildRecursive(A, B.cut(0, A.size()), this->tree_box_);
  assert(this->root_ != nullptr);
  return;
}

template <typename Point, typename SplitRule, uint_fast8_t kSkHeight,
          uint_fast8_t kImbaRatio>
template <typename Range, typename... Args>
void CoverTree<Point, SplitRule, kSkHeight, kImbaRatio>::Build(Range&& In,
                                                               Args&&... args) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  Build_(A, std::forward<Args>(args)...);
}
}  // namespace pstp

#endif  // PSTP_COVER_TREE_IMPL_COVER_BUILD_TREE_HPP_
