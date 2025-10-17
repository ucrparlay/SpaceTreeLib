#ifndef PSI_POINTER_VIEW_KD_TREE_IMPL_KD_BUILD_TREE_HPP_
#define PSI_POINTER_VIEW_KD_TREE_IMPL_KD_BUILD_TREE_HPP_

#include <parlay/range.h>
#include <parlay/slice.h>
#include <parlay/type_traits.h>

#include "dependence/tree_node.h"
#include "pointer_view/kd_tree.h"

namespace psi {
namespace pointer_view {

template <typename TypeTrait>
template <typename Range>
void KdTree<TypeTrait>::Build(Range&& In) {
  static_assert(parlay::is_random_access_range_v<Range>);
  static_assert(
      parlay::is_less_than_comparable_v<parlay::range_reference_type_t<Range>>);
  static_assert(std::is_constructible_v<parlay::range_value_type_t<Range>,
                                        parlay::range_reference_type_t<Range>>);

  Slice A = parlay::make_slice(In);
  Build_(A);
}

template <typename TypeTrait>
void KdTree<TypeTrait>::DivideRotate(Slice In, SplitterSeq& pivots,
                                     DimsType dim, BucketType idx,
                                     BoxSeq& box_seq, Box const& box) {
  if (idx > BT::kPivotNum) {
    box_seq[idx - BT::kBucketNum] = box;
    pivots[idx] = Splitter(0, idx - BT::kBucketNum);
    return;
  }
  size_t n = In.size();
  DimsType cutting_dim = split_rule_.FindCuttingDimension(box, dim);
  assert(cutting_dim < BT::kDim);

  pivots[idx] = split_rule_.SplitSample(In, cutting_dim, box);

  BoxCut box_cut(box, pivots[idx], true);

  cutting_dim = split_rule_.NextDimension(cutting_dim);
  DivideRotate(In.cut(0, n / 2), pivots, cutting_dim, 2 * idx, box_seq,
               box_cut.GetFirstBoxCut());
  DivideRotate(In.cut(n / 2, n), pivots, cutting_dim, 2 * idx + 1, box_seq,
               box_cut.GetSecondBoxCut());
  return;
}

template <typename TypeTrait>
void KdTree<TypeTrait>::PickPivots(Slice In, size_t const& n,
                                   SplitterSeq& pivots, DimsType const dim,
                                   BoxSeq& box_seq, Box const& box) {
  size_t size = std::min(n, static_cast<size_t>(32 * BT::kBucketNum));
  assert(size <= n);

  Points arr = Points::uninitialized(size);
  BT::SamplePoints(In, arr);

  DivideRotate(arr.cut(0, size), pivots, dim, 1, box_seq, box);
  return;
}

template <typename TypeTrait>
auto KdTree<TypeTrait>::SerialBuildRecursive(Slice In, Slice Out, DimsType dim,
                                             Box const& box) -> Node* {
  size_t n = In.size();

  if (n == 0) return AllocEmptyLeafNode<Slice, Leaf>();

  if (n <= BT::kLeaveWrap) return AllocNormalLeafNode<Slice, Leaf>(In);

  DimsType d = split_rule_.FindCuttingDimension(box, dim);
  auto [split_iter, split] = split_rule_.SplitInput(In, d, box);

  if (!split.has_value()) {
    if (In.end() == std::ranges::find_if_not(In, [&](Point const& p) {
          return p.SameDimension(In[0]);
        })) {
      if constexpr (IsAugPoint<Point>) {
        if constexpr (Point::IsNonTrivialAugmentation()) {
          return AllocFixSizeLeafNode<Slice, Leaf>(
              In, std::max(In.size(), static_cast<size_t>(BT::kLeaveWrap)));
        } else {
          return AllocDummyLeafNode<Slice, Leaf>(In);
        }
      } else {
        return AllocDummyLeafNode<Slice, Leaf>(In);
      }
    } else {
      return split_rule_.HandlingUndivide(*this, In, Out, box, dim);
    }
  }

  assert(std::ranges::all_of(In.begin(), split_iter, [&](Point& p) {
    return Num::Lt(p.pnt[split.value().second], split.value().first);
  }));
  assert(std::ranges::all_of(split_iter, In.end(), [&](Point& p) {
    return Num::Geq(p.pnt[split.value().second], split.value().first);
  }));

  BoxCut box_cut(box, split.value(), true);

  d = split_rule_.NextDimension(d);
  Node *L, *R;

  L = SerialBuildRecursive(In.cut(0, split_iter - In.begin()),
                           Out.cut(0, split_iter - In.begin()), d,
                           box_cut.GetFirstBoxCut());
  R = SerialBuildRecursive(In.cut(split_iter - In.begin(), n),
                           Out.cut(split_iter - In.begin(), n), d,
                           box_cut.GetSecondBoxCut());
  return AllocInteriorNode<Interior>(L, R, split.value());
}

template <typename TypeTrait>
auto KdTree<TypeTrait>::BuildRecursive(Slice In, Slice Out, DimsType dim,
                                       Box const& box) -> Node* {
  assert(In.size() == 0 || Geo::WithinBox(Geo::GetBox(In), box));

  if (In.size() <= BT::kSerialBuildCutoff) {
    return SerialBuildRecursive(In, Out, dim, box);
  }

  auto pivots = SplitterSeq::uninitialized(BT::kPivotNum + BT::kBucketNum + 1);
  auto box_seq = BoxSeq::uninitialized(BT::kBucketNum);
  parlay::sequence<BallsType> sums;

  PickPivots(In, In.size(), pivots, dim, box_seq, box);
  BT::Partition(In, Out, In.size(), pivots, sums);

  auto tree_nodes = parlay::sequence<Node*>::uninitialized(BT::kBucketNum);
  auto nodes_map = BucketSeq::uninitialized(BT::kBucketNum);
  BucketType zeros = std::ranges::count(sums, 0), cnt = 0;

  if (zeros == BT::kBucketNum - 1) {
    return SerialBuildRecursive(In, Out, dim, box);
  }

  for (BucketType i = 0; i < BT::kBucketNum; ++i) {
    if (!sums[i]) {
      tree_nodes[i] = AllocEmptyLeafNode<Slice, Leaf>();
    } else {
      nodes_map[cnt++] = i;
    }
  }

  for (BucketType i = 0; i < BT::kBuildDepthOnce; ++i) {
    dim = split_rule_.NextDimension(dim);
  }

  parlay::parallel_for(
      0, BT::kBucketNum - zeros,
      [&](BucketType i) {
        size_t start = 0;
        for (BucketType j = 0; j < nodes_map[i]; ++j) {
          start += sums[j];
        }

        tree_nodes[nodes_map[i]] =
            BuildRecursive(Out.cut(start, start + sums[nodes_map[i]]),
                           In.cut(start, start + sums[nodes_map[i]]), dim,
                           box_seq[nodes_map[i]]);
      },
      1);

  return BT::template BuildInnerTree<Leaf, Interior>(1, pivots, tree_nodes);
}

template <typename TypeTrait>
void KdTree<TypeTrait>::Build_(Slice A) {
  Points B = Points::uninitialized(A.size());
  this->tree_box_ = Geo::GetBox(A);
  this->root_ = BuildRecursive(A, B.cut(0, A.size()), 0, this->tree_box_);
  assert(this->root_ != nullptr);
  return;
}

}  // namespace pointer_view
}  // namespace psi

#endif  // PSI_POINTER_VIEW_KD_TREE_IMPL_KD_BUILD_TREE_HPP_
