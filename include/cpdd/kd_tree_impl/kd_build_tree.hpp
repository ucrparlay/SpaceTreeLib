#pragma once

#include <parlay/range.h>
#include <parlay/type_traits.h>
#include <parlay/slice.h>
#include "../kd_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point, typename SplitRule>
template<typename Range>
void KdTree<Point, SplitRule>::Build(Range&& In, int DIM) {
    static_assert(parlay::is_random_access_range_v<Range>);
    static_assert(parlay::is_less_than_comparable_v<
                  parlay::range_reference_type_t<Range>>);
    static_assert(
        std::is_constructible_v<parlay::range_value_type_t<Range>,
                                parlay::range_reference_type_t<Range>>);

    Slice A = parlay::make_slice(In);
    Build_(A, DIM);
}

template<typename Point, typename SplitRule>
Node* KdTree<Point, SplitRule>::SerialBuildRecursive(Slice In, Slice Out,
                                                     DimsType dim,
                                                     const DimsType DIM,
                                                     const Box& bx) {
    size_t n = In.size();

    if (n == 0) return AllocLeafNode<Point, Slice, AllocEmptyLeafTag>();

    if (n <= BT::kLeaveWrap)
        return AllocLeafNode<Point, Slice, AllocNormalLeafTag>(In);

    // DimsType d =
    //     (split_rule_ == kMaxStretchDim ? pick_max_stretch_dim(bx, DIM) :
    //     dim);
    DimsType d = split_rule_.FindCuttingDimension(bx, dim, DIM);
    PointsIter splitIter = BT::SerialPartition(In, d);
    PointsIter diffEleIter;

    Splitter split;

    if (splitIter <= In.begin() + n / 2) {  // NOTE: split is on left half
        split = splitter(In[n / 2].pnt[d], d);
    } else if (splitIter != In.end()) {  // NOTE: split is on right half
        auto minEleIter = std::ranges::min_element(
            splitIter, In.end(), [&](const Point& p1, const Point& p2) {
                return Num::Lt(p1.pnt[d], p2.pnt[d]);
            });
        split = splitter(minEleIter->pnt[d], d);
    } else if (In.end() ==
               (diffEleIter = std::ranges::find_if_not(In, [&](const Point& p) {
                    return p.sameDimension(In[0]);
                }))) {  // NOTE: check whether all elements are identical
        return AllocLeafNode<Point, Slice, AllocDummyLeafTag>(In);
    } else {  // NOTE: current dim d is same but other dims are not
        auto [new_box, new_dim] = split_rule_.SwitchDimension();
        SerialBuildRecursive(In, Out, new_dim, DIM, new_box);
    }

    assert(std::ranges::all_of(In.begin(), splitIter,
                               [&](Point& p) {
                                   return Num::Lt(p.pnt[split.second],
                                                  split.first);
                               }) &&
           std::ranges::all_of(splitIter, In.end(), [&](Point& p) {
               return Num::Geq(p.pnt[split.second], split.first);
           }));

    Box lbox(bx), rbox(bx);
    lbox.second.pnt[d] = split.first;  //* loose
    rbox.first.pnt[d] = split.first;

    d = (d + 1) % DIM;
    Node *L, *R;

    L = SerialBuildRecursive(In.cut(0, splitIter - In.begin()),
                             Out.cut(0, splitIter - In.begin()), d, DIM, lbox);
    R = SerialBuildRecursive(In.cut(splitIter - In.begin(), n),
                             Out.cut(splitIter - In.begin(), n), d, DIM, rbox);
    return alloc_interior_node(L, R, split);
}

template<typename Point, typename SplitRule>
Node* KdTree<Point, SplitRule>::BuildRecursive(Slice In, Slice Out,
                                               DimsType dim, const DimsType DIM,
                                               const Box& bx) {
    assert(In.size() == 0 || within_box(get_box(In), bx));

    // if ( In.size() ) {
    if (In.size() <= BT::kSerialBuildCutoff) {
        return SerialBuildRecursive(In, Out, dim, DIM, bx);
    }

    //* parallel partitons
    auto pivots =
        SplitterSeq::uninitialized(BT::kPivotNum + BT::kBucketNum + 1);
    auto boxs = BoxSeq::uninitialized(BT::kBucketNum);
    parlay::sequence<BallsType> sums;

    BT::PickPivots(In, In.size(), pivots, dim, DIM, boxs, bx);
    BT::Partition(In, Out, In.size(), pivots, sums);

    // NOTE: if random sampling failed to split points, re-partitions using
    // serail approach
    auto tree_nodes = parlay::sequence<Node*>::uninitialized(BT::kBucketNum);
    auto nodes_map =
        parlay::sequence<BucketType>::uninitialized(BT::kBucketNum);
    BucketType zeros = 0, cnt = 0;
    for (BucketType i = 0; i < BT::kBucketNum; ++i) {
        if (!sums[i]) {
            ++zeros;
            tree_nodes[i] = AllocLeafNode<Point, Slice, AllocEmptyLeafTag>();
        } else {
            nodes_map[cnt++] = i;
        }
    }

    if (zeros == BT::kBucketNum - 1) {  // NOTE: switch to seral
        // TODO: add parallelsim within this call
        return SeialBuildRecursive(In, Out, dim, DIM, bx);
    }

    dim = (dim + BT::kBuildDepthOnce) % DIM;

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
                               DIM, boxs[nodes_map[i]]);
        },
        1);
    return BuildInnerTree(1, pivots, tree_nodes);
}

template<typename Point, typename SplitRule>
void KdTree<Point, SplitRule>::Build_(Slice A, const DimsType DIM) {
    Points B = Points::uninitialized(A.size());
    this->tree_box_ = BT::GetBox(A);
    // this->root = build_recursive(A, B.cut(0, A.size()), 0, DIM, this->bbox);
    assert(this->root_ != nullptr);
    return;
}

}  // namespace cpdd
