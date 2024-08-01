#pragma once

#include <parlay/range.h>
#include <parlay/type_traits.h>
#include <parlay/slice.h>
#include "../kd_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point, typename SplitRule, uint8_t kBDO>
template<typename Range>
void KdTree<Point, SplitRule, kBDO>::Build(Range&& In, int DIM) {
    static_assert(parlay::is_random_access_range_v<Range>);
    static_assert(parlay::is_less_than_comparable_v<
                  parlay::range_reference_type_t<Range>>);
    static_assert(
        std::is_constructible_v<parlay::range_value_type_t<Range>,
                                parlay::range_reference_type_t<Range>>);

    Slice A = parlay::make_slice(In);
    Build_(A, DIM);
}

template<typename Point, typename SplitRule, uint8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::DivideRotate(Slice In, SplitterSeq& pivots,
                                                  DimsType dim, BucketType idx,
                                                  BucketType deep,
                                                  BucketType& bucket,
                                                  const DimsType DIM,
                                                  BoxSeq& boxs, const Box& bx) {
    if (deep > BT::kBuildDepthOnce) {
        // WARN: sometimes cut dimension can be -1
        //  never use pivots[idx].first to check whether it is in bucket;
        //  instead, use idx > PIVOT_NUM
        boxs[bucket] = bx;
        pivots[idx] = Splitter(-1, bucket++);
        return;
    }
    size_t n = In.size();
    uint_fast8_t d = split_rule_.FindCuttingDimension(bx, dim, DIM);
    assert(d < DIM);

#if __cplusplus <= 201703L
    std::nth_element(In.begin(), In.begin() + n / 2, In.end(),
                     [&](const Point& p1, const Point& p2) {
                         return Num::Lt(p1.pnt[d], p2.pnt[d]);
                     });
#else
    std::ranges::nth_element(In, In.begin() + n / 2,
                             [&](const Point& p1, const Point& p2) {
                                 return Num::Lt(p1.pnt[d], p2.pnt[d]);
                             });
#endif

    pivots[idx] = Splitter(In[n / 2].pnt[d], d);

    Box lbox(bx), rbox(bx);
    lbox.second.pnt[d] = pivots[idx].first;  // PERF: loose
    rbox.first.pnt[d] = pivots[idx].first;

    d = (d + 1) % DIM;
    DivideRotate(In.cut(0, n / 2), pivots, d, 2 * idx, deep + 1, bucket, DIM,
                 boxs, lbox);
    DivideRotate(In.cut(n / 2, n), pivots, d, 2 * idx + 1, deep + 1, bucket,
                 DIM, boxs, rbox);
    return;
}

// NOTE: starting at dimesion dim and pick pivots in a rotation manner
template<typename Point, typename SplitRule, uint8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::PickPivots(Slice In, const size_t& n,
                                                SplitterSeq& pivots,
                                                const DimsType dim,
                                                const DimsType DIM,
                                                BoxSeq& boxs, const Box& bx) {
    size_t size = std::min(n, static_cast<size_t>(32 * BT::kBucketNum));
    assert(size <= n);

    Points arr = Points::uninitialized(size);
    BT::SamplePoints(In, arr);

    // NOTE: pick pivots
    BucketType bucket = 0;
    DivideRotate(arr.cut(0, size), pivots, dim, 1, 1, bucket, DIM, boxs, bx);
    assert(bucket == BT::kBucketNum);
    return;
}

template<typename Point, typename SplitRule, uint8_t kBDO>
Node* KdTree<Point, SplitRule, kBDO>::SerialBuildRecursive(Slice In, Slice Out,
                                                           DimsType dim,
                                                           const DimsType DIM,
                                                           const Box& bx) {
    size_t n = In.size();

    if (n == 0) return AllocEmptyLeafNode<Slice, Leaf>();

    if (n <= BT::kLeaveWrap) return AllocNormalLeafNode<Slice, Leaf>(In);

    DimsType d = split_rule_.FindCuttingDimension(bx, dim, DIM);
    PointsIter splitIter = BT::SerialPartition(In, d);
    PointsIter diffEleIter;  // TODO: we can remove this iter

    Splitter split;

    if (splitIter <= In.begin() + n / 2) {  // NOTE: split is on left half
        split = Splitter(In[n / 2].pnt[d], d);
    } else if (splitIter != In.end()) {  // NOTE: split is on right half
#if __cplusplus <= 201703L
        auto minEleIter = std::min_element(
            splitIter, In.end(), [&](const Point& p1, const Point& p2) {
                return Num::Lt(p1.pnt[d], p2.pnt[d]);
            });
#else
        auto minEleIter = std::ranges::min_element(
            splitIter, In.end(), [&](const Point& p1, const Point& p2) {
                return Num::Lt(p1.pnt[d], p2.pnt[d]);
            });
#endif
        split = Splitter(minEleIter->pnt[d], d);
    } else if (In.end() ==
               (diffEleIter =
                    std::find_if_not(In.begin(), In.end(), [&](const Point& p) {
                        return p.sameDimension(In[0]);
                    }))) {  // NOTE: check whether all elements are identical
        return AllocDummyLeafNode<Slice, Leaf>(In.cut(0, 1));
    } else {  // NOTE: current dim d is same but other dims are not
        auto [new_box, new_dim] = split_rule_.SwitchDimension(In, d, DIM, bx);
        return SerialBuildRecursive(In, Out, new_dim, DIM, new_box);
    }

    assert(std::all_of(In.begin(), splitIter, [&](Point& p) {
        return Num::Lt(p.pnt[split.second], split.first);
    }));
    assert(std::all_of(splitIter, In.end(), [&](Point& p) {
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
    return AllocInteriorNode<Interior>(L, R, split, false);
}

template<typename Point, typename SplitRule, uint8_t kBDO>
Node* KdTree<Point, SplitRule, kBDO>::BuildRecursive(Slice In, Slice Out,
                                                     DimsType dim,
                                                     const DimsType DIM,
                                                     const Box& bx) {
    assert(In.size() == 0 || BT::WithinBox(BT::GetBox(In), bx));

    // if (In.size()) {
    if (In.size() <= BT::kSerialBuildCutoff) {
        return SerialBuildRecursive(In, Out, dim, DIM, bx);
    }

    auto pivots =
        SplitterSeq::uninitialized(BT::kPivotNum + BT::kBucketNum + 1);
    auto boxs = BoxSeq::uninitialized(BT::kBucketNum);
    parlay::sequence<BallsType> sums;

    PickPivots(In, In.size(), pivots, dim, DIM, boxs, bx);
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
            tree_nodes[i] = AllocEmptyLeafNode<Slice, Leaf>();
        } else {
            nodes_map[cnt++] = i;
        }
    }

    if (zeros == BT::kBucketNum - 1) {  // NOTE: switch to seral
        // TODO: add parallelsim within this call
        return SerialBuildRecursive(In, Out, dim, DIM, bx);
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

    return BT::template BuildInnerTree<Interior>(1, pivots, tree_nodes);
}

template<typename Point, typename SplitRule, uint8_t kBDO>
void KdTree<Point, SplitRule, kBDO>::Build_(Slice A, const DimsType DIM) {
    Points B = Points::uninitialized(A.size());
    this->tree_box_ = BT::GetBox(A);
    this->root_ =
        BuildRecursive(A, B.cut(0, A.size()), 0, DIM, this->tree_box_);
    assert(this->root_ != nullptr);
    return;
}

}  // namespace cpdd
