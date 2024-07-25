#pragma once

#include <parlay/range.h>
#include <parlay/type_traits.h>
#include <parlay/slice.h>
#include "../kd_tree.h"
#include "cpdd/dependence/tree_node.h"

namespace cpdd {
template<typename Point>
template<typename R>
void KdTree<Point>::Build(R&& In, int DIM) {
    static_assert(parlay::is_random_access_range_v<R>);
    static_assert(
        parlay::is_less_than_comparable_v<parlay::range_reference_type_t<R>>);
    static_assert(std::is_constructible_v<parlay::range_value_type_t<R>,
                                          parlay::range_reference_type_t<R>>);
    Slice A = parlay::make_slice(In);
    Build_(A, DIM);
}

template<typename Point>
Node* KdTree<Point>::BuildRecursive(Slice In, Slice Out, DimsType dim,
                                    const DimsType DIM, const Box& bx) {
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

template<typename Point>
void KdTree<Point>::Build_(Slice A, const DimsType DIM) {
    Points B = Points::uninitialized(A.size());
    this->tree_box_ = BT::GetBox(A);
    // this->root = build_recursive(A, B.cut(0, A.size()), 0, DIM, this->bbox);
    assert(this->root_ != nullptr);
    return;
}

}  // namespace cpdd
