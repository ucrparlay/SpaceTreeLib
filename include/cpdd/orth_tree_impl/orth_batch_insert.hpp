#pragma once

#include "../orth_tree.h"
#include "parlay/primitives.h"

namespace cpdd {

template<typename Point, typename SplitRule, uint_fast8_t kMD,
         uint_fast8_t kBDO>
template<typename Range>
void OrthTree<Point, SplitRule, kMD, kBDO>::BatchInsert(Range&& In) {
    static_assert(parlay::is_random_access_range_v<Range>);
    static_assert(parlay::is_less_than_comparable_v<
                  parlay::range_reference_type_t<Range>>);
    static_assert(
        std::is_constructible_v<parlay::range_value_type_t<Range>,
                                parlay::range_reference_type_t<Range>>);
    static_assert(BT::kBuildDepthOnce % kMD == 0);
    assert(kMD == BT::kDim);
    // TODO: handling the case that insert box is no in the tree box
    assert(BT::WithinBox(BT::GetBox(In), this->tree_box_));

    Slice A = parlay::make_slice(In);
    BatchInsert_(A);
}

template<typename Point, typename SplitRule, uint_fast8_t kMD,
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

template<typename Point, typename SplitRule, uint_fast8_t kMD,
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

    PointsIter split_iter = std::ranges::partition(In, [&](const Point& p) {
                                return Num::Lt(p.pnt[dim], mid);
                            }).begin();

    SerialSplitSkeleton(T, In.cut(0, split_iter - In.begin()), dim + 1, 2 * idx,
                        sums);
    SerialSplitSkeleton(T, In.cut(split_iter - In.begin(), In.size()), dim + 1,
                        2 * idx + 1, sums);
    return;
}

template<typename Point, typename SplitRule, uint_fast8_t kMD,
         uint_fast8_t kBDO>
Node* OrthTree<Point, SplitRule, kMD, kBDO>::UpdateInnerTreeByTag(
    BucketType idx, const NodeTagSeq& tags, parlay::sequence<Node*>& tree_nodes,
    BucketType& p) {
    if (tags[idx].second < BT::kBucketNum) {
        return tree_nodes[p++];
    }

    assert(tags[idx].second == BT::kBucketNum);
    assert(tags[idx].first != NULL);
    Nodes new_nodes;
    for (BucketType i = 0; i < kNodeRegions; ++i) {
        new_nodes[i] =
            UpdateInnerTreeByTag(idx * kNodeRegions + i, tags, tree_nodes, p);
    }
    BT::template UpdateInterior<Interior>(tags[idx].first, new_nodes);
    return tags[idx].first;
}

template<typename Point, typename SplitRule, uint_fast8_t kMD,
         uint_fast8_t kBDO>
Node* OrthTree<Point, SplitRule, kMD, kBDO>::BatchInsertRecursive(Node* T,
                                                                  Slice In,
                                                                  Slice Out) {
    size_t n = In.size();

    if (n == 0) return T;

    if (T->is_leaf) {
        Leaf* TL = static_cast<Leaf*>(T);
        if ((!TL->is_dummy && n + T->size <= BT::kLeaveWrap) ||
            (TL->is_dummy && parlay::all_of(In, [&](const Point& p) {
                 return p == TL->pts[0];
             }))) {
            return BT::template InsertPoints2Leaf<Leaf>(T, In);
        } else {
            return BT::template RebuildWithInsert<Leaf, Interior>(T, In);
        }
    }

    // if (n) {
    if (n <= BT::kSerialBuildCutoff) {
        parlay::sequence<BallsType> sums(kNodeRegions, 0);
        SerialSplitSkeleton(T, In, 0, 1, sums);
        assert(std::accumulate(sums.begin(), sums.end(), 0) == n);

        Nodes new_nodes;
        size_t start = 0;
        for (DimsType i = 0; i < kNodeRegions; ++i) {
            new_nodes[i] =
                BatchInsertRecursive(static_cast<Interior*>(T)->tree_nodes[i],
                                     In.cut(start, start + sums[i]),
                                     Out.cut(start, start + sums[i]));
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

    // IT.TagInbalanceNode();
    // assert(IT.tags_num > 0 && IT.tags_num <= BT::kBucketNum);
    auto tree_nodes = parlay::sequence<Node*>::uninitialized(IT.tags_num);

    // TODO: if none points are sieved into bucket, skip that nodes
    parlay::parallel_for(
        0, IT.tags_num,
        [&](size_t i) {
            size_t s = 0;
            for (int j = 0; j < i; j++) {
                s += IT.sums[j];
            }

            tree_nodes[i] = BatchInsertRecursive(IT.tags[IT.rev_tag[i]].first,
                                                 Out.cut(s, s + IT.sums[i]),
                                                 In.cut(s, s + IT.sums[i]));
        },
        1);

    BucketType beatles = 0;
    // TODO: rewrite to use the one in inner tree
    return UpdateInnerTreeByTag(1, IT.tags, tree_nodes, beatles);
}
}  // namespace cpdd
