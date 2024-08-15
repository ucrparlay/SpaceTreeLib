#pragma once

#include <functional>
#include <utility>
#include "base_tree.h"

namespace cpdd {

template<typename Point, typename SplitRule, uint8_t kMD = 2, uint8_t kBDO = 6>
class OrthTree : private BaseTree<Point, kBDO> {
 public:
    static constexpr size_t kSplitterNum = kMD;
    static constexpr size_t kNodeRegions = 1 << kMD;

    using BT = BaseTree<Point, kBDO>;

    using BucketType = BT::BucketType;
    using BallsType = BT::BallsType;
    using DimsType = BT::DimsType;
    using Coord = typename Point::Coord;
    using Coords = typename Point::Coords;
    using Num = Num_Comparator<Coord>;
    using Slice = BT::Slice;
    using Points = BT::Points;
    using PointsIter = BT::PointsIter;
    using Box = BT::Box;
    using BoxSeq = BT::BoxSeq;
    using Circle = BT::Circle;

    using HyperPlane = BT::HyperPlane;
    using HyperPlaneSeq = BT::HyperPlaneSeq;
    using Splitter = std::array<HyperPlane, kSplitterNum>;
    using SplitterSeq = parlay::sequence<Splitter>;
    using AugType = bool;

    struct KdInteriorNode;
    struct OrthInteriorNode;

    using SplitRuleType = SplitRule;
    using Leaf =
        LeafNode<Point, Slice, BT::kLeaveWrap, parlay::move_assign_tag>;
    using Interior = OrthInteriorNode;
    using Nodes = Interior::Nodes;

    // NOTE: expose basetree interface
    using BT::Expand2Binary;
    using BT::GetAveTreeHeight;
    using BT::GetBox;
    using BT::GetMaxTreeDepth;
    using BT::GetRoot;
    using BT::GetRootBox;
    using BT::SetRoot;
    using BT::Validate;

    // NOTE: functions
    template<typename Range>
    void Build(Range&& In);

    void DeleteTree() override;

    template<typename Range>
    void Flatten(Range&& Out);

    template<typename Range>
    auto KNN(Node* T, const Point& q, kBoundedQueue<Point, Range>& bq);

    size_t RangeCount(const Box& query_box);

    size_t RangeCount(const Circle& cl);

    template<typename Range>
    size_t RangeQuery(const Box& query_box, Range&& Out);

 private:
    void Build_(Slice In);

    void SerialSplit(Slice In, DimsType dim, DimsType idx, const Box& box,
                     const Splitter& split, parlay::sequence<BallsType>& sums,
                     BoxSeq& box_seq);

    void DivideRotate(HyperPlaneSeq& pivots, DimsType dim, BucketType idx,
                      DimsType deep, BucketType& bucket, BoxSeq& box_seq,
                      const Box& box);

    void PickPivots(Slice In, const size_t& n, HyperPlaneSeq& pivots,
                    const DimsType dim, BoxSeq& boxs, const Box& bx);

    Node* BuildRecursive(Slice In, Slice Out, const Box& bx);

    Node* SerialBuildRecursive(Slice In, Slice Out, const Box& bx,
                               bool checked_duplicate);

    static Node* QuadBuildInnerTree(BucketType idx, const HyperPlaneSeq& pivots,
                                    const parlay::sequence<Node*>& tree_nodes);

    SplitRule split_rule_;
    size_t alloc_dummy_num_ = 0;
    size_t alloc_empty_num_ = 0;
    size_t alloc_normal_num_ = 0;
    size_t alloc_interior_num_ = 0;
};

}  // namespace cpdd

#include "orth_tree_impl/orth_build_tree.hpp"
#include "orth_tree_impl/orth_inter_node.hpp"
#include "orth_tree_impl/orth_override.hpp"
