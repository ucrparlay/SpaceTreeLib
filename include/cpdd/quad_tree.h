#pragma once

#include <functional>
#include <utility>
#include "base_tree.h"

namespace cpdd {

template<typename Point, typename SplitRule, uint8_t kMD = 2, uint8_t kBDO = 6>
class QuadTree : private BaseTree<Point, kBDO> {
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

    struct QuadInteriorNode;

    using Leaf =
        LeafNode<Point, Slice, BT::kLeaveWrap, parlay::move_assign_tag>;
    using Interior = QuadInteriorNode;
    using Nodes = Interior::Nodes;

    // NOTE: expose basetree interface
    using BT::GetAveTreeHeight;
    using BT::GetBox;
    using BT::GetRoot;
    using BT::GetRootBox;
    using BT::Validate;

    // NOTE: functions
    template<typename Range>
    void Build(Range&& In, uint8_t DIM);

    void Build_(Slice In, const DimsType DIM);

    void DeleteTree() override;

    template<typename Range>
    void Flatten(Range&& Out);

    template<typename Range>
    void KNN(Node* T, const Point& q, const DimsType DIM,
             kBoundedQueue<Point, Range>& bq, const Box& bx,
             size_t& vis_node_num);

    size_t RangeCount(const Box& query_box);

    size_t RangeCount(const Circle& cl);

    template<typename Range>
    size_t RangeQuery(const Box& query_box, Range&& Out);

 private:
    void SerialSplit(Slice In, DimsType dim, DimsType DIM, DimsType idx,
                     const Box& box, const Splitter& split,
                     parlay::sequence<BallsType>& sums, BoxSeq& box_seq);

    void DivideRotate(HyperPlaneSeq& pivots, DimsType dim, BucketType idx,
                      DimsType deep, BucketType& bucket, const DimsType DIM,
                      BoxSeq& box_seq, const Box& box);

    void PickPivots(Slice In, const size_t& n, HyperPlaneSeq& pivots,
                    const DimsType dim, const DimsType DIM, BoxSeq& boxs,
                    const Box& bx);

    Node* BuildRecursive(Slice In, Slice Out, DimsType dim, const DimsType DIM,
                         const Box& bx);

    Node* SerialBuildRecursive(Slice In, Slice Out, DimsType dim,
                               const DimsType DIM, const Box& bx);

    static Node* QuadBuildInnerTree(BucketType idx, const HyperPlaneSeq& pivots,
                                    const parlay::sequence<Node*>& tree_nodes);

    SplitRule split_rule_;
};

}  // namespace cpdd

#include "quad_tree_impl/quad_build_tree.hpp"
#include "quad_tree_impl/quad_inter_node.hpp"
#include "quad_tree_impl/quad_override.hpp"