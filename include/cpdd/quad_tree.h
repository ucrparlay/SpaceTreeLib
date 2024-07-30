#pragma once

#include <functional>
#include <utility>
#include "base_tree.h"

namespace cpdd {

template<typename Point, typename SplitRule, uint8_t kBDO = 6>
class QuadTree : private BaseTree<Point, kBDO> {
 public:
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

    using Splitter = std::pair<Coord, DimsType>;
    using SplitterSeq = parlay::sequence<Splitter>;
    using SplitType = std::pair<Splitter, Splitter>;
    using AugType = bool;

    struct QuadInteriorNode;

    using Leaf =
        LeafNode<Point, Slice, BT::kLeaveWrap, parlay::move_assign_tag>;
    using Interior = QuadInteriorNode;

    // NOTE: expose basetree interface
    using BT::GetAveTreeHeight;
    using BT::GetBox;
    using BT::GetRoot;
    using BT::GetRootBox;
    using BT::Validate;

    // NOTE: functions
    template<typename Range>
    void Build(Range&& In, int DIM);

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
    void DivideRotate(Slice In, SplitterSeq& pivots, DimsType dim,
                      BucketType idx, BucketType deep, BucketType& bucket,
                      const DimsType DIM, BoxSeq& boxs, const Box& bx);

    void PickPivots(Slice In, const size_t& n, SplitterSeq& pivots,
                    const DimsType dim, const DimsType DIM, BoxSeq& boxs,
                    const Box& bx);

    Node* BuildRecursive(Slice In, Slice Out, DimsType dim, const DimsType DIM,
                         const Box& bx);

    Node* SerialBuildRecursive(Slice In, Slice Out, DimsType dim,
                               const DimsType DIM, const Box& bx);

    SplitRule split_rule_;
};

}  // namespace cpdd

#include "quad_tree_impl/quad_build_tree.hpp"
#include "quad_tree_impl/quad_inter_node.hpp"
#include "quad_tree_impl/quad_override.hpp"
